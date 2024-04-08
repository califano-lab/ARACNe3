#include <iostream>
#include <memory>
#include <numeric>

#include "analysis.hpp"
#include "apmi_nullmodel.hpp"
#include "stopwatch.hpp"
#include "subnet_operations.hpp"

int ARACNe3Analysis(
    const ARACNe3IOHandler &io, const std::string &version,
    const std::string &runid, const uint32_t seed, const uint8_t threads,
    const bool verbose, const float alpha, const float subsamp_pct,
    const uint32_t n_nulls, const std::string &method, const bool prune_alpha,
    const bool prune_MaxEnt, const bool save_subnets,
    const std::string &subnets_log_dir, uint16_t n_subnets,
    const bool adaptive, const uint16_t min_subnets, const uint16_t max_subnets,
    const uint16_t min_regulon_occpuancy, const bool consolidate_mode,
    const bool skip_consolidate, Logger *const aracne3_logger,
    const std::string &subnets_dir, const std::string &cached_dir) {

  // ---- Quick macros ----

  const std::string M = " -A3- ";
  auto qlog = [&](const std::string &cur_msg) {
    std::cout << M << cur_msg << std::endl;
    if (aracne3_logger)
      aracne3_logger->writeLineWithTime(cur_msg);
    return;
  };

  auto qexit = [&](const std::string &cur_err) {
    std::cout << M << cur_err << std::endl;
    if (aracne3_logger)
      aracne3_logger->writeLineWithTime(cur_err);
    std::exit(EXIT_FAILURE);
  };

  auto qlog_subnet = [&](const uint16_t subnet_number,
                         const uint32_t subnet_size) {
    const std::string cur_msg =
        "...subnetwork " + std::to_string(subnet_number) +
        " completed = " + std::to_string(subnet_size) + " edges returned...";
    qlog(cur_msg);
    return;
  };

  // ---- Begin ARACNe3 instance ----

  std::mt19937 rnd{seed};

  qlog("Beginning ARACNe3 instance...");

  // ---- Reading input files ----

  vv_float exp_mat;
  compression_map compressor;
  decompression_map decompressor;
  gene_to_geneset regulons;
  geneset genes, regulators;

  qlog("Processing inputs...");
  Stopwatch watch1{};

  try {
    if (aracne3_logger)
      aracne3_logger->writeLineWithTime("...processing expression matrix...");

    std::tie(exp_mat, genes, compressor, decompressor) =
        io.readExpMatrixAndCopulaTransform(rnd, aracne3_logger);

    if (aracne3_logger)
      aracne3_logger->writeLineWithTime("...processing regulators...");

    regulators = io.readRegList(compressor, aracne3_logger, verbose);
  } catch (const std::exception &e) {
    std::string err_msg = std::string("Error processing input: ") + e.what();

    std::cerr << err_msg << std::endl;
    if (aracne3_logger)
      aracne3_logger->writeLineWithTime(err_msg);

    throw; // re-throw the exception for natural program termination
  }

  qlog("Inputs processed. Time elapsed: " + watch1.getSeconds());

  // ---- Get the null model for mutual information ----

  qlog("Getting null model for mutual information by adaptive partitioning...");

  const uint16_t n_samps = exp_mat.at(0).size();
  const uint16_t n_subsamp = std::ceil(n_samps * subsamp_pct);

  constexpr uint32_t null_model_seed = 0u;
  const std::string cached_blob_name =
      cached_dir + "APMINullModel_" + std::to_string(n_subsamp) + "_" +
      std::to_string(n_nulls) + "_" + std::to_string(null_model_seed) + "_" +
      version + ".blob";

  APMINullModel apmi_null_model;

#ifdef _DEBUG // Avoid caches when debugging
  qlog("Debug build: Generating new null model for subsampled APMI...");
  watch1.reset();

  apmi_null_model = APMINullModel(n_subsamp, n_nulls, null_model_seed);

  qlog("Null model generated. Time elapsed: " + watch1.getSeconds());

  std::string null_model_dir = output_dir + "null-model/";

  qlog("Debug build: Printing null model in \'" + null_model_dir + "\'...");

  makeDirs(null_model_dir, aracne3_logger.get());

  std::ofstream mi_ofs(null_model_dir + "null-mis.csv");
  std::ofstream stats_ofs(null_model_dir + "null-model-stats.csv");

  if (!mi_ofs || !stats_ofs)
    qexit("Error: could not write to files in \'" + null_model_dir + "\'.");

  auto [null_mis, ols_m, ols_b, n_samps_nm, n_nulls_nm, seed_nm] =
      apmi_null_model.getModel();

  for (const float mi : null_mis)
    mi_ofs << mi << ',';

  stats_ofs << "ols_m" << ',' << ols_m << '\n';
  stats_ofs << "ols_b" << ',' << ols_b << '\n';
  stats_ofs << "n_samps" << ',' << n_samps_nm << '\n';
  stats_ofs << "n_nulls" << ',' << n_nulls_nm << '\n';
  stats_ofs << "seed" << ',' << seed_nm << '\n';

  mi_ofs.close();
  stats_ofs.close();

  qlog("Debug build: Null model printed.");

#else
  if (std::filesystem::exists(cached_blob_name)) {
    qlog("Cached null model for subsampled APMI found. Reading in null model "
         "(n = " +
         std::to_string(n_subsamp) + ")...");
    watch1.reset();

    apmi_null_model = APMINullModel::getCachedModel(cached_blob_name);

    qlog("Null model read. Time elapsed: " + watch1.getSeconds());

  } else {
    qlog("Generating new null model for subsampled APMI...");
    watch1.reset();

    apmi_null_model = APMINullModel(n_subsamp, n_nulls, null_model_seed);

    qlog("Null model generated. Time elapsed: " + watch1.getSeconds());

    qlog("Caching null model for subsampled APMI (n = " +
         std::to_string(n_subsamp) + ")...");
    watch1.reset();

    apmi_null_model.cacheModel(cached_blob_name);

    qlog("Null model cached. Time elapsed: " + watch1.getSeconds());
  }
#endif /* _DEBUG */

  // ---- Begin analysis ----

  std::vector<gene_to_gene_to_float> subnets;
  std::vector<uint32_t> subnet_sizes;
  std::vector<float> FPR_estimates;
  float FPR_estimate = 1.5e-4f;

  if (!consolidate_mode) {
    qlog("Creating subnetwork(s)...");
    watch1.reset();

    if (adaptive) {
      gene_to_geneset regulons(regulators.size());

      bool stoppingCriteriaMet = false;
      uint16_t cur_subnet_ct = 0u;

      while (!stoppingCriteriaMet) {
        vv_float subsample_exp_mat =
            sampleExpMatAndReCopulaTransform(exp_mat, n_subsamp, rnd);

        auto [subnet, FPR_estimate_subnet, subnet_size] = createARACNe3Subnet(
            subsample_exp_mat, regulators, genes, n_samps, n_subsamp,
            subnets.size() + 1u, prune_alpha, apmi_null_model, method, alpha,
            prune_MaxEnt, subnets_log_dir, threads, runid, decompressor,
            save_subnets, io);
        subnets.push_back(std::move(subnet));
        FPR_estimates.push_back(FPR_estimate_subnet);

        if (subnets[subnets.size() - 1u].size() == 0u)
          qexit("Abort: No edges returned. Empty subnetwork.");

        // add any new edges to the regulon_set
        for (const auto &[reg, regulon] : subnets[subnets.size() - 1u])
          for (const auto [tar, mi] : regulon)
            regulons[reg].insert(tar);

        // check minimum regulon size
        uint16_t min_regulon_size = 65535u;
        for (const auto &[reg, regulon] : regulons)
          if (regulons.at(reg).size() < min_regulon_size)
            min_regulon_size = regulons.at(reg).size();

        if (min_regulon_size >= min_regulon_occpuancy &&
            subnets.size() >= min_subnets)
          stoppingCriteriaMet = true;

        qlog_subnet(subnets.size(), subnet_size);

        if (subnets.size() >= max_subnets)
          qexit("Abort: Max subnets reached.");
      }
      n_subnets = subnets.size();
    } else if (!adaptive) {
      subnets = std::vector<gene_to_gene_to_float>(n_subnets);
      subnet_sizes = std::vector<uint32_t>(n_subnets);
      FPR_estimates = std::vector<float>(n_subnets);

      for (uint16_t i = 0u; i < n_subnets; ++i) {
        vv_float subsample_exp_mat =
            sampleExpMatAndReCopulaTransform(exp_mat, n_subsamp, rnd);

        std::tie(subnets.at(i), FPR_estimates.at(i), subnet_sizes.at(i)) =
            createARACNe3Subnet(subsample_exp_mat, regulators, genes, n_samps,
                                n_subsamp, i + 1u, prune_alpha, apmi_null_model,
                                method, alpha, prune_MaxEnt, subnets_log_dir,
                                threads, runid, decompressor, save_subnets, io);

        if (subnets.at(i).size() == 0)
          qexit("Abort: No edges returned. Empty subnetwork.");

        qlog_subnet(i + 1u, subnet_sizes.at(i));
      }
    }

    qlog("Subnetwork generation complete. Time elapsed: " +
         watch1.getSeconds());
    qlog("Total subnetworks generated: " + std::to_string(n_subnets));
  } else if (consolidate_mode) {
    qlog("Skip to consolidate step requested. Attempting to read subnetworks "
         "from directory \"" +
         subnets_dir + "\" and subnetwork logging information from \"" +
         subnets_log_dir +
         "\" (subnetwork logging information is required)...");
    watch1.reset();

    const auto &[subnet_filenames, subnet_log_filenames] =
        io.findSubnetFilesAndSubnetLogFiles(subnets_dir, subnets_log_dir,
                                            aracne3_logger);

    if (subnet_filenames.size() < n_subnets)
      qexit("Error: Too many subnets requested. Only " +
            std::to_string(subnet_filenames.size()) + " subnets found in \"" +
            subnets_dir + "\".");

    for (uint16_t subnet_idx = 0; subnet_idx < n_subnets; ++subnet_idx) {
      const auto &[subnet, FPR_estimate_subnet] =
          io.loadARACNe3SubnetsAndUpdateFPRFromLog(
              subnets_dir + subnet_filenames[subnet_idx],
              subnets_log_dir + subnet_log_filenames[subnet_idx], compressor,
              regulators, aracne3_logger);
      subnets.push_back(subnet);
      FPR_estimates.push_back(FPR_estimate_subnet);
    }

    n_subnets = subnets.size();

    qlog("Subnetworks read. Time elapsed: " + watch1.getSeconds());
    qlog("Total subnetworks read: " + std::to_string(n_subnets));
  }

  // set the FPR estimate
  FPR_estimate =
      std::accumulate(FPR_estimates.begin(), FPR_estimates.end(), 0.0f) /
      FPR_estimates.size();

  if (!skip_consolidate) {
    qlog("Consolidating subnetworks...");
    watch1.reset();

    std::vector<ARACNe3_df> final_df = consolidateSubnetsVec(
        subnets, FPR_estimate, exp_mat, regulators, genes, rnd);

    qlog("Consolidation complete. Time elapsed: " + watch1.getSeconds());
    qlog("Total subnetworks consolidated: " + std::to_string(n_subnets));

    qlog("Writing final network...");
    watch1.reset();

    io.writeARACNe3DF(final_df, decompressor);

  } else if (skip_consolidate) {
    qlog("No consolidation requested.");
  }

  qlog("Success!");

  return EXIT_SUCCESS;
}

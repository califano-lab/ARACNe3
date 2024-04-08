#include "config.h"

#include "apmi_nullmodel.hpp"
#include "cmdline_parser.hpp"
#include "fsio.hpp"
#include "logger.hpp"
#include "stopwatch.hpp"
#include "subnet_operations.hpp"

#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>

// ---- Fetch app information ----
#ifndef CACHE_PATH
#define CACHE_PATH "cache/aracne3"
#endif

#ifndef APP_VERSION
#define APP_VERSION "0.0.0"
#endif

/**
 * Wraps STL functions to create a full directory path; this is mainly a
 * convenience function
 *
 * @param dir_name A constant reference to a string representing the directory
 * path.
 */
bool makeDirs(const std::string &dir_name, Logger *const logger) {
  if (!std::filesystem::exists(dir_name)) {
    std::filesystem::create_directories(dir_name);
    if (std::filesystem::exists(dir_name)) {
      const std::string out_msg = "Directory Created: \"" + dir_name + "\".";

      std::cout << out_msg << std::endl;
      if (logger)
        logger->writeLineWithTime(out_msg);

      return true;
    } else {
      const std::string err_msg =
          "Failed to create directory: \"" + dir_name + "\".";

      std::cerr << err_msg << std::endl;
      if (logger)
        logger->writeLineWithTime(err_msg);

      return false;
    }
  }
  return true;
}

// ---- Begin ARACNe3 runtime ----

int main(int argc, char *argv[]) {

  // ---- Initialize ARACNe3 runtime variables ----

  uint32_t seed = static_cast<uint32_t>(std::time(nullptr));
  bool seed_provided = false;
  uint8_t threads = 1U;
  std::string runid = "defaultid";
  bool verbose = false, suppress_log = false;
  std::string cached_dir = CACHE_PATH;

  uint16_t n_subnets = 1u;
  float subsamp_pct = 1 - std::exp(-1);
  bool skip_consolidate = false;
  bool consolidate_mode = false;
  bool adaptive = false;
  float alpha = 0.05f;
  bool prune_alpha = true;
  bool prune_MaxEnt = true;
  bool save_subnets = false;
  uint16_t min_regulon_occpuancy = 30U;
  uint16_t min_subnets = 0u;
  uint16_t max_subnets = 65535u;
  std::string method = "FDR";

  std::unique_ptr<Logger> aracne3_logger;

  float mi_cutoff = 0.f;
  uint32_t n_nulls = 1'000'000u;

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

  // ---- Capturing command line arguments ----

  CmdLineParser clp(argc, argv);

  if (clp.optExists("-h") || clp.optExists("--help") || !clp.optExists("-e") ||
      !clp.optExists("-r") || !clp.optExists("-o")) {
    std::cout << M
              << "usage: " + ((std::string)argv[0]) +
                     " -e path/to/matrix.txt -r path/to/regulators.txt -o "
                     "path/to/output/directory"
              << std::endl;
    return EXIT_FAILURE;
  }

  const std::string exp_mat_file_path = clp.getOpt("-e");
  const std::string regulators_list_file_path = clp.getOpt("-r");
  std::string output_dir = clp.getOpt("-o");

  try {
    if (clp.optExists("--seed"))
      seed = std::stoi(clp.getOpt("--seed"));
    if (clp.optExists("--threads"))
      threads = std::stoi(clp.getOpt("--threads"));
    if (clp.optExists("--runid"))
      runid = clp.getOpt("--runid");
    if (clp.optExists("--verbose") || clp.optExists("-v"))
      verbose = true;
    if (clp.optExists("--suppress-log") | clp.optExists("--suppress-logs"))
      suppress_log = true;

    if (clp.optExists("--alpha"))
      alpha = std::stof(clp.getOpt("--alpha"));
    if (clp.optExists("--subsample"))
      subsamp_pct = std::stof(clp.getOpt("--subsample"));

    if (clp.optExists("-x"))
      n_subnets = min_regulon_occpuancy = std::stoi(clp.getOpt("-x"));

    if (clp.optExists("--threads"))
      threads = std::stoi(clp.getOpt("--threads"));

    if (clp.optExists("--skip-alpha"))
      prune_alpha = false;
    if (clp.optExists("--skip-maxent"))
      prune_MaxEnt = false;

    if (clp.optExists("--FDR"))
      method = "FDR";
    if (clp.optExists("--FWER"))
      method = "FWER";
    if (clp.optExists("--FPR"))
      method = "FPR";

    if (clp.optExists("--adaptive"))
      adaptive = true;
    if (clp.optExists("--min-subnets"))
      min_subnets = std::stoi(clp.getOpt("--min-subnets"));
    if (clp.optExists("--min-subnetworks"))
      min_subnets = std::stoi(clp.getOpt("--min-subnetworks"));
    if (clp.optExists("--max-subnets"))
      max_subnets = std::stoi(clp.getOpt("--max-subnets"));
    if (clp.optExists("--max-subnetworks"))
      max_subnets = std::stoi(clp.getOpt("--max-subnetworks"));
    if (clp.optExists("--save-subnetworks"))
      save_subnets = true;
    if (clp.optExists("--skip-consolidate"))
      skip_consolidate = save_subnets = true;
    if (clp.optExists("--consolidate-mode") || clp.optExists("--consolidate"))
      consolidate_mode = true;

    // ---- Developer options ----
    if (clp.optExists("--mithresh"))
      mi_cutoff = std::stof(clp.getOpt("--mithresh"));
    if (clp.optExists("--numnulls"))
      n_nulls = std::stoi(clp.getOpt("--numnulls"));

    // ---- Checking edge cases ----
    if (alpha > 1.f || alpha <= 0.f)
      throw std::runtime_error("--alpha must be on range [0,1]");
    if (subsamp_pct > 1.f || subsamp_pct <= 0.f)
      throw std::runtime_error("--subsample must be on range [0,1]");

    if (skip_consolidate && consolidate_mode)
      throw std::runtime_error(
          "Cannot simultaneously skip consolidate and enter consolidate mode");

    if (n_nulls < 0u)
      throw std::runtime_error("--numnulls be greater than 0");

  } catch (const std::exception &e) {
    std::string err_msg =
        std::string("Error parsing command line option: ") + e.what();

    std::cerr << err_msg << std::endl;

    throw; // re-throw the exception for natural program termination
  }

  // ---- Processing command line arguments ----

  if (output_dir.back() != '/')
    output_dir += '/';

  if (cached_dir.back() != '/')
    cached_dir += '/';

  makeDirs(output_dir, aracne3_logger.get());

  const std::string log_file_name = output_dir + "log_" + runid + ".txt";
  if (!suppress_log)
    aracne3_logger = std::make_unique<Logger>(log_file_name, argc, argv);

  makeDirs(cached_dir, aracne3_logger.get());

  const std::string subnets_dir = output_dir + "subnetworks/";
  const std::string subnets_log_dir = output_dir + "log-subnetworks/";

  if (save_subnets) {
    makeDirs(subnets_dir, aracne3_logger.get());
    makeDirs(subnets_log_dir, aracne3_logger.get());
  }

  // ---- Begin ARACNe3 instance ----

  std::mt19937 rnd{seed};

  FilesystemIOHandler io(exp_mat_file_path, regulators_list_file_path,
                         output_dir + "network_" + runid + ".tsv", subnets_dir,
                         runid, '\t');

  qlog("See logs and progress reports in \"" + log_file_name + "\"");

  // ---- Reading input files ----

  vv_float exp_mat;
  compression_map compressor;
  decompression_map decompressor;
  gene_to_geneset regulons;
  geneset genes, regulators;

  qlog("Reading input files...");
  Stopwatch watch1{};

  try {
    if (aracne3_logger)
      aracne3_logger->writeLineWithTime("...processing expression matrix...");
    std::tie(exp_mat, genes, compressor, decompressor) =
        io.readExpMatrixAndCopulaTransform(rnd, aracne3_logger.get());

    if (aracne3_logger)
      aracne3_logger->writeLineWithTime("...processing regulators...");
    regulators = io.readRegList(compressor, aracne3_logger.get(), verbose);
  } catch (const std::exception &e) {
    std::string err_msg = std::string("Error reading input files: ") + e.what();

    std::cerr << err_msg << std::endl;
    if (aracne3_logger)
      aracne3_logger->writeLineWithTime(err_msg);

    throw; // re-throw the exception for natural program termination
  }

  qlog("Input files read. Time elapsed: " + watch1.getSeconds());

  // ---- Get the null model for mutual information ----

  qlog("Getting null model for mutual information by adaptive partitioning...");

  const uint16_t n_samps = exp_mat.at(0).size();
  const uint16_t n_subsamp = std::ceil(n_samps * subsamp_pct);

  constexpr uint32_t null_model_seed = 0u;
  const std::string cached_blob_name =
      cached_dir + "APMINullModel_" + std::to_string(n_subsamp) + "_" +
      std::to_string(n_nulls) + "_" + std::to_string(null_model_seed) + "_" +
      APP_VERSION + ".blob";

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
                                            aracne3_logger.get());

    if (subnet_filenames.size() < n_subnets)
      qexit("Error: Too many subnets requested. Only " +
            std::to_string(subnet_filenames.size()) + " subnets found in \"" +
            subnets_dir + "\".");

    for (uint16_t subnet_idx = 0; subnet_idx < n_subnets; ++subnet_idx) {
      const auto &[subnet, FPR_estimate_subnet] =
          io.loadARACNe3SubnetsAndUpdateFPRFromLog(
              subnets_dir + subnet_filenames[subnet_idx],
              subnets_log_dir + subnet_log_filenames[subnet_idx], compressor,
              regulators, aracne3_logger.get());
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

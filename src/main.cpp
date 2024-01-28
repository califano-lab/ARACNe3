#include "config.h"

#include "ARACNe3.hpp"
#include "apmi_nullmodel.hpp"
#include "cmdline_parser.hpp"
#include "io.hpp"
#include "logger.hpp"
#include "stopwatch.hpp"
#include "subnet_operations.hpp"

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp> // for object caching
#include <boost/serialization/vector.hpp>

// ---- Fetch app information ----
#ifndef CACHE_PATH
#define CACHE_PATH "cache/aracne3"
#endif

#ifndef APP_VERSION
#define APP_VERSION "0.0.0"
#endif

// ---- Define APMINullModel serialization fn for Boost ----
namespace boost {
namespace serialization {

template <class Archive>
void serialize(Archive &ar, APMINullModel &model, const unsigned int version) {
  auto [null_mis, ols_m, ols_b, n_samps, n_nulls, seed] = model.getModel();

  ar &null_mis;
  ar &ols_m;
  ar &ols_b;
  ar &n_samps;
  ar &n_nulls;
  ar &seed;

  // For deserialization, set the values back into the model
  if (Archive::is_loading::value) {
    model = APMINullModel(null_mis, ols_m, ols_b, n_samps, n_nulls, seed);
  }
}
} // namespace serialization
} // namespace boost

// ---- Begin ARACNe3 runtime ----

extern std::vector<std::string> decompression_map;

int main(int argc, char *argv[]) {

  // ---- Initialize ARACNe3 runtime variables ----

  uint32_t seed = static_cast<uint32_t>(std::time(nullptr));
  bool seed_provided = false;
  uint8_t threads = 1U;
  std::string runid = "defaultid";
  bool verbose = false, suppress_logs = false;
  std::string cached_dir = CACHE_PATH;

  uint16_t n_subnets = 1u;
  float subsamp_pct = 1 - std::exp(-1);
  bool skip_consolidate = false;
  bool consolidate_mode = false;
  bool adaptive = false;
  float alpha = 0.05f;
  bool prune_alpha = true;
  bool prune_MaxEnt = true;
  uint16_t min_regulon_occpuancy = 30U;
  uint16_t min_subnets = 0u;
  std::string method = "FDR";

  std::unique_ptr<Logger> aracne3_logger;
  const std::string M = " -A3- ";

  float mi_cutoff = 0.f;
  uint32_t n_nulls = 1'000'000u;

  // ---- Quick logging macro ----

  auto qlog = [&](const std::string &cur_msg) {
    std::cout << M << cur_msg << std::endl;
    if (aracne3_logger)
      aracne3_logger->writeLineWithTime(cur_msg);
    return;
  };

  // ---- Capturing command line arguments ----

  CmdLineParser clp(argc, argv);

  if (clp.optExists("-h") || clp.optExists("--help") || !clp.optExists("-e") ||
      !clp.optExists("-r") || !clp.optExists("-o")) {
    std::cout << M
              << "usage: " + ((std::string)argv[0]) +
                     makeUnixDirectoryNameUniversal(
                         " -e path/to/matrix.txt -r path/to/regulators.txt -o "
                         "path/to/output/directory")
              << std::endl;
    return EXIT_FAILURE;
  }

  const std::string exp_mat_file_path = clp.getOpt("-e");
  const std::string reg_list_file_path = clp.getOpt("-r");
  std::string output_dir = clp.getOpt("-o");

  if (clp.optExists("--seed"))
    seed = std::stoi(clp.getOpt("--seed"));
  if (clp.optExists("--threads"))
    threads = std::stoi(clp.getOpt("--threads"));
  if (clp.optExists("--runid"))
    runid = clp.getOpt("--runid");
  if (clp.optExists("--verbose") || clp.optExists("-v"))
    verbose = true;
  if (clp.optExists("--suppress-logs"))
    suppress_logs = true;

  if (clp.optExists("--alpha"))
    alpha = std::stof(clp.getOpt("--alpha"));
  if (clp.optExists("--subsample"))
    subsamp_pct = std::stof(clp.getOpt("--subsample"));

  if (alpha > 1.f || alpha <= 0)
    alpha = 1.f;
  if (subsamp_pct > 1.f || subsamp_pct <= 0)
    subsamp_pct = 1.f;

  if (clp.optExists("-x"))
    n_subnets = min_regulon_occpuancy = std::stoi(clp.getOpt("-x"));

  if (clp.optExists("--threads"))
    threads = std::stoi(clp.getOpt("--threads"));

  if (clp.optExists("--noalpha"))
    prune_alpha = false;
  if (clp.optExists("--noMaxEnt"))
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
  if (clp.optExists("--skipconsolidate"))
    skip_consolidate = true;
  if (clp.optExists("--consolidate"))
    consolidate_mode = true;

  // ---- Developer options ----

  if (clp.optExists("--mithresh"))
    mi_cutoff = std::stof(clp.getOpt("--mithresh"));
  if (clp.optExists("--numnulls"))
    n_nulls = std::stoi(clp.getOpt("--numnulls"));

  if (n_nulls < 0u)
    n_nulls = 1'000'000u;

  // ---- Processing command line arguments ----

  if (output_dir.back() != '/')
    output_dir += '/';

  if (cached_dir.back() != '/')
    cached_dir += '/';

  makeDirs(output_dir, aracne3_logger.get());

  const std::string log_file_name = output_dir + "log_" + runid + ".txt";
  if (!suppress_logs)
    aracne3_logger = std::make_unique<Logger>(log_file_name, argc, argv);

  makeDirs(cached_dir, aracne3_logger.get());

  const std::string subnets_dir =
      makeUnixDirectoryNameUniversal(output_dir + "subnets/");
  const std::string subnets_log_dir =
      makeUnixDirectoryNameUniversal(output_dir + "subnets_log/");

  makeDirs(subnets_dir, aracne3_logger.get());
  makeDirs(subnets_log_dir, aracne3_logger.get());

  // ---- Begin ARACNe3 instance ----

  std::mt19937 rand{seed};

  std::string cur_msg = "Beginning ARACNe3 instance...";
  std::cout << M << cur_msg << std::endl;

  if (aracne3_logger) {
    std::cout << M
              << "See logs and progress reports in \"" + log_file_name + "\"."
              << std::endl;
    aracne3_logger->writeLineWithTime(cur_msg);
  }

  // ---- Reading input files ----

  qlog("Reading input files...");
  Stopwatch watch1{};

  auto data = readExpMatrixAndCopulaTransform(exp_mat_file_path, rand);
  const gene_to_floats &exp_mat = std::get<0>(data);
  const gene_to_shorts &ranks_mat = std::get<1>(data);
  const geneset &genes = std::get<2>(data);
  const uint16_t tot_n_samps = std::get<3>(data);

  uint16_t tot_n_subsample = std::ceil(subsamp_pct * tot_n_samps);
  if (tot_n_subsample >= tot_n_samps || tot_n_subsample < 0) {
    std::cerr
        << "Warning: subsample quantity invalid. All samples will be used."
        << std::endl;
    tot_n_subsample = tot_n_samps;
  }

  const geneset regulators = readRegList(reg_list_file_path, verbose);

  qlog("Input files read. Time elapsed: " + watch1.getSeconds());

  // ---- Get the null model for mutual information ----

  qlog("Getting null model for mutual information by adaptive partitioning...");

  const size_t n_samps = exp_mat.at(0).size();
  constexpr uint32_t null_model_seed = 0u;
  const std::string cached_blob_name =
      cached_dir + "APMINullModel_" + std::to_string(n_samps) + "_" +
      std::to_string(n_nulls) + "_" + std::to_string(null_model_seed) + "_" +
      APP_VERSION + ".blob";

  APMINullModel apmi_null_model;

#ifdef _DEBUG // Avoid caches when debugging
  qlog("Debug build: Generating new null model...");
  watch1.reset();

  apmi_null_model = APMINullModel(n_samps, n_nulls, seed);

  qlog("Null model generated. Time elapsed: " + watch1.getSeconds());
#else
  if (std::filesystem::exists(cached_blob_name)) {
    qlog("Cached null model found. Reading in null model...");
    watch1.reset();

    std::ifstream ifs(cached_blob_name, std::ios::binary);
    if (ifs) {
      boost::archive::binary_iarchive ia(ifs);
      ia >> apmi_null_model;
    }
    qlog("Null model read. Time elapsed: " + watch1.getSeconds());

  } else {
    qlog("Generating new null model...");
    watch1.reset();

    apmi_null_model = APMINullModel(n_samps, n_nulls, null_model_seed);

    qlog("Null model generated. Time elapsed: " + watch1.getSeconds());

    qlog("Caching null model...");
    watch1.reset();
    std::ofstream ofs(cached_blob_name, std::ios::binary);
    boost::archive::binary_oarchive oa(ofs);
    oa << apmi_null_model;
    qlog("Null model cached. Time elapsed: " + watch1.getSeconds());
  }
#endif /* _DEBUG */

  // ---- Begin subnetwork generation ----

  std::vector<gene_to_gene_to_float> subnets;
  std::vector<float> FPR_estimates;
  float FPR_estimate = 1.5e-4f;

  if (!consolidate_mode) {

    qlog("Creating subnetwork(s)...");
    watch1.reset();

    if (adaptive) {
      gene_to_geneset regulons(regulators.size());

      bool stoppingCriteriaMet = false;
      uint16_t cur_subnet_ct = 0;

      while (!stoppingCriteriaMet) {
        gene_to_floats subsample_exp_mat =
            sampleExpMatAndReCopulaTransform(exp_mat, tot_n_subsample, rand);

        const auto &[subnet, FPR_estimate_subnet] = createARACNe3Subnet(
            subsample_exp_mat, regulators, genes, tot_n_samps, tot_n_subsample,
            cur_subnet_ct, prune_alpha, apmi_null_model, method, alpha,
            prune_MaxEnt, output_dir, subnets_dir, subnets_log_dir, threads,
            runid);

        subnets.push_back(subnet);
        FPR_estimates.push_back(FPR_estimate_subnet);

        if (subnet.size() == 0) {
          std::cerr << "Abort: No edges left after all pruning steps. Empty "
                       "subnetwork."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }

        // add any new edges to the regulon_set
        for (const auto [reg, tar_mi] : subnet)
          for (const auto [tar, mi] : tar_mi)
            regulons[reg].insert(tar);

        // Check minimum regulon size
        uint16_t min_regulon_size = 65535U;
        for (const auto &[reg, regulon] : regulons)
          if (regulons[reg].size() < min_regulon_size)
            min_regulon_size = regulons[reg].size();

        ++cur_subnet_ct;

        if (min_regulon_size >= min_regulon_occpuancy &&
            cur_subnet_ct >= min_subnets)
          stoppingCriteriaMet = true;
      }
      n_subnets = subnets.size();
    } else if (!adaptive) {
      subnets = std::vector<gene_to_gene_to_float>(n_subnets);
      FPR_estimates = std::vector<float>(n_subnets);
      for (int i = 0; i < n_subnets; ++i) {
        gene_to_floats subsample_exp_mat =
            sampleExpMatAndReCopulaTransform(exp_mat, tot_n_subsample, rand);

        std::tie(subnets[i], FPR_estimates[i]) = createARACNe3Subnet(
            subsample_exp_mat, regulators, genes, tot_n_samps, tot_n_subsample,
            i, prune_alpha, apmi_null_model, method, alpha, prune_MaxEnt,
            output_dir, subnets_dir, subnets_log_dir, threads, runid);
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
         "\" (subnetwork logging information is required) ...");
    watch1.reset();

    const auto &[subnet_filenames, subnet_log_filenames] =
        findSubnetFilesAndSubnetLogFiles(subnets_dir, subnets_log_dir);

    if (subnet_filenames.size() < n_subnets) {
      std::cerr << "Error: Too many subnets requested. Only " +
                       std::to_string(subnet_filenames.size()) +
                       " subnets found in \"" + subnets_dir + "\"."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }

    for (uint16_t subnet_idx = 0; subnet_idx < n_subnets; ++subnet_idx) {
      const auto &[subnet, FPR_estimate_subnet] =
          loadARACNe3SubnetsAndUpdateFPRFromLog(
              subnets_dir + subnet_filenames[subnet_idx],
              subnets_log_dir + subnet_log_filenames[subnet_idx]);
      subnets.push_back(subnet);
      FPR_estimates.push_back(FPR_estimate_subnet);
    }

    n_subnets = subnets.size();

    qlog("Subnetworks read. Time elapsed: " + watch1.getSeconds());
    qlog("Total subnets read: " + std::to_string(n_subnets));
  }

  // set the FPR estimate
  FPR_estimate =
      std::accumulate(FPR_estimates.begin(), FPR_estimates.end(), 0.0f) /
      FPR_estimates.size();

  if (!skip_consolidate) {

    qlog("Consolidating subnetworks...");
    watch1.reset();

    std::vector<consolidated_df_row> final_df = consolidateSubnetsVec(
        subnets, FPR_estimate, exp_mat, regulators, genes, ranks_mat);

    qlog("Consolidation complete. Time elapsed: " + watch1.getSeconds());
    qlog("Total subnets consolidated: " + std::to_string(n_subnets));

    qlog("Writing final network...");
    watch1.reset();

    writeConsolidatedNetwork(final_df,
                             output_dir + "consolidated-net_" + runid + ".tsv");

  } else if (skip_consolidate) {

    qlog("No consolidation requested.");
  }

  qlog("Success!");

  return EXIT_SUCCESS;
}

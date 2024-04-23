#include "analysis.hpp"
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
  uint8_t threads = 1U;
  std::string runid = "defaultid";
  bool verbose = false, suppress_log = false;
  
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

  uint32_t n_nulls = 1'000'000u;

  // ---- Quick macros ----

  const std::string M = " -A3- ";
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
    if (clp.optExists("--suppress-log") || clp.optExists("--suppress-logs"))
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

  makeDirs(output_dir, aracne3_logger.get());

  const std::string log_file_name = output_dir + "log_" + runid + ".txt";
  if (!suppress_log)
    aracne3_logger = std::make_unique<Logger>(log_file_name, argc, argv);

  const std::string subnets_dir = output_dir + "subnetworks/";
  const std::string subnets_log_dir = output_dir + "log-subnetworks/";

  if (save_subnets) {
    makeDirs(subnets_dir, aracne3_logger.get());
    makeDirs(subnets_log_dir, aracne3_logger.get());
  }

  // ---- Begin ARACNe3 instance ----

  FilesystemIOHandler io(exp_mat_file_path, regulators_list_file_path,
                         output_dir + "network_" + runid + ".tsv", subnets_dir,
                         runid, '\t');

  if(!suppress_log)
    qlog("See logs and progress reports in \"" + log_file_name + "\"");

#if _DEBUG  // when debug, cache null model in output dir
  std::string cached_dir = output_dir;
#else
  std::string cached_dir = CACHE_PATH;
  makeDirs(cached_dir, aracne3_logger.get());
#endif

  ARACNe3Analysis(io, APP_VERSION, runid, seed, threads, verbose, alpha,
                  subsamp_pct, n_nulls, method, prune_alpha, prune_MaxEnt,
                  save_subnets, subnets_log_dir, n_subnets, adaptive,
                  min_subnets, max_subnets, min_regulon_occpuancy,
                  consolidate_mode, skip_consolidate, aracne3_logger.get(),
                  subnets_dir, cached_dir);

  return EXIT_SUCCESS;
}

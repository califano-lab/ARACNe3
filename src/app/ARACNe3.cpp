#include <omp.h>
#include <iostream>
#include <fstream>
#include "ARACNe3.hpp"
#include "cmdline_parser.hpp"
#include "stopwatch.hpp"
#include "apmi_nullmodel.hpp"
#include "subnet_operations.hpp"
#include "io.hpp"

uint16_t nthreads = 1U;

extern std::vector<std::string> decompression_map;

/*
 Main function is the command line executable; this primes the global variables and parses the command line.  It will also return usage notes if the user incorrectly calls ./ARACNe3.
 
 Example:
 ./ARACNe3 -e test/matrix.txt -r test/regulators.txt -o test/output --noAlpha -a 0.05 --alpha 0.05 --noMaxEnt --subsample 0.6321 --seed 1 --mithresh 0.2 --numnulls 1000000
 */
int main(int argc, char *argv[]) {

  //--------------------check requirements------------------------

  if (cmdOptionExists(argv, argv + argc, "-h") ||
      cmdOptionExists(argv, argv + argc, "--help") ||
      !cmdOptionExists(argv, argv + argc, "-e") ||
      !cmdOptionExists(argv, argv + argc, "-r") ||
      !cmdOptionExists(argv, argv + argc, "-o")) {
    std::cout << "usage: " + ((std::string)argv[0]) +
                     makeUnixDirectoryNameUniversal(
                         " -e path/to/matrix.txt -r path/to/regulators.txt -o "
                         "path/to/output/directory")
              << std::endl;
    return EXIT_FAILURE;
  }
  
  //--------------------initialize parameters---------------------

  uint16_t num_subnets = 1U;
  double subsampling_percent = 1 - std::exp(-1);
  bool do_not_consolidate = false;
  bool go_to_consolidate = false;
  bool adaptive = false;
  float alpha = 0.05f;
  bool prune_alpha = true;
  bool prune_MaxEnt = true;
  uint32_t seed = 0U;
  uint16_t num_subnets_to_consolidate = 0U;
  uint16_t targets_per_regulator = 30U;
  std::string runid = "defaultid";
  std::string method = "FDR";

  float DEVELOPER_mi_cutoff = 0.0f;
  uint32_t DEVELOPER_num_null_marginals = 1000000U;

  //--------------------parsing filesystem------------------------
  
  const std::string exp_mat_file = makeUnixDirectoryNameUniversal(
      (std::string)getCmdOption(argv, argv + argc, "-e"));
  const std::string reg_list_file = makeUnixDirectoryNameUniversal(
      (std::string)getCmdOption(argv, argv + argc, "-r"));
  std::string output_dir = makeUnixDirectoryNameUniversal(
      (std::string)getCmdOption(argv, argv + argc, "-o"));

  // make sure output_dir has a trailing slash
  if (output_dir.back() != directory_slash)
    output_dir += directory_slash;

  const std::string cached_dir = makeUnixDirectoryNameUniversal("./" + hiddenfpre + "ARACNe3_cached/");
  //--------------------parsing parameters------------------------

  if (cmdOptionExists(argv, argv + argc, "--alpha"))
    alpha = std::stof(getCmdOption(argv, argv + argc, "--alpha"));
  if (alpha > 1.f || alpha <= 0) {
    std::cout << "alpha not on range [0,1], setting to 1." << std::endl;
    alpha = 1.f;
  }

  if (cmdOptionExists(argv, argv + argc, "--seed"))
    seed = std::stoi(getCmdOption(argv, argv + argc, "--seed"));

  if (cmdOptionExists(argv, argv + argc, "--subsample"))
    subsampling_percent =
        std::stod(getCmdOption(argv, argv + argc, "--subsample"));

  if (subsampling_percent > 1.f || subsampling_percent <= 0) {
    std::cout << "Subsampling percent not on range (0,1]; setting to 1."
              << std::endl;
    subsampling_percent = 1.f;
  }

  if (cmdOptionExists(argv, argv + argc, "-x"))
    num_subnets = targets_per_regulator = num_subnets_to_consolidate =
        std::stoi(getCmdOption(argv, argv + argc, "-x"));

  if (cmdOptionExists(argv, argv + argc, "--threads"))
    nthreads = std::stoi(getCmdOption(argv, argv + argc, "--threads"));

  if (cmdOptionExists(argv, argv + argc, "--noalpha"))
    prune_alpha = false;
  if (cmdOptionExists(argv, argv + argc, "--noMaxEnt"))
    prune_MaxEnt = false;
  if (cmdOptionExists(argv, argv + argc, "--FDR"))
    method = "FDR";
  if (cmdOptionExists(argv, argv + argc, "--FWER"))
    method = "FWER";
  if (cmdOptionExists(argv, argv + argc, "--FPR"))
    method = "FPR";
  if (cmdOptionExists(argv, argv + argc, "--adaptive"))
    adaptive = true;
  if (cmdOptionExists(argv, argv + argc, "--noConsolidate"))
    do_not_consolidate = true;
  if (cmdOptionExists(argv, argv + argc, "--consolidate"))
    go_to_consolidate = true;
  if (cmdOptionExists(argv, argv + argc, "--runid"))
    runid = getCmdOption(argv, argv + argc, "--runid");

  //--------------------developer parameters----------------------

  if (cmdOptionExists(argv, argv + argc, "--mithresh"))
    DEVELOPER_mi_cutoff =
        std::stof(getCmdOption(argv, argv + argc, "--mithresh"));
  if (DEVELOPER_mi_cutoff < 0)
    DEVELOPER_mi_cutoff = 0.0f;

  if (cmdOptionExists(argv, argv + argc, "--numnulls"))
    DEVELOPER_num_null_marginals =
        std::stoi(getCmdOption(argv, argv + argc, "--numnulls"));
  if (DEVELOPER_num_null_marginals < 0) {
    std::cout
        << "Number of null marginals not on range (0,inf); setting to 1000000."
        << std::endl;
    DEVELOPER_num_null_marginals = 1000000;
  }

  //--------------------------------------------------------------

  std::string log_filename = "log_" + runid + ".txt";
  std::string log_file_path = output_dir + log_filename;

  std::ofstream log_output(log_file_path);
  const std::string subnets_dir = makeUnixDirectoryNameUniversal(output_dir + "subnets_" + runid + "/");
  const std::string subnets_log_dir = makeUnixDirectoryNameUniversal(output_dir + "subnets_log_" + runid + "/");

  makeDir(output_dir);
  makeDir(cached_dir);
  makeDir(subnets_dir);
  makeDir(subnets_log_dir);
  std::mt19937 rand{seed};

  // print the initial command to the log output
  for (uint16_t i = 0; i < argc; ++i)
    log_output << std::string(argv[i]) << " ";
  log_output << std::endl;

  std::time_t t = std::time(nullptr);
  std::cout << "\n---------" << std::put_time(std::localtime(&t), "%c %Z")
            << "---------" << std::endl;
  log_output << "\n---------" << std::put_time(std::localtime(&t), "%c %Z")
             << "---------" << std::endl;

  std::cout
      << "Beginning ARACNe3 instance.  See logs and progress reports in \"" + log_file_path + "\"."
      << std::endl;
  log_output << "Beginning ARACNe3 instance..." << std::endl;

  Watch watch1;
  watch1.reset();

  log_output << "\nGene expression matrix & regulators list read time: ";

  auto data = readExpMatrixAndCopulaTransform(exp_mat_file, subsampling_percent, rand);
  const gene_to_floats& exp_mat = std::get<0>(data);
  const gene_to_shorts& ranks_mat = std::get<1>(data);
  const geneset& genes = std::get<2>(data);
  const uint16_t tot_num_samps = std::get<3>(data);

  uint16_t tot_num_subsample = std::ceil(subsampling_percent * tot_num_samps);
  if (tot_num_subsample >= tot_num_samps || tot_num_subsample < 0) {
    std::cerr
        << "Warning: subsample quantity invalid. All samples will be used."
        << std::endl;
    tot_num_subsample = tot_num_samps;
  }

  std::cout << "\nTotal N Samples: " + std::to_string(tot_num_samps)
            << std::endl;
  std::cout << "Subsampled N Samples: " + std::to_string(tot_num_subsample)
            << std::endl;

  const geneset regulators = readRegList(reg_list_file);

  //-------time module-------
  log_output << watch1.getSeconds() << std::endl;
  log_output << "\nMutual Information null model calculation time: ";
  watch1.reset();
  //-------------------------

  APMINullModel nullmodel = APMINullModel(DEVELOPER_num_null_marginals, tot_num_subsample, cached_dir, rand);
  nullmodel.cacheNullModel(cached_dir);

  //-------time module-------
  log_output << watch1.getSeconds() << std::endl;
  //-------------------------

  // Must exist regardless of whether we skip to consolidation
  std::vector<gene_to_gene_to_float> subnets;
  std::vector<float> FPR_estimates;
  float FPR_estimate = 1.5E-4f;

  if (!go_to_consolidate) {

    //-------time module-------
    log_output << "\nCreating subnetwork(s) time: ";
    watch1.reset();
    //-------------------------

    if (adaptive) {
      gene_to_geneset regulons;
      for (const uint16_t reg : regulators)
        regulons[reg];

      int max_subnets = 65536; // arbitrary
      int max_subnets_when_criteria_filled = 0; // needed since multithreading
      bool stoppingCriteriaMet = false;

#pragma omp parallel for num_threads(nthreads)
      for (int subnet_num = 0; subnet_num < max_subnets; ++subnet_num) {
        if (stoppingCriteriaMet && subnet_num >= max_subnets_when_criteria_filled)
          continue; // skip loop iteration if stopping condition is met
        gene_to_floats subsample_exp_mat;
#pragma omp critical(randObjectAccess)
        {
          subsample_exp_mat = sampleExpMatAndReCopulaTransform(
              exp_mat, tot_num_subsample, rand);
        }

        const auto &[subnet, FPR_estimate_subnet] = ARACNe3_subnet(
            subsample_exp_mat, regulators, genes, tot_num_samps,
            tot_num_subsample, subnet_num, prune_alpha, method, alpha,
            prune_MaxEnt, output_dir, subnets_dir, subnets_log_dir);

#pragma omp critical(stopConditionCheck)
        {
          subnets.push_back(subnet);
          FPR_estimates.push_back(FPR_estimate_subnet);
          uint16_t min_regulon_size = 65535U;
          // add any new edges to the regulon_set
          for (const auto [reg, tar_mi] : subnet) {
            for (const auto [tar, mi] : tar_mi)
              regulons[reg].insert(tar);
            if (regulons[reg].size() < min_regulon_size)
              min_regulon_size = regulons[reg].size();
          }
          if (min_regulon_size >= targets_per_regulator) {
            stoppingCriteriaMet = true;
            // adaptive + race conditions can cause overestimation of subnets
            if (subnet_num > max_subnets_when_criteria_filled)
              max_subnets_when_criteria_filled = subnet_num; 
          }
        }
      }
      subnets.resize(max_subnets_when_criteria_filled);
      num_subnets = subnets.size();
    } else if (!adaptive) {
      subnets = std::vector<gene_to_gene_to_float>(num_subnets);
      FPR_estimates = std::vector<float>(num_subnets);
#pragma omp parallel for num_threads(nthreads)
      for (int i = 0; i < num_subnets; ++i) {
        gene_to_floats subsample_exp_mat;
#pragma omp critical(randObjectAccess)
        {
          subsample_exp_mat = sampleExpMatAndReCopulaTransform(
              exp_mat, tot_num_subsample, rand);
        }
        const auto &[subnet, FPR_estimate_subnet] = ARACNe3_subnet(
            subsample_exp_mat, regulators, genes, tot_num_samps,
            tot_num_subsample, i, prune_alpha, method, alpha, prune_MaxEnt,
            output_dir, subnets_dir, subnets_log_dir);
        subnets[i] = subnet;
        FPR_estimates[i] = FPR_estimate_subnet;
      }
    }

    //-------time module-------
    log_output << watch1.getSeconds() << std::endl;
    //-------------------------

    log_output << "Total subnetworks generated: " + std::to_string(num_subnets)
               << std::endl;
  } else if (go_to_consolidate) {

    //-------time module-------
    log_output << "\nConsolidation requested." << std::endl;
    log_output << "Reading subnetwork(s) time: ";
    watch1.reset();
    //-------------------------

    for (uint16_t subnet_num = 0; subnet_num < num_subnets_to_consolidate;
         ++subnet_num) {
      const std::string subnets_file_path = subnets_dir + "subnet" + std::to_string(subnet_num+1) + ".tsv";
      const std::string subnets_log_file_path = subnets_log_dir + "log_subnet" + std::to_string(subnet_num+1) + ".txt";

      try {
        const auto &[subnet, FPR_estimate_subnet] = loadARACNe3SubnetsAndUpdateFPRFromLog(subnets_file_path, subnets_log_file_path);
        subnets.push_back(subnet);
        FPR_estimates.push_back(FPR_estimate_subnet);
      } catch (TooManySubnetsRequested e) {
        std::cout << "Warning: " + std::string(e.what()) << std::endl;
        break;
      }
    }
    num_subnets = subnets.size();

    //-------time module-------
    log_output << watch1.getSeconds() << std::endl;
    //-------------------------

    log_output << "Total subnets read: " + std::to_string(num_subnets)
               << std::endl;
  }

  // set the FPR estimate
  FPR_estimate =
      std::accumulate(FPR_estimates.begin(), FPR_estimates.end(), 0.0f) /
      FPR_estimates.size();

  if (!do_not_consolidate) {

    //-------time module-------
    log_output << "\nConsolidating subnetwork(s) time: ";
    watch1.reset();
    //-------------------------

    std::vector<consolidated_df_row> final_df = consolidate_subnets_vec(subnets, exp_mat, regulators, genes, ranks_mat);

    //-------time module-------
    log_output << watch1.getSeconds() << std::endl;
    log_output << "Subnetworks consolidated: " << std::to_string(num_subnets) << std::endl;
    log_output << "\nWriting final network..." << std::endl;
    //-------------------------

    writeConsolidatedNetwork(final_df, output_dir + "consolidated-net_" + runid + ".tsv");

  } else if (do_not_consolidate) {

    //-------time module-------
    log_output << "\nNo consolidation requested." << std::endl;
    //-------------------------

  }

  using namespace std::string_literals;
  const char *success_A3 =
      R"(

                |
                |
                ;                            
                ;                            
                |                            
           ,  / | \   ,
         , ;_/ ,L-, `_;  ,
         \._/.ARACNe3.\_./
           \_./(::)\._/                      
                ''


SUCCESS!
)";
  std::cout << success_A3 << std::endl;
  log_output << success_A3 << std::endl;

  return EXIT_SUCCESS;
}

#pragma once

#include "logger.hpp"

#include <cstdint>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// ---- Typedefs ----

using gene_id = uint16_t;
using compression_map = std::unordered_map<std::string, gene_id>;
using decompression_map = std::vector<std::string>;

using geneset = std::unordered_set<gene_id>;
using gene_to_geneset = std::unordered_map<gene_id, geneset>;
using gene_to_floats = std::unordered_map<gene_id, std::vector<float>>;
using gene_to_shorts = std::unordered_map<gene_id, std::vector<uint16_t>>;
using gene_to_float = std::unordered_map<gene_id, float>;
using gene_to_gene_to_float = std::unordered_map<gene_id, gene_to_float>;

using vv_float = std::vector<std::vector<float>>;

using pair_string_vecs =
    std::pair<std::vector<std::string>, std::vector<std::string>>;

struct ARACNe3_df {
  const gene_id regulator;
  const gene_id target;
  const float final_mi;
  const float final_scc;
  const uint16_t num_subnets_incident;
  const float final_log_p;

  ARACNe3_df(const gene_id r, const gene_id t, const float mi, const float scc,
             const uint16_t n, const float lp)
      : regulator(r), target(t), num_subnets_incident(n), final_mi(mi),
        final_scc(scc), final_log_p(lp){};
};

// ---- Functions ----

const std::string getSystemTemporaryDirectory();

/**
 * Creates a new full directory path using the Standard Template Library (STL)
 * filesystem functions.
 *
 * @param dir_name A constant reference to a string representing the directory
 * path.
 */
bool makeDirs(const std::string &dir_name, Logger *const logger);

std::tuple<vv_float, geneset, compression_map, decompression_map>
readExpMatrixAndCopulaTransform(const std::string &exp_mat_file_path,
                                std::mt19937 &rnd, Logger *const logger);

/*
 Create a subsampled gene_to_floats.  Requires that exp_mat and
n_subsample are set.
 */
const vv_float sampleExpMatAndReCopulaTransform(const vv_float &exp_mat,
                                                const uint16_t n_subsample,
                                                std::mt19937 &rnd);

geneset readRegList(const std::string &regulators_list_file_path,
                    const compression_map &defined_genes, Logger *const logger,
                    bool verbose);

/*
 Function that prints the Regulator, Target, and MI to the output_dir given the
 output_suffix.  Does not print to the console.  The data structure input is a
 gene_to_edge_tars, which is defined in "ARACNe3.hpp".
 */
void writeNetworkRegTarMI(const std::string &output_file_name, const char sep,
                          const gene_to_gene_to_float &network,
                          const decompression_map &decompressor);

void writeARACNe3DF(const std::string &output_file_name, const char sep,
                    const std::vector<ARACNe3_df> &output_df,
                    const decompression_map &decompressor);

/*
 * Lists the files in the provided subnets dir and matches them to their log
 * files, assuming heterogeneous runids.  Returns a pair of string vectors that
 * are index-matched according to subnet, so the correct log statistics are used
 * in later computation
 */
pair_string_vecs
findSubnetFilesAndSubnetLogFiles(const std::string &subnets_dir,
                                 const std::string &subnets_log_dir);

/*
 Reads a subnet file and then updates the FPR_estimates vector defined in
 "subnet_operations.cpp"
 */
std::pair<gene_to_gene_to_float, float>
loadARACNe3SubnetsAndUpdateFPRFromLog(const std::string &subnet_file_path,
                                      const std::string &subnet_log_file_path,
                                      const compression_map &defined_genes,
                                      const geneset &regulators,
                                      Logger *const logger);

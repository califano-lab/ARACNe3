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

/**
 * Wraps STL functions to create a full directory path; this is mainly a
 * convenience function
 *
 * @param dir_name A constant reference to a string representing the directory
 * path.
 */
bool makeDirs(const std::string &dir_name, Logger *const logger);

/**
 * @brief Read an expression matrix file, copula-transform the data, and return
 * the parsed data.
 *
 * This function reads an expression matrix file, parses its contents, and
 * copula-transforms the expression values. It expects the input file to be in a
 * tab-separated format, where the first row contains the sample names, and each
 * subsequent row represents a gene, with the first column containing the gene
 * name and the remaining columns containing the expression values for each
 * sample.
 *
 * The function performs the following steps:
 * 1. Reads the file line by line, skipping the first row (header).
 * 2. For each line:
 *    - Extracts the gene name from the first column.
 *    - Parses the expression values from the remaining columns.
 *    - Checks if the number of samples is consistent with the previous rows.
 *    - Copula-transforms the expression values.
 *    - Stores the gene in the decompression map and the gene index in the
 * compression map.
 *    - Adds the copula-transformed expression values to the expression matrix.
 * 3. Returns a tuple containing the copula-transformed expression matrix, all
 * genes, the compression mapping, and the decompression maping.
 *
 * @param exp_mat_file_path The path to the expression matrix file.
 * @param rnd A reference to a std::mt19937 random number engine used for copula
 * transformation.
 * @param logger A pointer to a Logger object for logging errors. If nullptr, no
 * logging is performed.
 *
 * @return A tuple containing:
 *         - vv_float: The copula-transformed expression matrix.
 *         - geneset: A set of gene IDs.
 *         - compression_map: A map from gene names to uint16_t indices.
 *         - decompression_map: A vector of gene names, where the index
 * corresponds to the uint16_t index in the compression map.
 *
 * @note If the expression matrix file cannot be opened or if the genes in the
 * matrix do not have equal samples, the function terminates the program with an
 * error message.
 */
std::tuple<vv_float, geneset, compression_map, decompression_map>
readExpMatrixAndCopulaTransform(const std::string &exp_mat_file_path,
                                std::mt19937 &rnd, Logger *const logger);

/**
 * @brief Sample a subset of the expression matrix, re-copula-transform the
 * subsampled data, and return the result.
 *
 * This function takes an expression matrix, randomly samples a subset of its
 * columns (samples), re-copula-transforms the subsampled data, and returns the
 * resulting subsampled expression matrix.
 *
 * @param exp_mat The input expression matrix to be subsampled.
 * @param n_subsamp The number of samples to be subsampled from the expression
 * matrix.
 * @param rnd A reference to a std::mt19937 random number engine used for
 * sampling and copula transformation.
 *
 * @return The subsampled and re-copula-transformed expression matrix.
 *
 * @note The subsampled expression matrix has the same number of rows (genes) as
 * the input matrix, but the number of columns (samples) is reduced to
 * n_subsamp.
 */
const vv_float sampleExpMatAndReCopulaTransform(const vv_float &exp_mat,
                                                const uint16_t n_subsample,
                                                std::mt19937 &rnd);

/**
 * @brief Read a list of regulators from a file and filter out genes not defined
 * in the expression matrix.
 *
 * This function reads a file containing regulator genes (one per line) and
 * filters out genes not present in the provided compression map
 * (defined_genes). It returns the set of those gene IDs (compressed) which are
 * regulators.
 *
 * @param regulators_list_file_path The path to the file containing the list of
 * regulator genes.
 * @param defined_genes A compression map containing the mapping of gene names
 * to their IDs in the expression matrix.
 * @param logger A pointer to a Logger object for logging messages. If nullptr,
 * no logging is performed.
 * @param verbose If true, all warning messages are printed; otherwise, messages
 * are suppressed after the first 3.
 *
 * @return A set of gene IDs representing the regulators found in both the
 * regulator list file and the defined genes.
 *
 * @note The function assumes that the file contains one regulator gene per
 * line. If the file cannot be opened or no regulators are found after
 * filtering, the function terminates the program.
 */
geneset readRegList(const std::string &regulators_list_file_path,
                    const compression_map &defined_genes, Logger *const logger,
                    bool verbose);

/**
 * @brief Write a network to a file, with regulator-target pairs and their
 * mutual information (MI) values.
 *
 * @param output_file_name The name of the output file to write the network to.
 * @param sep The separator character used to delimit the fields in the output
 * file.
 * @param network A map of maps representing the gene regulatory network, where
 * the keys are regulator gene IDs and the values are regulons, represented as
 * individual maps of target gene IDs to their MI values.
 * @param decompressor A decompression map that maps gene IDs to their
 * corresponding gene names.
 */
void writeNetworkRegTarMI(const std::string &output_file_name, const char sep,
                          const gene_to_gene_to_float &network,
                          const decompression_map &decompressor);

/**
 * @brief Write ARACNe3 data frame to a file, with regulator-target pairs and
 * their associated metrics.
 *
 * @param output_file_name The name of the output file to write the ARACNe3 data
 * frame to.
 * @param sep The separator character used to delimit the fields in the output
 * file.
 * @param output_df A vector of ARACNe3_df structs, representing the ARACNe3
 * data frame to be written.
 * @param decompressor A decompression map that maps gene IDs to their
 * corresponding gene names.
 */
void writeARACNe3DF(const std::string &output_file_name, const char sep,
                    const std::vector<ARACNe3_df> &output_df,
                    const decompression_map &decompressor);

/**
 * @brief Find subnet files and their corresponding log files in the specified
 * directories.
 *
 * This function searches for subnet files in the specified `subnets_dir`
 * directory and attempts to find their corresponding log files in the
 * `subnets_log_dir` directory. It assumes that the subnet files have a `.tsv`
 * extension and the log files have a `.txt` extension with a prefix of "log-"
 * followed by the subnet file name.
 *
 * @param subnets_dir The directory path containing the subnet files.
 * @param subnets_log_dir The directory path containing the subnet log files.
 * @param logger A pointer to a Logger object for logging messages. If nullptr,
 * no logging is performed.
 *
 * @return A pair of vectors containing the matched subnet file names and their
 * corresponding log file names.
 *
 * @note The function expects the subnet files to have a `.tsv` extension and
 * the log files to have a `.txt` extension with a prefix of "log-" followed by
 * the subnet file name.
 */
pair_string_vecs
findSubnetFilesAndSubnetLogFiles(const std::string &subnets_dir,
                                 const std::string &subnets_log_dir,
                                 Logger *const logger);
/**
 * @brief Load ARACNe3 subnets from files and update the false positive rate
 * (FPR) estimate.
 *
 * This function reads an ARACNe3 subnet file and its corresponding log file to
 * load the gene regulatory network and update the FPR estimate. The subnet and
 * log file should follow the exact conventions of the outputs from
 * `writeNetworkRegTarMI` and the subnet logging routine defined in
 * `subnet_operations.cpp`. In particular, the subnet file contains
 * regulator-target pairs and their mutual information (MI) values. The log file
 * contains information about the pruning method, alpha value, and the number of
 * edges after threshold and MaxEnt pruning, in a format that is output by the
 * function `createARACNe3Subnet`.
 *
 * @param subnet_file_path The path to the ARACNe3 subnet file.
 * @param subnet_log_file_path The path to the corresponding log file for the
 * subnet.
 * @param defined_genes A compression map containing the mapping of gene names
 * to their IDs.
 * @param regulators A set of regulator gene IDs.
 * @param logger A pointer to a Logger object for logging messages. If nullptr,
 * no logging is performed.
 *
 * @return A pair containing the loaded gene regulatory network and the updated
 * FPR estimate.
 *
 * @note This function is only intended to run as a consolidation of previously
 * generated ARACNe3 subnetworks and subnetwork log files. This is usually
 * performed when individual subnetworks take hour(s) to generate, and the user
 * intends to offload jobs onto an HPC cluster with a separate consolidation
 * step to follow. If the subnet files contain gene names not defined in the
 * expression matrix, the function terminates the program with an error message.
 */
std::pair<gene_to_gene_to_float, float>
loadARACNe3SubnetsAndUpdateFPRFromLog(const std::string &subnet_file_path,
                                      const std::string &subnet_log_file_path,
                                      const compression_map &defined_genes,
                                      const geneset &regulators,
                                      Logger *const logger);

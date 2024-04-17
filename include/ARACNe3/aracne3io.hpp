#pragma once

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

using str_vec = std::vector<std::string>;
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

class ARACNe3IOHandler {
public:
  /**
   * @brief Read an expression matrix, copula-transform the data, and
   * return the parsed data.
   *
   * @param rnd A reference to a std::mt19937 random number engine used for
   * copula transform.
   *
   * @return A tuple containing:
   *         - vv_float: The copula-transformed expression matrix.
   *         - geneset: A set of gene IDs.
   *         - compression_map: A map from gene names to uint16_t indices.
   *         - decompression_map: A vector of gene names, where the index
   * corresponds to the uint16_t index in the compression map.
   */
  virtual std::tuple<vv_float, geneset, compression_map, decompression_map>
  readExpMatrixAndCopulaTransform(std::mt19937 &rnd) const = 0;

  /**
   * @brief Filters regulators based on expression data, returning valid gene
   * IDs and warnings.
   *
   * @param defined_genes Map of gene names to IDs, set by
   * readExpMatrixAndCopulaTransform.
   *
   * @return Pair of gene ID set from regulators with expression data, and
   * warning messages.
   */
  virtual std::pair<geneset, str_vec>
  readRegList(const compression_map &defined_genes) const = 0;

  /**
   * @brief Write ARACNe3 data frame as regulator-target pairs and
   * their associated metrics.
   *
   * @param output_df A vector of ARACNe3_df structs, representing the ARACNe3
   * data frame to be written.
   * @param decompressor A decompression map that maps gene IDs to their
   * corresponding gene names.
   */
  virtual void writeARACNe3DF(const std::vector<ARACNe3_df> &output_df,
                              const decompression_map &decompressor) const = 0;

  /**
   * @brief Write a network, with regulator-target pairs and their
   * mutual information (MI) values.
   *
   * @param output_file_name The name of the output file to write the network
   * to.
   * @param sep The separator character used to delimit the fields in the output
   * file.
   * @param network A map of maps representing the gene regulatory network,
   * where the keys are regulator gene IDs and the values are regulons,
   * represented as individual maps of target gene IDs to their MI values.
   * @param decompressor A decompression map that maps gene IDs to their
   * corresponding gene names.
   *
   * @note This functionality is currently only available for filesystem
   * implementations.
   */
  virtual void
  writeNetworkRegTarMI(const uint16_t subnet_number,
                       const gene_to_gene_to_float &network,
                       const decompression_map &decompressor) const {
    throw std::runtime_error(
        "Single subnetwork inference not supported in this implementation.");
  };

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
   *
   * @return A pair of vectors containing the matched subnet file names and
   * their corresponding log file names.
   *
   * @note The function expects the subnet files to have a `.tsv` extension and
   * the log files to have a `.txt` extension with a prefix of "log-" followed
   * by the subnet file name.
   *
   * @note This functionality is currently only available for filesystem
   * implementations.
   */
  virtual pair_string_vecs
  findSubnetFilesAndSubnetLogFiles(const std::string &subnets_dir,
                                   const std::string &subnets_log_dir) const {
    throw std::runtime_error(
        "Consolidate mode not supported in this implementation.");
  };

  /**
   * @brief Load ARACNe3 subnets from files and update the false positive rate
   * (FPR) estimate.
   *
   * This function reads an ARACNe3 subnet file and its corresponding log file
   * to load the gene regulatory network and update the FPR estimate. The subnet
   * and log file should follow the exact conventions of the outputs from
   * `writeNetworkRegTarMI` and the subnet logging routine defined in
   * `subnet_operations.cpp`. In particular, the subnet file contains
   * regulator-target pairs and their mutual information (MI) values. The log
   * file contains information about the pruning method, alpha value, and the
   * number of edges after threshold and MaxEnt pruning, in a format that is
   * output by the function `createARACNe3Subnet`.
   *
   * @param subnet_file_path The path to the ARACNe3 subnet file.
   * @param subnet_log_file_path The path to the corresponding log file for the
   * subnet.
   * @param defined_genes A compression map containing the mapping of gene names
   * to their IDs.
   * @param regulators A set of regulator gene IDs.
   *
   * @return A pair containing the loaded gene regulatory network and the
   * updated FPR estimate.
   *
   * @note This function is only intended to run as a consolidation of
   * previously generated ARACNe3 subnetworks and subnetwork log files. This is
   * usually performed when individual subnetworks take hour(s) to generate, and
   * the user intends to offload jobs onto an HPC cluster with a separate
   * consolidation step to follow. If the subnet files contain gene names not
   * defined in the expression matrix, the function terminates the program with
   * an error message.
   *
   * @note This functionality is currently only available for filesystem
   * implementations.
   */
  virtual std::pair<gene_to_gene_to_float, float>
  loadARACNe3SubnetsAndUpdateFPRFromLog(const std::string &subnet_file_path,
                                        const std::string &subnet_log_file_path,
                                        const compression_map &defined_genes,
                                        const geneset &regulators) const {
    throw std::runtime_error(
        "Consolidate mode not supported in this implementation.");
  };
};

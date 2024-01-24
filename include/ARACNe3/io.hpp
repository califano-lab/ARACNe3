#pragma once

#include "logger.hpp"
#include "ARACNe3.hpp"

#include <random>
#include <string>

#if defined __linux__ || defined __APPLE__
const std::string hiddenfpre = ".";
const char directory_slash = '/';
#elif defined _WIN32 // in windows you must set a hidden file via properties
const char directory_slash = '\\';
const std::string hiddenfpre = "";
#else
const std::string hiddenfpre = "";
#endif /* __linux__ || __APPLE__ */

typedef struct consolidated_df_row {
  const gene_id regulator;
  const gene_id target;
  const float final_mi;
  const float final_scc;
  const uint16_t num_subnets_incident;
  const double final_log_p;
  consolidated_df_row(const gene_id &r, const gene_id &t, const float &mi,
                      const float &scc, const uint16_t &n, const double &lp)
      : regulator(r), target(t), num_subnets_incident(n), final_mi(mi),
        final_scc(scc), final_log_p(lp){};
} consolidated_df_row;

typedef std::pair<std::vector<std::string>, std::vector<std::string>>
    pair_string_vecs;

/**
 * Transforms a Unix-style directory name into a universal format suitable for
 * both Unix and Windows. For Windows, it replaces '/' with '\'. Intended for
 * display purposes, adapting to the system-specific directory separator.
 *
 * @param dir_name A string representing the directory name to be transformed.
 * @return A string with the universal directory name.
 */
std::string makeUnixDirectoryNameUniversal(std::string dir_name);

/**
 * Creates a new full directory path using the Standard Template Library (STL)
 * filesystem functions.
 *
 * @param dir_name A constant reference to a string representing the directory
 * path.
 */
bool makeDirs(const std::string &dir_name, Logger *const logger);

std::tuple<const gene_to_floats, const gene_to_shorts, const geneset,
           const uint16_t>
readExpMatrixAndCopulaTransform(const std::string &filename,
                                std::mt19937 &rand);
const geneset readRegList(const std::string &filename, const bool verbose);
gene_to_floats
sampleExpMatAndReCopulaTransform(const gene_to_floats &exp_mat,
                                 const uint16_t &tot_num_subsample,
                                 std::mt19937 &rand);

void writeNetworkRegTarMI(const gene_to_gene_to_float& network,
                          const std::string &file_path);

void writeConsolidatedNetwork(const std::vector<consolidated_df_row> &final_df,
                              const std::string &file_path);

pair_string_vecs
findSubnetFilesAndSubnetLogFiles(const std::string &subnets_dir,
                                 const std::string &subnets_log_dir);

std::pair<gene_to_gene_to_float, float>
loadARACNe3SubnetsAndUpdateFPRFromLog(const std::string &subnet_file_path,
                                      const std::string &subnet_log_file_path);

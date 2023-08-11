#pragma once

#include "ARACNe3.hpp"
#include <random>
#include <string>

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

std::string makeUnixDirectoryNameUniversal(std::string &dir_name);
std::string makeUnixDirectoryNameUniversal(std::string &&dir_name);
void makeDir(const std::string &dir_name);

std::tuple<const gene_to_floats, const gene_to_shorts, const geneset,
           const uint16_t>
readExpMatrixAndCopulaTransform(const std::string &filename,
                                const float &subsampling_percent,
                                std::mt19937 &rand);
const geneset readRegList(const std::string &filename);
gene_to_floats
sampleExpMatAndReCopulaTransform(const gene_to_floats &exp_mat,
                                 const uint16_t &tot_num_subsample,
                                 std::mt19937 &rand);

void writeNetworkRegTarMI(gene_to_gene_to_float &network,
                          const std::string &file_path);

void writeConsolidatedNetwork(const std::vector<consolidated_df_row> &final_df,
                              const std::string &file_path);

std::pair<gene_to_gene_to_float, float>
loadARACNe3SubnetsAndUpdateFPRFromLog(const std::string &subnet_file_path,
                                      const std::string &subnet_log_file_path);

#pragma once

#include <vector>
#include "ARACNe3.hpp"
#include "io.hpp"

gene_to_edge_tars ARACNe3_subnet(gene_to_floats subnet_matrix, const uint16_t subnet_num,
                       const bool prune_alpha, const std::string &method,
                       const float alpha, const bool prune_MaxEnt,
                       const std::string &output_dir,
                       const std::string &subnets_dir,
                       const std::string &log_dir, const uint16_t nthreads);
std::vector<consolidated_df_row> consolidate_subnets_vec(std::vector<gene_to_edge_tars> &subnets);

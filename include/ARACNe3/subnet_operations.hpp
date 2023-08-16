#pragma once

#include "ARACNe3.hpp"
#include "apmi_nullmodel.hpp"
#include "io.hpp"
#include <vector>

std::pair<gene_to_gene_to_float, float> createARACNe3Subnet(
    const gene_to_floats &subsample_exp_mat, const geneset &regulators,
    const geneset &genes, const uint16_t tot_num_samps,
    const uint16_t tot_num_subsample, const uint16_t cur_subnet_ct,
    const bool prune_alpha, const APMINullModel &nullmodel,
    const std::string &method, const float alpha, const bool prune_MaxEnt,
    const std::string &output_dir, const std::string &subnets_dir,
    const std::string &subnet_log_dir);

const std::vector<consolidated_df_row>
consolidateSubnetsVec(const std::vector<gene_to_gene_to_float> &subnets,
                      const float FPR_estimate, const gene_to_floats &exp_mat,
                      const geneset &regulators, const geneset &genes,
                      const gene_to_shorts &ranks_mat);

class TooManySubnetsRequested : public std::exception {
public:
  const char *what() {
    return "You tried to consolidate too many subnetworks. Subnet files and "
           "subnet log files must exist to all subnetworks specified, and "
           "must be organized exactly according to ARACNe3 output.";
  }
};

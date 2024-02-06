#pragma once

#include "apmi_nullmodel.hpp"
#include "aracne3io.hpp"

#include <vector>

std::tuple<gene_to_gene_to_float, float, uint32_t> createARACNe3Subnet(
    const vv_float &subsample_exp_mat, const geneset &regulators,
    const geneset &genes, const uint16_t tot_num_samps,
    const uint16_t tot_num_subsample, const uint16_t cur_subnet_ct,
    const bool prune_alpha, const APMINullModel &nullmodel,
    const std::string &method, const float alpha, const bool prune_MaxEnt,
    const std::string &output_dir, const std::string &subnets_dir,
    const std::string &subnet_log_dir, const uint16_t nthreads,
    const std::string &runid, const decompression_map &decompressor,
    const bool save_subnet);

const std::vector<ARACNe3_df>
consolidateSubnetsVec(const std::vector<gene_to_gene_to_float> &subnets,
                      const float FPR_estimate, const vv_float &exp_mat,
                      const geneset &regulators, const geneset &genes,
                      std::mt19937 &rnd);

class TooManySubnetsRequested : public std::exception {
public:
  const char *what() {
    return "You tried to consolidate too many subnetworks. Subnet files and "
           "subnet log files must exist to all subnetworks specified, and "
           "must be organized exactly according to ARACNe3 output.";
  }
};

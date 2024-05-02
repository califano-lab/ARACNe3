#pragma once

#include "apmi_nullmodel.hpp"
#include "aracne3io.hpp"

#include <cstdint>
#include <stdexcept>
#include <vector>

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

std::tuple<gene_to_gene_to_float, float, uint32_t> createARACNe3Subnet(
    const vv_float &subsample_exp_mat, const geneset &regulators,
    const geneset &genes, const uint16_t tot_num_samps,
    const uint16_t tot_num_subsample, const uint16_t subnet_number,
    const bool prune_alpha, const APMINullModel &nullmodel,
    const std::string &method, const float alpha, const bool prune_MaxEnt,
    const std::string &subnets_log_dir, const uint8_t threads,
    const std::string &runid, const decompression_map &decompressor,
    const bool save_subnet, const ARACNe3IOHandler &io);

const std::vector<ARACNe3_df>
consolidateSubnetsVec(const std::vector<gene_to_gene_to_float> &subnets,
                      const float FPR_estimate, const vv_float &exp_mat,
                      const geneset &regulators, const geneset &genes);

class TooManySubnetsRequestedException : public std::runtime_error {
public:
  explicit TooManySubnetsRequestedException()
      : std::runtime_error(
            "You tried to consolidate too many subnetworks. Subnet files and "
            "subnet log files must exist to all subnetworks specified, and "
            "must be organized exactly according to ARACNe3 output.") {}
};

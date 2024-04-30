#pragma once

#include "aracne3io.hpp"
#include "logger.hpp"

class SubnetLogger : public Logger {
public:
  /**
   * Constructor that opens a log file and writes the time.
   *
   * @param log_file_path Path to the log file.
   */
  SubnetLogger(const std::string &log_file_name);

  /**
   * Initializes and logs subnetwork details. Designed to be called
   * post-construction for modular design.
   *
   * @param runid Identifier of the analysis.
   * @param subnet_number The subnet number.
   * @param regulators Regulator genes.
   * @param genes Target genes.
   * @param tot_num_samps Total number of samples in the dataset.
   * @param tot_num_subsample Number of samples after subsampling.
   * @param method Method of error control for first pruning step
   * (FDR;FPR;FWER).
   * @param alpha Alpha parameter for control by the pruning method.
   * @param prune_MaxEnt Indicates if MaxEnt pruning was performed.
   */
  void initSubnetLog(const std::string &runid, const uint16_t subnet_number,
                     const geneset &regulators, const geneset &genes,
                     const uint16_t tot_num_samps,
                     const uint16_t tot_num_subsample,
                     const std::string &method, const float alpha,
                     const bool prune_MaxEnt);
};

#include "subnet_logger.hpp"

SubnetLogger::SubnetLogger(const std::string &log_file_name)
    : Logger(log_file_name) {}

void SubnetLogger::initSubnetLog(
    const std::string &runid, const uint16_t subnet_number,
    const geneset &regulators, const geneset &genes,
    const uint16_t tot_num_samps, const uint16_t tot_num_subsample,
    const std::string &method, const float alpha, const bool prune_MaxEnt) {
  writeLine("\nRun ID: " + runid);
  writeLine("Subnetwork #: " + std::to_string(subnet_number));
  writeLine("Total # regulators (defined in gexp mat): " +
            std::to_string(regulators.size()));
  writeLine("Total # targets: " + std::to_string(genes.size()));
  writeLine("Total # samples: " + std::to_string(tot_num_samps));
  writeLine("Subsampled quantity: " + std::to_string(tot_num_subsample));
  writeLine("Total possible edges: " +
            std::to_string(regulators.size() * (genes.size() - 1)));
  writeLine("Method of first pruning step: " + method);
  writeLine("Alpha: " + std::to_string(alpha));
  writeLine("MaxEnt Pruning: " + std::string(prune_MaxEnt ? "true" : "false"));
  writeLine("\n-----------Begin Network Generation-----------");
}

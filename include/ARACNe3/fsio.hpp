#pragma once

#include "aracne3io.hpp"

#include <iostream>

class FilesystemIOHandler : public ARACNe3IOHandler {
public:
  FilesystemIOHandler(const std::string &exp_mat_file_path,
                      const std::string &regulators_list_file_path,
                      const std::string &final_output_file_name,
                      const std::string &subnets_dir, const std::string &runid,
                      const char sep);

  std::tuple<vv_float, geneset, compression_map, decompression_map>
  readExpMatrixAndCopulaTransform(std::mt19937 &rnd) const override;
  geneset readRegList(const compression_map &defined_genes,
                      Logger *const logger, bool verbose) const override;
  void writeARACNe3DF(const std::vector<ARACNe3_df> &output_df,
                      const decompression_map &decompressor) const override;
  void
  writeNetworkRegTarMI(const uint16_t subnet_number,
                       const gene_to_gene_to_float &network,
                       const decompression_map &decompressor) const override;
  pair_string_vecs findSubnetFilesAndSubnetLogFiles(
      const std::string &subnets_dir,
      const std::string &subnets_log_dir) const override;
  std::pair<gene_to_gene_to_float, float> loadARACNe3SubnetsAndUpdateFPRFromLog(
      const std::string &subnet_file_path,
      const std::string &subnet_log_file_path,
      const compression_map &defined_genes,
      const geneset &regulators) const override;

private:
  const std::string exp_mat_file_path, regulators_list_file_path,
      final_output_file_name, subnets_dir, runid;
  const char sep;
};

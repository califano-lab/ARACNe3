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
  readExpMatrixAndCopulaTransform(std::mt19937 &rnd,
                                  Logger *const logger) const override;
  geneset readRegList(const compression_map &defined_genes,
                      Logger *const logger, bool verbose) const override;
  void writeARACNe3DF(const std::vector<ARACNe3_df> &output_df,
                      const decompression_map &decompressor) const override;
  void
  writeNetworkRegTarMI(const uint16_t subnet_number,
                       const gene_to_gene_to_float &network,
                       const decompression_map &decompressor) const override;
  pair_string_vecs
  findSubnetFilesAndSubnetLogFiles(const std::string &subnets_dir,
                                   const std::string &subnets_log_dir,
                                   Logger *const logger) const override;
  std::pair<gene_to_gene_to_float, float>
  loadARACNe3SubnetsAndUpdateFPRFromLog(const std::string &subnet_file_path,
                                        const std::string &subnet_log_file_path,
                                        const compression_map &defined_genes,
                                        const geneset &regulators,
                                        Logger *const logger) const override;

private:
  const std::string exp_mat_file_path, regulators_list_file_path,
      final_output_file_name, subnets_dir, runid;
  const char sep;

  /**
   * @brief Check if a file stream is open, and exit the program if it's not.
   *
   * This function checks if the provided input file stream (std::ifstream) is
   * open. If the file stream is not open, it prints an error message to the
   * standard error stream (std::cerr) and logs the error message using the
   * provided Logger object (if not nullptr). Then, it terminates the program
   * with a failure status.
   *
   * @param ifs The input file stream to check.
   * @param file_path The path of the file associated with the file stream.
   * @param logger A pointer to a Logger object for logging the error message.
   *               If nullptr, no logging is performed.
   *
   * @return void
   *
   * @note This function does not return if the file stream is not open. It
   *       terminates the program using std::exit(EXIT_FAILURE).
   */
  void exitIfFileNotOpen(std::ifstream &ifs, const std::string &file_path,
                         Logger *const logger) const {
    if (!ifs.is_open()) {
      const std::string err_msg =
          "Error: could not open file \"" + file_path + "\".";

      std::cerr << err_msg << std::endl;
      if (logger)
        logger->writeLineWithTime(err_msg);

      std::exit(EXIT_FAILURE);
    }

    return;
  }
};

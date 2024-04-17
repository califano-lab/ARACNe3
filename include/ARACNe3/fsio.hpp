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
  std::pair<geneset, str_vec>
  readRegList(const compression_map &defined_genes) const override;
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

class FileNotOpenException : public std::runtime_error {
public:
  explicit FileNotOpenException(const std::string &filename)
      : std::runtime_error("Could not open file: " + filename) {}
};

class DirectoryReadException : public std::runtime_error {
public:
  explicit DirectoryReadException(const std::string &directoryPath)
      : std::runtime_error("Error reading directory: " + directoryPath) {}
};

class NoRegulatorsException : public std::runtime_error {
public:
  explicit NoRegulatorsException(const std::string &filename)
      : std::runtime_error("No regulators found in file: " + filename) {}
};

class MatrixDataException : public std::runtime_error {
public:
  explicit MatrixDataException(const std::string &message)
      : std::runtime_error(message) {}
};

class FileNotFoundException : public std::runtime_error {
public:
  explicit FileNotFoundException(const std::string &filename,
                                 const std::string &msg = "")
      : std::runtime_error("File not found: " + filename + ". " + msg) {}
};

class GeneInSubnetNotInMatrixException : public std::runtime_error {
public:
  explicit GeneInSubnetNotInMatrixException(const std::string &gene)
      : std::runtime_error(
            "Subnetwork file cannot contain genes gene names undefined in the "
            "expression matrix. You should only consolidate subnetworks from "
            "the expression matrix used to generate them. (" +
            gene + ")") {}
};

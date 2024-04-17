#include "fsio.hpp"
#include "algorithms.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>

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

std::ifstream getReadStream(const std::string &filePath) {
  std::ifstream ifs{filePath};
  if (!ifs.is_open())
    throw FileNotOpenException(filePath);
  return ifs;
}

std::ofstream getWriteStream(const std::string &filePath) {
  std::ofstream ofs{filePath};
  if (!ofs.is_open())
    throw FileNotOpenException(filePath);
  return ofs;
}

FilesystemIOHandler::FilesystemIOHandler(const std::string &emfp,
                                         const std::string &rlfp,
                                         const std::string &fofn,
                                         const std::string &sd,
                                         const std::string &r, const char s)
    : exp_mat_file_path(emfp), regulators_list_file_path(rlfp),
      final_output_file_name(fofn), subnets_dir(sd), runid(r), sep(s){};

std::tuple<vv_float, geneset, compression_map, decompression_map>
FilesystemIOHandler::readExpMatrixAndCopulaTransform(std::mt19937 &rnd) const {

  auto ifs = getReadStream(exp_mat_file_path);

  // pop the first line
  std::string cur_line;
  getline(ifs, cur_line, '\n');
  if (cur_line.back() == '\r')
    cur_line.pop_back();

  // parse the exp_mat
  vv_float gexp_matrix;
  geneset genes;
  compression_map compressor;
  decompression_map decompressor;

  while (std::getline(ifs, cur_line, '\n')) {
    if (cur_line.back() == '\r')
      cur_line.pop_back();
    std::vector<float> expr_vec;

    if (gexp_matrix.size() > 0u)
      expr_vec.reserve(gexp_matrix.at(0).size());

    // gene name is first column of matrix
    std::size_t cur_pos = 0u, cur_end = cur_line.find_first_of("\t", cur_pos);
    const std::string gene = cur_line.substr(cur_pos, cur_end - cur_pos);
    cur_pos = cur_end + 1;

    // add gene to all 3 data structures
    decompressor.push_back(gene);
    compressor.insert({gene, decompressor.size() - 1u});
    genes.insert(compressor.at(gene));

    // fill gene vector
    while ((cur_end = cur_line.find_first_of("\t", cur_pos)) !=
           std::string::npos) {
      if (cur_end > cur_pos)
        expr_vec.emplace_back(
            stof(cur_line.substr(cur_pos, cur_end - cur_pos)));
      cur_pos = cur_end + 1;
    }
    expr_vec.emplace_back(stof(cur_line.substr(cur_pos, std::string::npos)));

    // if this vector is bigger than the first, the data are invalid
    if (gexp_matrix.size() > 0u && expr_vec.size() != gexp_matrix[0].size()) {
      std::string err_msg = "Error: genes in expression matrix do not "
                            "all have equal samples. Check that \"" +
                            gene + "\" has same number of samples as " +
                            decompressor[0];

      throw MatrixDataException(err_msg);
    }

    // copula-transform expr_vec
    expr_vec = copulaTransform(expr_vec, rnd);

    gexp_matrix.push_back(expr_vec);
  }

  return {gexp_matrix, genes, compressor, decompressor};
}

std::pair<geneset, str_vec>
FilesystemIOHandler::readRegList(const compression_map &defined_genes) const {

  geneset regulators;
  str_vec warn_list;

  auto ifs = getReadStream(regulators_list_file_path);

  std::string cur_reg;
  while (std::getline(ifs, cur_reg, '\n')) {
    if (cur_reg.back() == '\r')
      cur_reg.pop_back();

    if (defined_genes.find(cur_reg) == defined_genes.end())
      // Warn if not found in exp_mat
      warn_list.push_back(std::string(
          "Warning: \"" + cur_reg +
          "\" found in regulators, but not in expression matrix. Removing "
          "from analysis."));
    else
      regulators.insert(defined_genes.at(cur_reg));
  }

  if (regulators.empty())
    throw NoRegulatorsException(regulators_list_file_path);

  return {regulators, warn_list};
}

void FilesystemIOHandler::writeNetworkRegTarMI(
    const uint16_t subnet_number, const gene_to_gene_to_float &network,
    const decompression_map &decompressor) const {

  const std::string output_file_name = subnets_dir + "subnetwork-" +
                                       std::to_string(subnet_number) + "_" +
                                       runid + ".tsv";

  std::ofstream ofs = getWriteStream(output_file_name);

  ofs << "regulator.values" << sep << "target.values" << sep << "mi.values"
      << std::endl;
  for (const auto &[reg, regulon] : network)
    for (const auto [tar, mi] : regulon)
      ofs << decompressor[reg] << sep << decompressor[tar] << sep << mi << '\n';

  ofs.close();
}

void FilesystemIOHandler::writeARACNe3DF(
    const std::vector<ARACNe3_df> &final_output_df,
    const decompression_map &decompressor) const {

  std::ofstream ofs(final_output_file_name);

  ofs << "regulator.values" << sep << "target.values" << sep << "mi.values"
      << sep << "scc.values" << sep << "count.values" << sep << "log.p.values"
      << '\n';

  for (std::size_t i = 0u; i < final_output_df.size(); ++i) {
    const ARACNe3_df &row = final_output_df[i];
    ofs << decompressor.at(row.regulator) << sep << decompressor.at(row.target)
        << sep << row.final_mi << sep << row.final_scc << sep
        << row.num_subnets_incident << sep << row.final_log_p << '\n';
  }

  ofs.close();
}

pair_string_vecs FilesystemIOHandler::findSubnetFilesAndSubnetLogFiles(
    const std::string &subnets_dir, const std::string &subnets_log_dir) const {
  std::vector<std::string> subnet_filenames, subnet_log_filenames;
  try {
    for (const auto &entry : std::filesystem::directory_iterator(subnets_dir)) {
      if (entry.is_regular_file()) {
        std::string subnet_filename = entry.path().filename().string();
        subnet_filenames.push_back(subnet_filename);

        // Construct the expected log file name, replacing .tsv with .txt
        std::string subnet_log_filename = "log-" + subnet_filename;
        subnet_log_filename.replace(subnet_log_filename.size() - 4, 4, ".txt");

        const std::string full_log_path = subnets_log_dir + subnet_log_filename;
        if (std::filesystem::exists(full_log_path))
          subnet_log_filenames.push_back(subnet_log_filename);
        else {
          const std::string err_msg = "Accompanying log file required for \"" +
                                      subnets_dir + subnet_filename + "\".";

          throw FileNotFoundException(full_log_path, err_msg);
        }
      }
    }
  } catch (std::filesystem::filesystem_error &e) {
    const std::string err_msg =
        "Error reading directory: " + std::string(e.what()) +
        "\nCheck that your output directory has the subdirectories "
        "\"subnets/\" and \"subnets_log/\"";

    throw DirectoryReadException(subnets_dir);
  }
  return {subnet_filenames, subnet_log_filenames};
}

std::pair<gene_to_gene_to_float, float>
FilesystemIOHandler::loadARACNe3SubnetsAndUpdateFPRFromLog(
    const std::string &subnet_file_path,
    const std::string &subnet_log_file_path,
    const compression_map &defined_genes, const geneset &regulators) const {

  // read in the subnet file
  auto subnet_ifs = getReadStream(subnet_file_path);

  // discard the first line (header)
  std::string line;
  std::getline(subnet_ifs, line, '\n');
  if (line.back() == '\r')
    line.pop_back();

  gene_to_gene_to_float subnet;
  while (std::getline(subnet_ifs, line, '\n')) {
    if (line.back() == '\r')
      line.pop_back();

    std::size_t prev = 0u, pos = line.find_first_of("\t", prev);
    const std::string reg = line.substr(prev, pos - prev);
    prev = pos + 1u;

    pos = line.find_first_of("\t", prev);
    const std::string tar = line.substr(prev, pos - prev);
    prev = pos + 1u;

    for (const std::string &gene : {reg, tar})
      if (defined_genes.find(gene) == defined_genes.end())
        throw GeneInSubnetNotInMatrixException(gene);

    const float mi = std::stof(line.substr(prev, std::string::npos));

    subnet[defined_genes.at(reg)][defined_genes.at(tar)] = mi;
  }

  // read in the log file
  auto log_ifs = getReadStream(subnet_log_file_path);

  // discard 9 lines
  std::string discard;
  for (uint8_t l = 0u; l < 9u; ++l)
    std::getline(log_ifs, discard, '\n');

  // next line contains the method
  std::string method;
  std::getline(log_ifs, line, '\n');
  if (line.back() == '\r')
    line.pop_back();
  if (line.find("FDR") != std::string::npos)
    method = "FDR";
  else if (line.find("FWER") != std::string::npos)
    method = "FWER";
  else if (line.find("FPR") != std::string::npos)
    method = "FPR";

  // next line contains alpha
  float alpha;
  std::getline(log_ifs, line, '\n');
  if (line.back() == '\r')
    line.pop_back();
  std::stringstream line_stream(line);
  line_stream >> discard >> alpha;

  // next line contains whether we have MaxEnt pruning
  std::string prune_MaxEnt_str;
  std::getline(log_ifs, line, '\n');
  if (line.back() == '\r')
    line.pop_back();
  line_stream = std::stringstream(line);
  line_stream >> discard >> discard >> prune_MaxEnt_str;
  bool prune_MaxEnt;
  prune_MaxEnt = (prune_MaxEnt_str == "true") ? true : false;

  // skip 8 lines
  for (uint8_t l = 0u; l < 8u; ++l)
    std::getline(log_ifs, discard, '\n');
  std::getline(log_ifs, line, '\n');
  if (line.back() == '\r')
    line.pop_back();
  line_stream = std::stringstream(line);
  uint32_t num_edges_after_threshold_pruning = 0u;
  line_stream >> discard >> discard >> discard >>
      num_edges_after_threshold_pruning;

  // skip 3 lines, the 4th contains edges after MaxEnt pruning
  float FPR_estimate_subnet;
  const uint32_t tot_possible_edges =
      regulators.size() * (defined_genes.size() - 1u);
  if (prune_MaxEnt) {
    uint32_t num_edges_after_MaxEnt_pruning = 0u;
    for (uint8_t l = 0u; l < 3u; ++l)
      std::getline(log_ifs, discard, '\n');
    std::getline(log_ifs, line, '\n');
    if (line.back() == '\r')
      line.pop_back();
    line_stream = std::stringstream(line);
    line_stream >> discard >> discard >> discard >>
        num_edges_after_MaxEnt_pruning;
    FPR_estimate_subnet = estimateFPRWithMaxEnt(
        alpha, method, num_edges_after_threshold_pruning,
        num_edges_after_MaxEnt_pruning, tot_possible_edges);
  } else {
    FPR_estimate_subnet = estimateFPRNoMaxEnt(
        alpha, method, num_edges_after_threshold_pruning, tot_possible_edges);
  }

  return std::make_pair(subnet, FPR_estimate_subnet);
}

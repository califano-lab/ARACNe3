#include "io.hpp"
#include "algorithms.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

std::string makeUnixDirectoryNameUniversal(std::string dir_name) {
  std::replace(dir_name.begin(), dir_name.end(), '/', directory_slash);
  return dir_name;
}

bool makeDirs(const std::string &dir_name, Logger *const logger) {
  if (!std::filesystem::exists(dir_name)) {
    std::filesystem::create_directories(dir_name);
    if (std::filesystem::exists(dir_name)) {
      const std::string out_msg = "Directory Created: \"" +
                                  makeUnixDirectoryNameUniversal(dir_name) +
                                  "\".";

      std::cout << out_msg << std::endl;
      if (logger)
        logger->writeLineWithTime(out_msg);

      return true;
    } else {
      const std::string err_msg = "Failed to create directory: \"" +
                                  makeUnixDirectoryNameUniversal(dir_name) +
                                  "\".";

      std::cerr << err_msg << std::endl;
      if (logger)
        logger->writeLineWithTime(err_msg);

      return false;
    }
  }
  return true;
}

void exitIfFileNotOpen(std::ifstream &ifs, const std::string &file_path,
                       Logger *const logger) {
  if (!ifs.is_open()) {
    const std::string err_msg = "Error: could not open file \"" +
                                makeUnixDirectoryNameUniversal(file_path) +
                                "\".";

    std::cerr << err_msg << std::endl;
    if (logger)
      logger->writeLineWithTime(err_msg);

    std::exit(EXIT_FAILURE);
  }

  return;
}

std::tuple<vv_float, geneset, compression_map, decompression_map>
readExpMatrixAndCopulaTransform(const std::string &exp_mat_file_path,
                                std::mt19937 &rnd, Logger *const logger) {

  std::ifstream ifs{exp_mat_file_path};
  exitIfFileNotOpen(ifs, exp_mat_file_path, logger);

  // pop the first line
  std::string cur_line;
  getline(ifs, cur_line, '\n');
  if (cur_line.back() == '\r') /* Alert! We have a Windows dweeb! */
    cur_line.pop_back();

  // parse the exp_mat
  vv_float gexp_matrix;
  geneset genes;
  compression_map compressor;
  decompression_map decompressor;

  while (std::getline(ifs, cur_line, '\n')) {
    if (cur_line.back() == '\r') /* Alert! We have a Windows dweeb! */
      cur_line.pop_back();
    std::vector<float> expr_vec;

    if (gexp_matrix.size() > 0u)
      expr_vec.reserve(gexp_matrix.at(0).size());

    // gene name is first column of matrix
    std::size_t cur_pos = 0U, cur_end = cur_line.find_first_of("\t", cur_pos);
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
    if (gexp_matrix.size() > 0U && expr_vec.size() != gexp_matrix[0].size()) {
      const std::string err_msg = "Error: genes in expression matrix do not "
                                  "all have equal samples. Check that \"" +
                                  gene + "\" has same number of samples as " +
                                  decompressor[0];

      std::cerr << err_msg << std::endl;
      if (logger)
        logger->writeLineWithTime(err_msg);

      std::exit(EXIT_FAILURE);
    }

    // copula-transform expr_vec
    expr_vec = copulaTransform(expr_vec, rnd);

    gexp_matrix.push_back(expr_vec);
  }

  return {gexp_matrix, genes, compressor, decompressor};
}

/*
 Create a subsampled gene_to_floats.  Requires that exp_mat and
n_subsample are set.
 */
const vv_float sampleExpMatAndReCopulaTransform(const vv_float &exp_mat,
                                                const uint16_t n_subsamp,
                                                std::mt19937 &rnd) {
  std::vector<uint16_t> idxs(exp_mat[0].size());
  std::iota(idxs.begin(), idxs.end(), 0U);

  std::vector<uint16_t> fold(n_subsamp);
  std::sample(idxs.begin(), idxs.end(), fold.begin(), n_subsamp, rnd);

  vv_float subsample_exp_mat(exp_mat.size(),
                             std::vector<float>(n_subsamp, 0.f));
  for (gene_id gene = 0u; gene < exp_mat.size(); ++gene) {
    for (uint16_t i = 0u; i < n_subsamp; ++i)
      subsample_exp_mat[gene][i] = exp_mat[gene][fold[i]];

    subsample_exp_mat[gene] = copulaTransform(subsample_exp_mat[gene], rnd);
  }

  return subsample_exp_mat;
}

geneset readRegList(const std::string &regulators_list_file_path,
                    const compression_map &defined_genes, Logger *const logger,
                    bool verbose) {

  std::ifstream ifs{regulators_list_file_path};
  exitIfFileNotOpen(ifs, regulators_list_file_path, logger);

  std::string cur_reg;
  geneset regulators;
  uint32_t num_warnings = 0U;
  while (std::getline(ifs, cur_reg, '\n')) {
    if (cur_reg.back() == '\r') /* Alert! We have a Windows dweeb! */
      cur_reg.pop_back();

    // Warn if not found in exp_mat once.
    if (defined_genes.find(cur_reg) == defined_genes.end()) {
      const std::string warning_msg =
          "Warning: \"" + cur_reg +
          "\" found in regulators, but not in expression matrix. Removing "
          "from analysis.";

      if (logger)
        logger->writeLineWithTime(warning_msg);

      if (num_warnings < 3u || verbose)
        std::cerr << warning_msg << std::endl;

      ++num_warnings;
    } else {
      regulators.insert(defined_genes.at(cur_reg));
    }
  }

  if (num_warnings > 3u && !verbose)
    std::cerr << "..." << num_warnings - 3U << " warnings suppressed ..."
              << std::endl;

  if (regulators.empty()) {
    const std::string abort_msg = "Abort: no regulators to analyze.";

    std::cerr << abort_msg << std::endl;
    if (logger)
      logger->writeLineWithTime(abort_msg);

    std::exit(EXIT_FAILURE);
  }

  return regulators;
}

/*
 Function that prints the Regulator, Target, and MI to the output_dir given the
 output_suffix.  Does not print to the console.  The data structure input is a
 gene_to_edge_tars, which is defined in "ARACNe3.hpp".
 */
void writeNetworkRegTarMI(const std::string &output_file_name, const char sep,
                          const gene_to_gene_to_float &network,
                          const decompression_map &decompressor) {
  std::ofstream ofs{output_file_name};
  if (!ofs) {
    std::cerr << "error: could not write to file: " << output_file_name << "."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  ofs << "regulator.values" << sep << "target.values" << sep << "mi.values"
      << std::endl;
  for (const auto &[reg, regulon] : network)
    for (const auto [tar, mi] : regulon)
      ofs << decompressor[reg] << sep << decompressor[tar] << sep << mi << '\n';

  ofs.close();
}

void writeARACNe3DF(const std::string &output_file_name, const char sep,
                    const std::vector<ARACNe3_df> &output_df,
                    const decompression_map &decompressor) {

  std::ofstream ofs(output_file_name);

  ofs << "regulator.values" << sep << "target.values" << sep << "mi.values"
      << sep << "scc.values" << sep << "count.values" << sep << "log.p.values"
      << '\n';

  for (size_t i = 0U; i < output_df.size(); ++i) {
    const ARACNe3_df &row = output_df[i];
    ofs << decompressor.at(row.regulator) << sep << decompressor.at(row.target)
        << sep << row.final_mi << sep << row.final_scc << sep
        << row.num_subnets_incident << sep << row.final_log_p << '\n';
  }

  ofs.close();
}

/*
 * Lists the files in the provided subnets dir and matches them to their log
 * files, assuming heterogeneous runids.  Returns a pair of string vectors that
 * are index-matched according to subnet, so the correct log statistics are used
 * in later computation
 */
pair_string_vecs
findSubnetFilesAndSubnetLogFiles(const std::string &subnets_dir,
                                 const std::string &subnets_log_dir) {
  std::vector<std::string> subnet_filenames, subnet_log_filenames;
  try {
    for (const auto &entry : std::filesystem::directory_iterator(subnets_dir)) {
      if (entry.is_regular_file()) {
        std::string subnet_filename = entry.path().filename().string();
        subnet_filenames.push_back(subnet_filename);

        // Construct the expected log file name, replacing .tsv with .txt
        std::string subnet_log_filename = "log-" + subnet_filename;
        size_t pos = subnet_log_filename.rfind(".tsv");
        if (pos != std::string::npos)
          subnet_log_filename.replace(pos, 4, ".txt");

        if (std::filesystem::exists(subnets_log_dir + subnet_log_filename))
          subnet_log_filenames.push_back(subnet_log_filename);
        else {
          std::cerr << "Fatal: expected \"" + subnets_log_dir +
                           subnet_log_filename +
                           "\" to exist based on the file \"" + subnets_dir +
                           subnet_filename +
                           "\", but the log file was not found."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
    }
  } catch (std::filesystem::filesystem_error &e) {
    std::cerr << "Error reading directory: " << e.what() << std::endl;
    std::cerr << "Check that your output directory has the subdirectories "
                 "\"subnets/\" and \"subnets_log/\""
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return {subnet_filenames, subnet_log_filenames};
}

/*
 Reads a subnet file and then updates the FPR_estimates vector defined in
 "subnet_operations.cpp"
 */
std::pair<gene_to_gene_to_float, float>
loadARACNe3SubnetsAndUpdateFPRFromLog(const std::string &subnet_file_path,
                                      const std::string &subnet_log_file_path,
                                      const compression_map &defined_genes,
                                      const geneset &regulators,
                                      Logger *const logger) {
  /*
   Read in the subnet file
   */
  std::ifstream subnet_ifs{subnet_file_path};
  exitIfFileNotOpen(subnet_ifs, subnet_file_path, logger);

  // discard the first line (header)
  std::string line;
  std::getline(subnet_ifs, line, '\n');
  if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
    line.pop_back();

  gene_to_gene_to_float subnet;
  while (std::getline(subnet_ifs, line, '\n')) {
    if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
      line.pop_back();

    std::size_t prev = 0U, pos = line.find_first_of("\t", prev);
    const std::string reg = line.substr(prev, pos - prev);
    prev = pos + 1;

    pos = line.find_first_of("\t", prev);
    const std::string tar = line.substr(prev, pos - prev);
    prev = pos + 1;

    for (const std::string &gene : {reg, tar})
      if (defined_genes.find(gene) == defined_genes.end()) {
        const std::string err_msg = "Error: subnetwork file(s) contain gene "
                                    "names undefined in expression matrix.";

        std::cerr << err_msg << std::endl;
        if (logger)
          logger->writeLineWithTime(err_msg);

        std::exit(EXIT_FAILURE);
      }

    const float mi = std::stof(line.substr(prev, std::string::npos));

    subnet[defined_genes.at(reg)][defined_genes.at(tar)] = mi;
  }

  /*
   Read in the log file
   */
  std::ifstream log_ifs{subnet_log_file_path};
  exitIfFileNotOpen(log_ifs, subnet_log_file_path, logger);

  // discard 9 lines
  std::string discard;
  for (uint8_t l = 0; l < 9; ++l)
    std::getline(log_ifs, discard, '\n');

  // next line contains the method
  std::string method;
  std::getline(log_ifs, line, '\n');
  if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
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
  if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
    line.pop_back();
  std::stringstream line_stream(line);
  line_stream >> discard >> alpha;

  // next line contains whether we have MaxEnt pruning
  std::string prune_MaxEnt_str;
  std::getline(log_ifs, line, '\n');
  if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
    line.pop_back();
  line_stream = std::stringstream(line);
  line_stream >> discard >> discard >> prune_MaxEnt_str;
  bool prune_MaxEnt;
  prune_MaxEnt = (prune_MaxEnt_str == "true") ? true : false;

  // skip 8 lines
  for (uint8_t l = 0U; l < 8U; ++l)
    std::getline(log_ifs, discard, '\n');
  std::getline(log_ifs, line, '\n');
  if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
    line.pop_back();
  line_stream = std::stringstream(line);
  uint32_t num_edges_after_threshold_pruning = 0U;
  line_stream >> discard >> discard >> discard >>
      num_edges_after_threshold_pruning;

  // skip 3 lines, the 4th contains edges after MaxEnt pruning
  float FPR_estimate_subnet;
  const uint32_t tot_possible_edges =
      regulators.size() * (defined_genes.size() - 1);
  if (prune_MaxEnt) {
    uint32_t num_edges_after_MaxEnt_pruning = 0U;
    for (uint8_t l = 0U; l < 3U; ++l)
      std::getline(log_ifs, discard, '\n');
    std::getline(log_ifs, line, '\n');
    if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
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

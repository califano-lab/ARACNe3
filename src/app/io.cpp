#include "io.hpp"
#include "ARACNe3.hpp"
#include "algorithms.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

std::vector<std::string> decompression_map;
static std::unordered_map<std::string, uint16_t> compression_map;

std::string makeUnixDirectoryNameUniversal(std::string &dir_name) {
  std::replace(dir_name.begin(), dir_name.end(), '/', directory_slash);
  return dir_name;
}

std::string makeUnixDirectoryNameUniversal(std::string &&dir_name) {
  std::replace(dir_name.begin(), dir_name.end(), '/', directory_slash);
  return dir_name;
}

/*
 Will automatically checked if there already is a directory.  Then creates the
 output directory.
 */
void makeDir(const std::string &dir_name) {
  if (!std::filesystem::exists(dir_name)) {
    if (std::filesystem::create_directory(dir_name)) {
      std::cout << "Directory Created: \"" + dir_name + "\"." << std::endl;
    } else {
      std::cerr
          << "Failed to create directory: \"" + dir_name +
                 "\". Make sure you have permissions over the output directory."
          << std::endl;
      std::exit(2);
    }
  }
  return;
}

/*
 Create a subsampled gene_to_floats.  Requires that exp_mat and
 tot_num_subsample are set.
 */
gene_to_floats
sampleExpMatAndReCopulaTransform(const gene_to_floats &exp_mat,
                                 const uint16_t &tot_num_subsample,
                                 std::mt19937 &rand) {
  std::vector<uint16_t> idxs(exp_mat.cbegin()->size());
  std::iota(idxs.begin(), idxs.end(), 0U);

  std::vector<uint16_t> fold(tot_num_subsample);
  std::sample(idxs.begin(), idxs.end(), fold.begin(), tot_num_subsample, rand);

  gene_to_floats subsample_exp_mat(exp_mat.size(),
                                   std::vector<float>(tot_num_subsample, 0.f));
  for (gene_id gene = 0U; gene < exp_mat.size(); ++gene) {
    for (uint16_t i = 0U; i < tot_num_subsample; ++i)
      subsample_exp_mat[gene][i] = exp_mat[gene][fold[i]];

    std::vector<uint16_t> idx_ranks =
        rankIndices(subsample_exp_mat[gene], rand);
    for (uint16_t r = 0U; r < tot_num_subsample; ++r)
      subsample_exp_mat[gene][idx_ranks[r]] =
          (r + 1) / ((float)tot_num_subsample + 1);
  }
  return subsample_exp_mat;
}

/* Reads a normalized (CPM, TPM) tab-separated (G+1)x(N+1) gene expression
 * matrix and outputs a pair containing the gene_to_floats for the entire
 * expression matrix (non-subsampled) as well as a subsampled version for every
 * subnetwork.
 */
std::tuple<const gene_to_floats, const gene_to_shorts, const geneset,
           const uint16_t>
readExpMatrixAndCopulaTransform(const std::string &filename,
                                const float &subsampling_percent,
                                std::mt19937 &rand) {
  std::ifstream ifs{filename};
  if (!ifs.is_open()) {
    std::cerr << "error: file open failed " << filename << "." << std::endl;
    std::exit(1);
  }

  uint16_t tot_num_samps = 0U;
  geneset genes;

  // for the first line, we simply want to count the number of samples
  std::string line;
  std::getline(ifs, line, '\n');
  if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
    line.pop_back();

  // count samples from number of columns in first line
  for (size_t pos = 0U;
       (pos = line.find_first_of("\t", pos)) != std::string::npos; ++pos)
    ++tot_num_samps;

  uint32_t linesread = 1U;
  gene_to_floats exp_mat;
  gene_to_shorts ranks_mat;
  while (std::getline(ifs, line, '\n')) {
    ++linesread;
    if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
      line.pop_back();
    std::vector<float> expr_vec;
    expr_vec.reserve(tot_num_samps);
    std::vector<uint16_t> expr_ranks_vec(tot_num_samps, 0U);

    std::size_t prev = 0U, pos = line.find_first_of("\t", prev);
    std::string gene = line.substr(prev, pos - prev);
    prev = pos + 1;
    while ((pos = line.find_first_of("\t", prev)) != std::string::npos) {
      if (pos > prev) {
        expr_vec.emplace_back(stof(line.substr(prev, pos - prev)));
      }
      prev = pos + 1;
    }
    expr_vec.emplace_back(stof(line.substr(prev, std::string::npos)));

    if (expr_vec.size() != tot_num_samps) {
      std::cerr
          << "Fatal: line " + std::to_string(linesread) +
                 " length is not equal to line 1 length. Rows should have the "
                 "same number of delimiters. Check that header row contains "
                 "N+1 columns (empty corner + N sample names))."
          << std::endl;
      std::exit(1);
    }

    // copula-transform expr_vec values
    {
      std::vector<uint16_t> idx_ranks = rankIndices(expr_vec, rand);
      for (uint16_t r = 0; r < tot_num_samps; ++r) {
        expr_vec[idx_ranks[r]] = (r + 1) / ((float)tot_num_samps + 1);
        expr_ranks_vec[idx_ranks[r]] = r + 1;
      }
    }

    // create compression scheme from exp_mat
    if (compression_map.find(gene) == compression_map.end()) {
      compression_map[gene] = decompression_map.size(); // str -> idx
      decompression_map.push_back(gene);                // idx -> str
      genes.insert(compression_map[gene]);

      // assumes index is same as the compression_map[gene]
      exp_mat.emplace_back(expr_vec);
      ranks_mat.emplace_back(expr_ranks_vec);
    } else {
      std::cerr << "Fatal: 2 rows corresponding to " + gene + " detected."
                << std::endl;
      std::exit(1);
    }
  }

  return std::make_tuple(exp_mat, ranks_mat, genes, tot_num_samps);
}

/*
 Reads a newline-separated regulator list and sets the decompression mapping, as
 well as the compression mapping, as file static variables hidden to the rest of
 the app.
 */
const geneset readRegList(const std::string &filename) {
  std::ifstream ifs{filename};
  if (!ifs.is_open()) {
    std::cerr << "error: file open failed \"" << filename << "\"." << std::endl;
    std::exit(1);
  }
  geneset regulators;

  std::string reg;
  while (std::getline(ifs, reg, '\n')) {
    if (reg.back() == '\r') /* Alert! We have a Windows dweeb! */
      reg.pop_back();
    if (compression_map.find(reg) == compression_map.end())
      std::cerr << "Warning: " + reg +
                       " found in regulators list, but no entry in expression "
                       "matrix. Ignoring in network generation."
                << std::endl;
    else
      regulators.insert(compression_map[reg]);
  }

  return regulators;
}

/*
 Function that prints the Regulator, Target, and MI to the output_dir given the
 output_suffix.  Does not print to the console.  The data structure input is a
 gene_to_edge_tars, which is defined in "ARACNe3.hpp".
 */
void writeNetworkRegTarMI(gene_to_gene_to_float &network,
                          const std::string &file_path) {
  std::ofstream ofs{file_path};
  if (!ofs) {
    std::cerr << "error: could not write to file: " << file_path << "."
              << std::endl;
    std::cerr << "Try making the output directory subdirectory of the working "
                 "directory. Example \"-o " +
                     makeUnixDirectoryNameUniversal("./run1") + "\"."
              << std::endl;
    std::exit(2);
  }

  ofs << "regulator.values\ttarget.values\tmi.values" << std::endl;
  for (const auto &[reg, tar_mi] : network)
    for (const auto [tar, mi] : tar_mi)
      ofs << decompression_map[reg] << '\t' << decompression_map[tar] << '\t'
          << mi << '\n';
}

void writeConsolidatedNetwork(const std::vector<consolidated_df_row> &final_df,
                              const std::string &file_path) {
  std::ofstream ofs{file_path};
  if (!ofs) {
    std::cerr << "error: could not write to file: " << file_path << "."
              << std::endl;
    std::cerr << "Try making the output directory subdirectory of the working "
                 "directory. Example \"-o " +
                     makeUnixDirectoryNameUniversal("./runs") + "\"."
              << std::endl;
    std::exit(2);
  }
  ofs << "regulator.values\ttarget.values\tmi.values\tscc.values\tcount."
         "values\tlog.p.values\n";
  for (const auto &edge : final_df)
    ofs << decompression_map[edge.regulator] << '\t'
        << decompression_map[edge.target] << '\t' << edge.final_mi << '\t'
        << edge.final_scc << '\t' << edge.num_subnets_incident << '\t'
        << edge.final_log_p
        << '\n'; // using '\n' over std::endl, better for performance
}

/*
 This function will add genes to the compression scheme in any order.  It's use
 is currently only when reading subnets.
 */
void addToCompressionVecs(const std::string &gene) {
  // If we have a new gene, put it in compression scheme
  if (compression_map.find(gene) == compression_map.end()) {
    compression_map[gene] = decompression_map.size(); // str -> idx (hashmap)
    decompression_map.push_back(gene);                // idx -> str (vector)
  }
}

/*
 * Lists the files in the provided subnets dir and matches them to their log
 * files, assuming heterogeneous runids.  Returns a pair of string vectors that
 * are index-matched according to subnet, so the correct log statistics are used
 * in later computation
 */
pair_string_vecs
findSubnetFilesAndSubnetLogFiles(const std::string &subnets_dir, const std::string &subnets_log_dir) {
  std::vector<std::string> subnet_filenames, subnet_log_filenames;
  try {
    for (const auto &entry : std::filesystem::directory_iterator(subnets_dir)) {
      if (entry.is_regular_file()) {
        std::string subnet_filename = entry.path().filename().string();
        subnet_filenames.push_back(subnet_filename);

        // Construct the expected log file name
        std::string subnet_log_filename = "log_" + subnet_filename;

        // Replace the extension .tsv with .txt
        size_t pos = subnet_log_filename.rfind(".tsv");
        if (pos != std::string::npos)
          subnet_log_filename.replace(pos, 4, ".txt");

        if (std::filesystem::exists(subnets_log_dir + subnet_log_filename))
          subnet_log_filenames.push_back(subnet_log_filename);
        else {
          std::cerr << "Fatal: expected \"" + subnets_log_dir +
                           subnet_log_filename +
                           "\" to exist based on the file \"" +
                           subnets_dir + subnet_filename +
                           "\", but the log file was not found."
                    << std::endl;
          std::exit(2);
        }
      }
    }
  } catch (std::filesystem::filesystem_error &e) {
    std::cerr << "Error reading directory: " << e.what() << std::endl;
    std::cerr << "Check that your output directory has the subdirectories "
                 "\"subnets/\" and \"subnets_log/\""
              << std::endl;
    std::exit(2);
  }
  return std::make_pair(subnet_filenames, subnet_log_filenames);
}

/*
 Reads a subnet file and then updates the FPR_estimates vector defined in
 "subnet_operations.cpp"
 */
std::pair<gene_to_gene_to_float, float>
loadARACNe3SubnetsAndUpdateFPRFromLog(const std::string &subnet_file_path,
                                      const std::string &subnet_log_file_path) {
  geneset regulators, genes;

  /*
   Read in the subnet file
   */
  std::ifstream subnet_ifs{subnet_file_path};
  if (!subnet_ifs) {
    std::cerr << "error: could not read from subnet file: " << subnet_file_path
              << "." << std::endl;
    std::cerr << "Subnet files must follow the output structure of ARACNe3. "
                 "Example \"-o " +
                     makeUnixDirectoryNameUniversal("./output") +
                     "\" will contain a subdirectory \"" +
                     makeUnixDirectoryNameUniversal("subnets/") +
                     "\", which has subnet files formatted *and named* exactly "
                     "how ARACNe3 outputs subnet files."
              << std::endl;
    std::exit(2);
  }

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

    const float mi = std::stof(line.substr(prev, std::string::npos));

    addToCompressionVecs(reg);
    addToCompressionVecs(tar);

    regulators.insert(compression_map[reg]);
    genes.insert(compression_map[tar]);

    subnet[compression_map[reg]][compression_map[tar]] = mi;
  }

  genes.insert(regulators.begin(), regulators.end());

  /*
   Read in the log file
   */
  std::ifstream log_ifs{subnet_log_file_path};
  if (!log_ifs) {
    std::cerr << "error: could read from subnet log file: "
              << subnet_log_file_path << "." << std::endl;
    std::cerr << "Subnet log files must follow the output structure of "
                 "ARACNe3. Example \"-o " +
                     makeUnixDirectoryNameUniversal("./output") +
                     "\" will contain a subdirectory \"" +
                     makeUnixDirectoryNameUniversal("subnets_log/") +
                     "\", which has subnet log files formatted *and named* "
                     "exactly how ARACNe3 outputs subnet log files."
              << std::endl;
      std::exit(2);
  }

  // discard 8 lines
  std::string discard;
  for (uint8_t l = 0; l < 8; ++l)
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
  uint32_t num_edges_after_MaxEnt_pruning = 0U;
  if (prune_MaxEnt) {
    for (uint8_t l = 0U; l < 3U; ++l)
      std::getline(log_ifs, discard, '\n');
    std::getline(log_ifs, line, '\n');
    if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
      line.pop_back();
    line_stream = std::stringstream(line);
    line_stream >> discard >> discard >> discard >>
        num_edges_after_MaxEnt_pruning;
  }

  float FPR_estimate_subnet;
  if (prune_MaxEnt) {
    if (method == "FDR")
      FPR_estimate_subnet = (alpha * num_edges_after_MaxEnt_pruning) /
                            (regulators.size() * (genes.size() - 1) -
                             (1 - alpha) * num_edges_after_threshold_pruning);
    else if (method == "FWER")
      FPR_estimate_subnet = (alpha / (regulators.size() * (genes.size() - 1))) *
                            (num_edges_after_MaxEnt_pruning) /
                            (num_edges_after_threshold_pruning);
    else if (method == "FPR")
      FPR_estimate_subnet = alpha * num_edges_after_MaxEnt_pruning /
                            num_edges_after_threshold_pruning;
  } else {
    if (method == "FDR")
      FPR_estimate_subnet = (alpha * num_edges_after_threshold_pruning) /
                            (regulators.size() * (genes.size() - 1) -
                             (1 - alpha) * num_edges_after_threshold_pruning);
    else if (method == "FWER")
      FPR_estimate_subnet = alpha / (regulators.size() * (genes.size() - 1));
    else if (method == "FPR")
      FPR_estimate_subnet = alpha;
  }

  return std::make_pair(subnet, FPR_estimate_subnet);
}

#include "io.hpp"
#include <filesystem>
#include <iostream>
#include <fstream>
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
 Will automatically checked if there already is a directory.  Then creates the output directory.
 */
void makeDir(const std::string &dir_name) {
	if (!std::filesystem::exists(dir_name)) {
		if (std::filesystem::create_directory(dir_name)) {
			std::cout << "Directory Created: \"" + dir_name + "\"." << std::endl;
		} else {
			std::cerr << "Failed to create directory: \"" + dir_name + "\". Make sure you have permissions over the output directory." << std::endl;
			std::exit(2);
		}
	}
	return;
}

/** @brief Ranks indices based on the values in vec. 
 * 
 * This function sorts the indices in the range [1, size) based on the values
 * of vec[index-1]. The returned vector represents the rank of each element in
 * vec with the smallest element ranked as 1 and the largest as size. If two
 * elements in vec have the same value, their corresponding indices in the
 * ranking are randomly shuffled.
 *
 * @param vec The input vector for which the ranking should be formed.
 * @param rand A Mersenne Twister pseudo-random generator of 32-bit numbers
 * with a state size of 19937 bits. Used to shuffle indices corresponding to
 * equal values in vec.
 * 
 * @return A vector representing the rank of indices in the input vector vec.
 *
 * @example vec = {9.2, 3.5, 7.4, 3.5} The function returns {4, 1, 3, 2}. Note
 * that the ranks for the elements with the same value 3.5 (indices 1 and 3)
 * may be shuffled differently in different runs.
 */
std::vector<uint16_t> rankIndices(const std::vector<float>& vec, std::mt19937 &rand) {
	std::vector<uint16_t> idx_ranks(vec.size());
	std::iota(idx_ranks.begin(), idx_ranks.end(), 0U); /* 0, 1, ..., size-1 */
	std::sort(idx_ranks.begin(), idx_ranks.end(), [&vec](const uint16_t &num1, const uint16_t &num2) -> bool { return vec[num1] < vec[num2]; }); /* sort ascending */
	for (uint16_t r = 0U; r < idx_ranks.size();) {
		uint16_t same_range = 1U;
		while (r + same_range < idx_ranks.size() && vec[idx_ranks[r]] == vec[idx_ranks[r+same_range]])
			++same_range; // same_range is off-end index
		if (same_range > 1U) {
			std::shuffle(idx_ranks.begin()+r, idx_ranks.begin()+r+same_range, rand);
			r = r + same_range;
		} else {
			++r;
		}
	}
	return idx_ranks;
}

/*
 Create a subsampled gene_to_floats.  Requires that exp_mat and tot_num_subsample are set.
 */
gene_to_floats sampleExpMatAndReCopulaTransform(const gene_to_floats &exp_mat, const uint16_t &tot_num_subsample, std::mt19937 &rand) {
	std::vector<uint16_t> idxs(exp_mat.cbegin()->second.size());
	std::iota(idxs.begin(), idxs.end(), 0U);

	std::vector<uint16_t> fold(tot_num_subsample);
	std::sample(idxs.begin(), idxs.end(), fold.begin(), tot_num_subsample, rand);
	
	gene_to_floats subsample_exp_mat;
	subsample_exp_mat.reserve(exp_mat.size());	
	for (const auto &[gene_id, expr_vec] : exp_mat) {
		subsample_exp_mat[gene_id] = std::vector<float>(tot_num_subsample, 0.0f);
		
		for (uint16_t i = 0U; i < tot_num_subsample; ++i)
			subsample_exp_mat[gene_id][i] = expr_vec[fold[i]];
		
		std::vector<uint16_t> idx_ranks = rank_indexes(subsample_exp_mat[gene_id], rand);
		for (uint16_t r = 0; r < tot_num_subsample; ++r)
			subsample_exp_mat[gene_id][idx_ranks[r]] = (r + 1)/((float)tot_num_subsample + 1); 
	}
	return subsample_exp_mat;
}

/* Reads a normalized (CPM, TPM) tab-separated (G+1)x(N+1) gene expression matrix and outputs a pair containing the gene_to_floats for the entire expression matrix (non-subsampled) as well as a subsampled version for every subnetwork. 
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
  for (size_t pos = 0;
       (pos = line.find_first_of("\t, ", pos)) != std::string::npos; ++pos)
    ++tot_num_samps;

  uint32_t linesread = 1U;
  gene_to_floats exp_mat;
  gene_to_shorts ranks_mat;
  while (std::getline(ifs, line, '\n')) {
    ++linesread;
    if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
      line.pop_back();
    std::vector<float> expr_vec;
    std::vector<uint16_t> expr_ranks_vec(tot_num_samps, 0U);

    expr_vec.reserve(tot_num_samps);

    std::size_t prev = 0U, pos = line.find_first_of("\t, ", prev);
    std::string gene = line.substr(prev, pos - prev);
    prev = pos + 1;
    while ((pos = line.find_first_of("\t, ", prev)) != std::string::npos) {
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
                 "N+1 columns (N sample names and the empty corner))."
          << std::endl;
      std::exit(1);
    }

    // copula-transform expr_vec values
    {
      std::vector<uint16_t> idx_ranks = rank_indexes(expr_vec, rand);
      for (uint16_t r = 0; r < tot_num_samps; ++r) {
        expr_vec[idx_ranks[r]] = (r + 1) / ((float)tot_num_samps + 1);
        expr_ranks_vec[idx_ranks[r]] = r + 1;
      }
    }

    // create compression scheme from exp_mat
    if (compression_map.find(gene) == compression_map.end()) {
      decompression_map.push_back(gene);
      compression_map[gene] = decompression_map.size() - 1; // str -> #
      genes.insert(compression_map[gene]);

      // the last index of decompression_vec is the new uint16_t
      exp_mat[compression_map[gene]] = expr_vec;
      // ranks of exp are stored for SCC later
      ranks_mat[compression_map[gene]] = expr_ranks_vec;
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
                       "matrix. Ignoring."
                << std::endl;
    else
      regulators.insert(compression_map[reg]);
  }

  return regulators;
}

/*
 Function that prints the Regulator, Target, and MI to the output_dir given the output_suffix.  Does not print to the console.  The data structure input is a gene_to_edge_tars, which is defined in "ARACNe3.hpp".
 */
void writeNetworkRegTarMI(gene_to_gene_to_float &network,
                          const std::string &output_dir,
                          const std::string &output_suffix) {
  const std::string filename = output_dir + output_suffix + ".tsv";
  std::ofstream ofs{filename};
  if (!ofs) {
    std::cerr << "error: could not write to file: " << filename << "."
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

void writeConsolidatedNetwork(const std::vector<consolidated_df_row>& final_df, std::string filename) {
	makeUnixDirectoryNameUniversal(filename);
	std::ofstream ofs{filename};
	if (!ofs) {
		std::cerr << "error: could not write to file: " << filename << "." << std::endl;
		std::cerr << "Try making the output directory subdirectory of the working directory. Example \"-o " + makeUnixDirectoryNameUniversal("./runs") + "\"." << std::endl;
		std::exit(2);
	}
	ofs << "regulator.values\ttarget.values\tmi.values\tscc.values\tcount.values\tp.values\n";
	for (const auto& edge : final_df)
		ofs << 
		decompression_map[edge.regulator] << '\t' <<
		decompression_map[edge.target] << '\t' <<
		edge.final_mi << '\t' <<
		edge.final_scc << '\t' <<
		edge.num_subnets_incident << '\t' <<
		edge.final_p << '\n'; // using '\n' over std::endl, better for performance
}

/*
 This function will add genes to the compression scheme in any order.  It's use is currently only when reading subnets.
 */
void addToCompressionVecs(const std::string &gene) {
  // If we have a new gene, put it in compression scheme
	if (compression_map.find(gene) == compression_map.end()) {
		decompression_map.push_back(gene);
		compression_map[gene] = decompression_map.size()-1;
	}
}

/*
 Reads a subnet file and then updates the FPR_estimates vector defined in "subnet_operations.cpp"
 */
gene_to_gene_to_float readARACNe3Subnet(const std::string output_dir, const uint16_t subnet_num) {
	std::string subnet_filename = output_dir + "subnets/output_subnet" + std::to_string(subnet_num) + ".txt";
	std::string log_filename = output_dir + "log/log_subnet" + std::to_string(subnet_num) + ".txt";
	
	makeUnixDirectoryNameUniversal(subnet_filename);
	makeUnixDirectoryNameUniversal(log_filename);
	
	/*
	 Read in the subnet file
	 */
	std::ifstream subnet_ifs{subnet_filename};
	if (!subnet_ifs) {
		std::cerr << "error: could read from implied subnet file: " << subnet_filename << "." << std::endl;
		std::cerr << "Try verifying that subnet files follow the output structure of ARACNe3. Example \"-o " + makeUnixDirectoryNameUniversal("./output") + "\" will contain a subdirectory \"" + makeUnixDirectoryNameUniversal("subnets/") + "\", which has subnet files formatted exactly how ARACNe3 outputs subnet files." << std::endl;
		throw TooManySubnetsRequested();
	}
	// discard the first line (header)
	std::string line;
	getline(subnet_ifs, line, '\n');
	if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
		line.pop_back();
	gene_to_edge_tars subnet;
	while(std::getline(subnet_ifs, line, '\n')) {
		if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
			line.pop_back();
		
		std::size_t prev = 0U, pos = line.find_first_of("\t", prev);
		const std::string reg = line.substr(prev, pos-prev);
		prev = pos + 1;
		
		pos = line.find_first_of("\t", prev);
		const std::string tar = line.substr(prev, pos-prev);
		prev = pos + 1;
		
		const float mi = std::stof(line.substr(prev, std::string::npos));
		
		subnet[compression_map[reg]-1].emplace_back(compression_map[tar]-1, mi);
	}
	
	/* 
	 Read in the log file
	 */
	std::ifstream log_ifs{log_filename};
	if (!log_ifs) {
		std::cerr << "error: could read from implied subnet log file: " << log_filename << "." << std::endl;
		std::cerr << "Try verifying that subnet log files follow the output structure of ARACNe3. Example \"-o " + makeUnixDirectoryNameUniversal("./output") + "\" will contain a subdirectory \"" + makeUnixDirectoryNameUniversal("log/") + "\", which has subnet log files formatted exactly how ARACNe3 outputs subnet log files." << std::endl;
		throw TooManySubnetsRequested();
	}
	// discard 8 lines
	std::string discard;
	for (uint8_t l = 0; l < 8; ++l)
		getline(log_ifs, discard, '\n');
	// next line contains the method 
	getline(log_ifs, line, '\n');
	if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
		line.pop_back();
	if (line.find("FDR") != std::string::npos)
		method = "FDR";
	else if (line.find("FWER") != std::string::npos)
		method = "FWER";
	else if (line.find("FPR") != std::string::npos)
		method = "FPR";
	// next line contains alpha
	getline(log_ifs, line, '\n');
	if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
		line.pop_back();
	std::stringstream line_stream(line);
	line_stream >> discard >> alpha;
	// next line contains whether we have MaxEnt pruning
	getline(log_ifs, line, '\n');
	if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
		line.pop_back();
	line_stream = std::stringstream(line);
	line_stream >> discard >> discard >> prune_MaxEnt;
	// skip 10 lines, the 11th contains edges after threshold pruning
	uint32_t num_edges_after_threshold_pruning = 0U;
	for (uint8_t l = 0; l < 10; ++l)
		getline(log_ifs, discard, '\n');
	getline(log_ifs, line, '\n');
	if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
		line.pop_back();
	line_stream = std::stringstream(line);
	line_stream >> discard >> discard >> discard >> num_edges_after_threshold_pruning;
	// skip 4 lines, the 5th contains edges after MaxEnt pruning if done so
	uint32_t num_edges_after_MaxEnt_pruning = 0U;
	if (prune_MaxEnt) {
		for (uint8_t l = 0; l < 4; ++l)
			getline(log_ifs, discard, '\n');
		getline(log_ifs, line, '\n');
		if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
			line.pop_back();
		line_stream = std::stringstream(line);
		line_stream >> discard >> discard >> discard >> num_edges_after_MaxEnt_pruning;
	}
	
	if(prune_MaxEnt) {
		if (method == "FDR")
			FPR_estimates.emplace_back((alpha*num_edges_after_MaxEnt_pruning)/(regulators.size()*genes.size()-(1-alpha)*num_edges_after_threshold_pruning));
		else if (method == "FWER")
			FPR_estimates.emplace_back((alpha/(regulators.size()*(genes.size()-1)))*(num_edges_after_MaxEnt_pruning)/(num_edges_after_threshold_pruning));
		else if (method == "FPR")
			FPR_estimates.emplace_back(alpha*num_edges_after_MaxEnt_pruning/num_edges_after_threshold_pruning);
	} else {
		if (method == "FDR")
			FPR_estimates.emplace_back((alpha*num_edges_after_threshold_pruning)/(regulators.size()*genes.size()-(1-alpha)*num_edges_after_threshold_pruning));
		else if (method == "FWER")
			FPR_estimates.emplace_back(alpha/(regulators.size()*(genes.size()-1)));
		else if (method == "FPR")
			FPR_estimates.emplace_back(alpha);
	}
	
	return subnet;
} 

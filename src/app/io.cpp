#include "MINDy3.hpp"
#include "algorithms.hpp"

uint16_t tot_num_samps = 0U;
std::set<gene_id> regulators, targets, modulators, genes; 
std::vector<std::string> decompression_map;

static std::unordered_map<std::string, uint16_t> compression_map;

std::string makeUnixDirectoryNameUniversal(std::string dir_name) {
	std::replace(dir_name.begin(), dir_name.end(), '/', directory_slash);
	return dir_name;
}

/*
 Will automatically checked if there already is a directory.  Then creates the output directory.
 */
void makeDir(std::string dir_name) {
	makeUnixDirectoryNameUniversal(dir_name);
	if (!std::filesystem::exists(dir_name)) {
		if (std::filesystem::create_directory(dir_name)) {
			std::cout << "Directory Created: \"" + dir_name + "\"." << std::endl;
		} else {
			std::cerr << "failed to create directory: \"" + dir_name + "\"." << std::endl;
			std::exit(1);
		}
	}
	return;
}

/* Reads a normalized (CPM, TPM) expression matrix and defines the gene expression matrix in an unordered map of gene_id -> vector.  
 */
std::pair<hash_float_vec, hash_short_vec> readExpMatrixAndCopulaTransform(std::string &filename, std::mt19937 &rand) {
	makeUnixDirectoryNameUniversal(filename);
	std::ifstream ifs{filename};
	if (!ifs.is_open()) { std::cerr << "error: file open failed " << filename << "." << std::endl; std::exit(1); }

	// for the first line, we simply want to count the number of samples
	std::string line;
	getline(ifs, line, '\n');
	if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
		line.pop_back();
	
	/*
	 Count number of samples from the number of columns in the first line
	 */ 
	for (size_t pos = 0; (pos = line.find_first_of("\t, ", pos)) != std::string::npos; ++pos)
		++tot_num_samps;
	
	// now, we can more efficiently load
	uint32_t linesread = 1U;
	hash_float_vec gexp_matrix;
	hash_short_vec gexp_ranks_matrix;
	while(std::getline(ifs, line, '\n')) {
		++linesread;
		if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
			line.pop_back();
		std::vector<float> expr_vec;
		std::vector<uint16_t> expr_ranks_vec(tot_num_samps, 0U); // must initialize because we don't replace later
		
		expr_vec.reserve(tot_num_samps);
		
		std::size_t prev = 0U, pos = line.find_first_of("\t, ", prev);
		std::string gene = line.substr(prev, pos-prev);
		prev = pos + 1;
		while ((pos = line.find_first_of("\t, ", prev)) != std::string::npos) {
			if (pos > prev) {
				expr_vec.emplace_back(stof(line.substr(prev, pos-prev)));
			}
			prev = pos + 1;
		}
		expr_vec.emplace_back(stof(line.substr(prev, std::string::npos)));
		
		/* This means that a user has inputted a matrix with unequal row 1 vs row 2 length
		 */
		if (expr_vec.size() > tot_num_samps) {
			std::cerr << "\nWARNING: ROW " + std::to_string(linesread) + " LENGTH DOES NOT EQUAL ROW 1 LENGTH.  ROWS MUST SHARE THE SAME NUMBER OF DELIMITERS.  SEGMENTATION FAULT MAY OCCUR.  CHECK THAT HEADER ROW CONTAINS G+1 COLUMNS.\n" << std::endl;
			tot_num_samps = expr_vec.size();
		}
		
		// copula-transform expr_vec values
		{
			std::vector<uint16_t> idx_ranks = rankIndices(expr_vec, rand);
			for (uint16_t r = 0; r < tot_num_samps; ++r) {
				expr_vec[idx_ranks[r]] = (r + 1)/((float)tot_num_samps + 1); // switches out expr_vec with copula-transformed values
				expr_ranks_vec[idx_ranks[r]] = r + 1;
			}
		}
		
		// Add to compression scheme if not already
		if (compression_map.find(gene) == compression_map.end()) {
			// we must have a target
			decompression_map.push_back(gene);
			compression_map[gene] = decompression_map.size()-1; 
			genes.insert(compression_map[gene]);
			
			// the last index of decompression_vec is the new uint16_t
			gexp_matrix[compression_map[gene]] = expr_vec;
			gexp_ranks_matrix[compression_map[gene]] = expr_ranks_vec;
		} else {
			gexp_matrix[compression_map[gene]] = expr_vec;
			gexp_ranks_matrix[compression_map[gene]] = expr_ranks_vec;
		}
	}
	
	return std::make_pair(gexp_matrix, gexp_ranks_matrix);
}


/* Reads a network file with significant regulator-target interactions as decided by preprocessing.  The compression scheme is initially defined here.  All unique regulators and targets are also identified here.
 */
hash_network readNetworkFile(std::string &filename) {
	makeUnixDirectoryNameUniversal(filename);
	std::ifstream ifs{filename};
	if (!ifs.is_open()) { std::cerr << "error: file open failed " << filename << "." << std::endl; std::exit(1); }
	
	// We assume that the first line contains headers
	std::string line;
	getline(ifs, line, '\n');
	if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
		line.pop_back();
	
	/* Read all regulator-target interactions into vector and set.  The set is
	 for storing unique indentifiers, and the vectors are for preserving the edge
	 information.
	 */
	std::vector<std::string> reg_vec, tar_vec;
	while(std::getline(ifs, line, '\n')) {
		if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
			line.pop_back();
		
		std::size_t prev = 0U, pos = line.find_first_of("\t, ", prev);
		std::string reg = line.substr(prev, pos-prev);
		prev = pos + 1;
		pos = line.find_first_of("\t, ", prev);
		std::string tar = line.substr(prev, pos-prev);
		
		reg_vec.emplace_back(reg);
		tar_vec.emplace_back(tar);
	}
	
	for (const auto &reg : reg_vec) {
		if (compression_map.find(reg) == compression_map.end()) {
			std::cerr << "ABORT: " + reg + " NOT IN GEXP MATRIX BUT IN NETWORK FILE.  DID YOU INCLUDE A HEADER ROW IN THE GEXP MATRIX?" << std::endl;
			std::exit(1);
		}
		regulators.insert(compression_map[reg]);
	}
	for (const auto &tar : tar_vec) {
		if (compression_map.find(tar) == compression_map.end()) {
			std::cerr << "ABORT: " + tar + " NOT IN GEXP MATRIX BUT IN NETWORK FILE.  DID YOU INCLUDE A HEADER ROW IN THE GEXP MATRIX?" << std::endl;
			std::exit(1);
		}
		targets.insert(compression_map[tar]);
	}
	
	// Store edge information 
	hash_network network;
	for (size_t i = 0U; i < reg_vec.size(); ++i)
		network[compression_map[reg_vec[i]]].insert(compression_map[tar_vec[i]]);
	return network;
}


/* Reads a newline-separated modulator list
 */
std::set<gene_id> readModList(std::string &filename) {
	makeUnixDirectoryNameUniversal(filename);
	std::ifstream ifs{filename};
	
	if (!ifs.is_open()) { std::cerr << "error: file open failed \"" << filename << "\"." << std::endl; std::exit(1); }
	
	std::string mod;
	while (std::getline(ifs, mod, '\n')) {
		if (mod.back() == '\r') /* Alert! We have a Windows dweeb! */
			mod.pop_back();
		if (compression_map.find(mod) == compression_map.end()) {
			std::cerr << "ABORT: " + mod + " NOT IN GEXP MATRIX BUT IN MODULATORS FILE" << std::endl;
			std::exit(1);
		}
		modulators.insert(compression_map[mod]);
	}
	
	return modulators;
}

/* Outputs a modulator_df vector to the specified output stream.  Returns EXIT_SUCCESS on success, EXIT_FAILURE on failure. 
 */
int outputModulatorDf(const std::vector<modulator_df> &output_df, std::ofstream &ofs) {
	
	return EXIT_SUCCESS;
}


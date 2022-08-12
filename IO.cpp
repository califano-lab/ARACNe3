#include "ARACNe3.hpp"

/*
 Compression (gene (string) to uint16_t) and decompression (uint16_t back to gene (string)) mapping.  When doing the matrix/regulator list IO, this application automatically compresses all gene identifiers into unsigned short (2B) 0-65535, as this substantially can decrease the memory load of data structures which must copy the initial values (strings, for the gene identifiers these are anywhere from 5B-25B), or refer to them through pointers (8B).  Reading less memory also can speed up computation.
 */
static std::unordered_map<std::string, uint16_t> compression_map;
static std::vector<std::string> decompression_map;

/*
 Global variables are passed from ARACNe3.cpp, which are the user-defined parameters.
 */
uint16_t tot_num_samps = 0U;
uint16_t tot_num_subsample = 0U;
uint16_t tot_num_regulators, defined_regulators = 0U;
genemap global_gm;
genemap_r global_gm_r;

extern uint32_t global_seed;
extern uint16_t num_subnets;
extern double subsampling_percent;
extern uint16_t nthreads;
extern bool adaptive;

extern float alpha;
extern std::string method;
extern std::vector<float> FPR_estimates;

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
void makeDir(std::string &dir_name) {
	makeUnixDirectoryNameUniversal(dir_name);
	if (!std::filesystem::exists(dir_name)) {
		if (std::filesystem::create_directory(dir_name)) {
			std::cout << "Directory Created: \"" + dir_name + "\"." << std::endl;
		} else {
			std::cerr << "failed to create directory: \"" + dir_name + "\"." << std::endl;
			std::exit(2);
		}
	}
	return;
}

/*
 A ranking is formed in the following way.  Indices index = [1,size) are sorted based on the ranking of vec[index-1], so that we get some new sorted set of indexes (5, 2, 9, ... ) that is the rank of each element in vec.  Rank is 1 for smallest, size for largest.
 
 For a lambda function, brackets indicate the scope of the function.
 */
std::vector<uint16_t> rank_indexes(const std::vector<float>& vec) {
	static std::mt19937 rand{global_seed++};
	std::vector<uint16_t> idx_ranks(vec.size());
	std::iota(idx_ranks.begin(), idx_ranks.end(), 0U); /* 0, 1, ..., size-1 */
	std::sort(idx_ranks.begin(), idx_ranks.end(), [&vec](const uint16_t &num1, const uint16_t &num2) -> bool { return vec[num1] < vec[num2];}); /* sort ascending */
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
 Reads a newline-separated regulator list and sets the decompression mapping, as well as the compression mapping, as file static variables hidden to the rest of the app.
 */
void readRegList(std::string &filename) {
	makeUnixDirectoryNameUniversal(filename);
	std::fstream f {filename};
	std::vector<std::string> regs;
	
	if (!f.is_open()) {
		std::cerr << "error: file open failed \"" << filename << "\"." << std::endl;
		std::exit(2);
	}
	
	std::string reg;
	
	while (std::getline(f, reg, '\n')) {
		if (reg.back() == '\r') /* Alert! We have a Windows dweeb! */
			reg.pop_back();
		regs.push_back(reg);
	}
	
	compression_map.reserve(regs.size());
	/* NOTE** This map starts from values i+1 because we are only using it to make the compression step below faster.  The compression map is redundant for every uint16_t greater than the number of regulators (i.e., contents are emptied after readExpMatrix
	*/
	for (uint16_t i = 0; i < regs.size(); ++i)
		compression_map[regs[i]] = i+1; //NOTE** it's "regulator" -> i+1
	
	/* The original regs string vector is also the decompression map.  We index this vector with uint16_t -> "gene", and hence this is how we decompress.
	 */
	tot_num_regulators = regs.size();
	decompression_map = regs;
	return;
}

/*
 Create a subsampled genemap.  Requires that global_gm and tot_num_subsample are set by readExpMat(), which must occur on program launch anyway
 */
genemap sampleFromGlobalGenemap() {
	static std::mt19937 rand{global_seed++};
	std::vector<uint16_t> idxs(tot_num_samps);
	std::iota(idxs.begin(), idxs.end(), 0U);
	std::vector<uint16_t> fold(tot_num_subsample);
	std::sample(idxs.begin(), idxs.end(), fold.begin(), tot_num_subsample, rand);
	
	genemap subsample_gm;
	subsample_gm.reserve(global_gm.size());	
	for (const auto &[gene, expr_vec] : global_gm) {
		subsample_gm[gene] = std::vector<float>(tot_num_subsample, 0.0f);
		
		for (uint16_t i = 0U; i < tot_num_subsample; ++i)
			subsample_gm[gene][i] = expr_vec[fold[i]];
		
		std::vector<uint16_t> idx_ranks = rank_indexes(subsample_gm[gene]);
		for (uint16_t r = 0; r < tot_num_subsample; ++r)
			subsample_gm[gene][idx_ranks[r]] = (r + 1)/((float)tot_num_subsample + 1); 
	}
	
	return subsample_gm;
}

/* Reads a normalized (CPM, TPM) tab-separated (G+1)x(N+1) gene expression matrix and outputs a pair containing the genemap for the entire expression matrix (non-subsampled) as well as a subsampled version for every subnetwork. 
 */
void readExpMatrix(std::string &filename) {
	makeUnixDirectoryNameUniversal(filename);
	std::fstream f{filename};
	genemap gm;
	genemap_r gm_r; //to store ranks of gexp values
	if (!f.is_open()) {
		std::cerr << "error: file open failed " << filename << "." << std::endl;
		std::exit(2);
	}

	// for the first line, we simply want to count the number of samples
	std::string line;
	getline(f, line, '\n');
	if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
		line.pop_back();
	
	/*
	 Count number of samples from the number of columns in the first line
	 */ 
	for (size_t pos = 0; (pos = line.find_first_of("\t, ", pos)) != std::string::npos; ++pos)
		++tot_num_samps;
	
	// find subsample number
	tot_num_subsample = std::ceil(subsampling_percent * tot_num_samps);
	if (tot_num_subsample >= tot_num_samps || tot_num_subsample < 0)
		tot_num_subsample = tot_num_samps;
	
	// now, we can more efficiently load
	uint32_t linesread = 1U;
	while(std::getline(f, line, '\n')) {
		++linesread;
		if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
			line.pop_back();
		std::vector<float> expr_vec;
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
			std::cerr << std::endl << "WARNING: ROW " + std::to_string(linesread) + " LENGTH DOES NOT EQUAL ROW 1 LENGTH.  THIS MAY BE BECAUSE THE NUMBER OF HEADER COLUMNS DOES NOT EQUAL THE NUMBER OF COLUMS IN THE MATRIX.  SEGMENTATION FAULT MAY OCCUR.  CHECK THAT SIZE OF MATRIX IS G+1 x N+1." << std::endl;
			tot_num_samps = expr_vec.size();
			tot_num_subsample = std::ceil(subsampling_percent * tot_num_samps);
			if (tot_num_subsample >= tot_num_samps || tot_num_subsample < 0)
				tot_num_subsample = tot_num_samps;
		}
				
		// copula-transform expr_vec values
		std::vector<uint16_t> idx_ranks = rank_indexes(expr_vec);
		
		std::vector<uint16_t> expr_vec_ranked(tot_num_samps);
		for (uint16_t r = 0; r < tot_num_samps; ++r) {
			expr_vec_ranked[idx_ranks[r]] = r + 1; // ranks the values of expr_vec
			expr_vec[idx_ranks[r]] = (r + 1)/((float)tot_num_samps + 1); // switches out expr_vec with copula-transformed values
		}
			
		
		/*
		 This compression works as follows.  When you input a key (gene) not in the table, it is immediately value initialized to uint16_t = 0.  However, no values are 0 in the table, as we added 1 to the index (see NOTE** above).  Note that *as soon as* we try to check if there exists 'gene' as a KEY, it is instantaneously made into a "key" with its own bin.
		 */
		if (compression_map[gene] == 0) {
			// we must have a target
			decompression_map.push_back(gene);
			compression_map[gene] = decompression_map.size(); // note: it's i+1
			
			// the last index of decompression_vec is the new uint16_t
			gm[decompression_map.size()-1] = expr_vec;
			gm_r[decompression_map.size()-1] = expr_vec_ranked; //store ranks of idx's for SCC later
		} else {
			/* we already mapped this regulator, so we must use the std::string map to find its compression value.  We do -1 because of NOTE** above */
			gm[compression_map[gene]-1] = expr_vec;
			gm_r[compression_map[gene]-1] = expr_vec_ranked; //store ranks of idx's for SCC later
		}

	}
	
	// now we must determine how many regulators are actually defined in the expression profile
	for (gene_id_t reg = 0; reg < tot_num_regulators; ++reg)
		if (global_gm.contains(reg))
			++defined_regulators;
		else
			std::cout << "WARNING: REGULATOR " + decompression_map[reg] + " DOES NOT HAVE A DEFINED GENE EXPRESSION PROFILE." << std::endl;
	
	std::cout << std::endl << "Initial Num Samples: " + std::to_string(tot_num_samps) << std::endl;
	std::cout << "Sampled Num Samples: " + std::to_string(tot_num_subsample) << std::endl << std::endl;
	
	global_gm = gm;
	global_gm_r = gm_r;
	return;
}

/*
 Function that prints the Regulator, Target, and MI to the output_dir given the output_suffix.  Does not print to the console.  The data structure input is a reg_web, which is defined in "ARACNe3.hpp".
 */
void writeNetworkRegTarMI(const reg_web &network, std::string &output_dir, const std::string &output_suffix) {
	makeUnixDirectoryNameUniversal(output_dir);
	const std::string filename = output_dir + "output_" + output_suffix + ".txt";
	std::ofstream ofs{filename};
	if (!ofs) {
		std::cerr << "error: could not write to file: " << filename << "." << std::endl;
		std::cerr << "Try making the output directory subdirectory of the working directory. Example \"-o " + makeUnixDirectoryNameUniversal("./runs") + "\"." << std::endl;
		std::exit(2);
	}
	
	ofs << "regulator.values\ttarget.values\tmi.values" << std::endl;
	for (auto it = network.begin(); it != network.end(); ++it) {
		for (auto &edge_tar : it->second) {
				ofs << decompression_map[it->first] << '\t' << decompression_map[edge_tar.target] << '\t' << edge_tar.mi << '\n'; // using '\n' over std::endl, better for performance
		}
	}
}

void writeConsolidatedNetwork(const std::vector<consolidated_df>& final_df, std::string filename) {
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
	if (compression_map[gene] == 0) {
		// we must have a new gene
		decompression_map.push_back(gene);
		// the last index of decompression_vec is the new uint16_t
		compression_map[gene] = decompression_map.size();
	}
}

/*
 Reads a subnet file and then updates the FPR_estimates vector defined in "Consolidator.cpp"
 */
reg_web readSubNetAndUpdateFPRFromLog(const std::string &output_dir, const uint16_t subnet_num) {
	std::string subnet_filename = output_dir + "subnets/output_subnet" + std::to_string(subnet_num) + ".txt";
	std::string log_filename = output_dir + "log/log_subnet" + std::to_string(subnet_num) + ".txt";
	
	makeUnixDirectoryNameUniversal(subnet_filename);
	makeUnixDirectoryNameUniversal(log_filename);
	
	std::ifstream subnet_ifs{subnet_filename};
	if (!subnet_ifs) {
		std::cerr << "error: could read from implied subnet file: " << subnet_filename << "." << std::endl;
		std::cerr << "Try verifying that subnet files follow the output structure of ARACNe3. Example \"-o " + makeUnixDirectoryNameUniversal("./output") + "\" will contain a subdirectory \"" + makeUnixDirectoryNameUniversal("subnets/") + "\", which has subnet files formatted exactly how ARACNe3 outputs subnet files." << std::endl;
		std::exit(2);
	}
	// discard the first line (header)
	std::string line;
	getline(subnet_ifs, line, '\n');
	if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
		line.pop_back();
	reg_web subnet;
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
	
	uint32_t num_edges_after_threshold_pruning, num_edges_after_MaxEnt_pruning;
	uint16_t defined_regulators;
	std::ifstream log_ifs{log_filename};
	if (!log_ifs) {
		std::cerr << "error: could read from implied subnet log file: " << log_filename << "." << std::endl;
		std::cerr << "Try verifying that subnet log files follow the output structure of ARACNe3. Example \"-o " + makeUnixDirectoryNameUniversal("./output") + "\" will contain a subdirectory \"" + makeUnixDirectoryNameUniversal("log/") + "\", which has subnet log files formatted exactly how ARACNe3 outputs subnet log files." << std::endl;
		std::exit(2);
	}
	// discard 3 lines
	for (uint8_t l = 0; l < 3; ++l) {
		getline(log_ifs, line, '\n');
		if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
			line.pop_back();
	}
	// next line contains the # of defined regulators
	getline(log_ifs, line, '\n');
	if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
		line.pop_back();
	
	if (line.find("FDR", 0) != std::string::npos)
		method = "FDR";
	else
		method = "FWER";
	// next line contains alpha
	getline(log_ifs, line, '\n');
	if (line.back() == '\r') /* Alert! We have a Windows dweeb! */
		line.pop_back();
	std::stringstream linestr(line);
	std::string discard;
	linestr >> discard >> alpha;
	// skip 10 lines, the 11th contains edges after threshold pruning
	
		
	
//TODO: You must have a way to interpret log files and generate FPR estimates from them.
	
#if 0
	method
	alpha 
	num_edges_after_threshold_pruning)
	tot_num_regulators
	tot_num_targets
	
	if(prune_maxent) {
		if (method == "FDR")
			FPR_estimates.emplace_back((alpha*num_edges_after_MaxEnt_pruning)/(defined_regulators*global_gm.size()-(1-alpha)*num_edges_after_threshold_pruning));
		else if (method == "FWER")
			FPR_estimates.emplace_back((alpha/(defined_regulators*(global_gm.size()-1)))*(num_edges_after_MaxEnt_pruning)/(num_edges_after_threshold_pruning));
	} else {
		if (method == "FDR")
			FPR_estimates.emplace_back((alpha*num_edges_after_threshold_pruning)/(defined_regulators*global_gm.size()-(1-alpha)*num_edges_after_threshold_pruning));
		else if (method == "FWER")
			FPR_estimates.emplace_back(alpha/(defined_regulators*(global_gm.size()-1)));
	}
	
	FPR_estimates.push_back();
#endif
	
	return subnet;
} 

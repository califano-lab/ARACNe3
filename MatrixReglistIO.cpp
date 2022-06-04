#include "ARACNe3.hpp"

using namespace std;

extern uint16_t tot_num_samps;
extern uint16_t tot_num_regulators;
extern bool verbose;

static unordered_map<string, uint16_t> compression_map;
static vector<string> decompression_map;

/*
 Will automatically checked if there already is an output directory.  Then creates the output directory.
 */
void makeOutputDir(const std::string &output_dir) {
	if (!std::filesystem::exists(output_dir)) {
		if (std::filesystem::create_directory(output_dir)) {
			if (verbose) std::cout << "Directory Created: " + output_dir << std::endl;
		} else {
			std::cerr << "failed to create directory: " + output_dir << std::endl;
			std::exit(2);
		}
	}
	return;
}

/*
 Will automatically checked if there already is a cached directory.  Then creates the cached directory.
 */
void makeCachedDir(const std::string &cached_dir) {
	if (!std::filesystem::exists(cached_dir)) {
		if (std::filesystem::create_directory(cached_dir)) {
			if (verbose) std::cout << "Directory Created: " + cached_dir << std::endl;
		} else {
			std::cerr << "failed to create directory: " + cached_dir << std::endl;
			std::exit(2);
		}
	}
	return;
}

/*
 Reads a newline-separated regulator list and sets the decompression mapping, as well as the compression mapping, as file static variables hidden to the rest of the app.
 */
void readRegList(string filename) {
	fstream f {filename};
	vector<string> regs;
	if (!f.is_open()) {
        	cerr << "error: file open failed " << filename << ".\n";
		std::exit(2);
	}
	string reg;
	while (getline(f, reg, '\n')) {
		regs.push_back(reg);
	}
	
	compression_map.reserve(regs.size());
	/* NOTE** This map starts from values i+1 because we are only using it to make the compression step below faster.  The compression map is redundant for every uint16_t greater than the number of regulators (i.e., contents are emptied after readExpMatrix
	*/
	for (uint16_t i = 0; i < regs.size(); ++i) {
		compression_map[regs[i]] = i+1; //NOTE** it's "regulator" -> i+1
	}
	
	/* The original regs string vector is also the decompression map.  We index this vector with uint16_t -> "gene", and hence this is how we decompress.
	 */
	tot_num_regulators = regs.size();
	decompression_map = regs;
	return;
}

/*
 * Reads a copula-transformed tab-separated (G+1)x(N+1) gene expression matrix and 
 * outputs an unordered hash table corresponding to the {gene, expression} 
 * values
 */
genemap readExpMatrix(string filename, double subsampling_percent) {
	fstream f {filename};
	genemap gm;
	if (!f.is_open()) {
        	cerr << "error: file open failed " << filename << ".\n";
		std::exit(2);
	}

	// for the first line, we simply want to count the number of samples
	string line;
	getline(f, line, '\n');
	for (size_t i = 0; (i = line.find('\t', i)) != string::npos; ++i) {
		++tot_num_samps;
	}
	
	// find subsample number **NOTE** must update tot_num_samps after subsampling
	uint16_t subsample_quant = std::round(subsampling_percent * tot_num_samps);
	std::vector<uint16_t> samps_idx(tot_num_samps), fold(subsample_quant);
	std::iota(samps_idx.begin(), samps_idx.end(), 0U); /* 0, 1, ..., size-1 */
	// now, fold is a vector with subsample_quant indices sampled from [0,tot_num_samps).
	std::sample(samps_idx.begin(), samps_idx.end(), fold.begin(), subsample_quant, std::mt19937{std::random_device{}()});
	
	// now, we can more efficiently load
	string gene;
	while(getline(f, line, '\n')) {
		vector<float> expr_vec;
		expr_vec.reserve(tot_num_samps);
		
		size_t prev = 0U, pos = line.find_first_of("\t", prev);
		gene = line.substr(prev, pos-prev);
		prev = pos + 1;
		while ((pos = line.find_first_of("\t", prev)) != string::npos) {
			if (pos > prev) {
				expr_vec.emplace_back(stof(line.substr(prev, pos-prev)));
			}
			prev = pos + 1;
		}
		expr_vec.emplace_back(stof(line.substr(prev, string::npos)));

		// subsample
		vector <float> expr_vec_sampled;
		expr_vec_sampled.reserve(subsample_quant);
		for (uint16_t i = 0U; i < subsample_quant; ++i)
			expr_vec_sampled.emplace_back(expr_vec[fold[i]]);
		
		/*
		 A ranking is formed in the following way.  Indices index = [0,subsample_quant) are sorted based on the ranking of expr_vec_sampled[index], so that we get some new sorted set of indexes (5, 2, 9, ... ) that is the rank of each element in expr_vec_sampled
		 
		 For a lambda function, brackets indicate the scope of the function.
		 */
		std::vector<uint16_t> rank_vec(subsample_quant);
		std::iota(rank_vec.begin(), rank_vec.end(), 0U); /* 0, 1, ..., size-1 */
		std::sort(rank_vec.begin(), rank_vec.end(), [&expr_vec_sampled](const uint16_t &num1, const uint16_t &num2) -> bool { return expr_vec_sampled[num1] < expr_vec_sampled[num2];}); /* sort ascending */
		
		/*
		 Lambda function to copula transform rank_vec, store in expr_vec_sampled.  Note that rank_vec goes from 0...size-1, so we must add 1 and then divide by size + 1
		 */
		std::transform(rank_vec.begin(), rank_vec.end(), expr_vec_sampled.begin(), [&subsample_quant](const uint16_t &rank) -> const float { return (rank+1)/((float) subsample_quant+1); });
		
		/*
		 This compression works as follows.  When you input a key (gene) not in the table, it is immediately value initialized to uint16_t = 0.  However, no values are 0 in the table, as we added 1 to the index (see NOTE** above).  Note that *as soon as* we try to check if tehre exists 'gene' as a KEY, it is instantaneously made into a "key" with its own bin.
		 */
		if (compression_map[gene] == 0) {
			// we must have a target
			decompression_map.push_back(gene);
			// the last index of decompression_vec is the new uint16_t
			gm[decompression_map.size()-1] = expr_vec_sampled;
		} else {
			/* we already mapped this regulator, so we must use the string map to find its compression value.  We do -1 because of NOTE** above */
			gm[compression_map[gene]-1] = expr_vec_sampled;
		}

	}
	
	if (verbose) {
		std::cout << "Initial Num Samples: " + std::to_string(tot_num_samps) << std::endl;
		std::cout << "Sampled Num Samples: " + std::to_string(subsample_quant) << std::endl;
	}
	
	// update tot_num_samps because of **NOTE** above
	tot_num_samps = subsample_quant;
	return gm;
}

/*
 Function that prints the Regulator, Target, and MI to the output_dir given the output_suffix.  Does not print to the console.  The data structure input is a reg_web, which is defined in "ARACNe3.hpp".
 */
void writeNetworkRegTarMI(const reg_web &network, const std::string &output_dir, const std::string &output_suffix) {
	const string filename = output_dir + "output_" + output_suffix + ".txt";
	ofstream ofs{filename};
	if (!ofs) {
		cerr << "error: could not write to file: " << filename << ".\n";
		cerr << "Try using the working directory. Example \"-o output.text\"." << endl;
		std::exit(2);
	}
	auto cout_buff = cout.rdbuf();
	cout.rdbuf(ofs.rdbuf());
	
	cout << "REGULATOR\tTARGET\tMI\t" << endl;
	for (auto it = network.begin(); it != network.end(); ++it) {
		for (auto &edge_tar : it->second) {
				cout << decompression_map[it->first] << '\t' << decompression_map[edge_tar.target] << '\t' << edge_tar.mi << '\t' << endl;
		}
	}
	
	cout.rdbuf(cout_buff);
}

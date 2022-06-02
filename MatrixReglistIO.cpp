#include "ARACNe3.hpp"

using namespace std;

extern uint16_t tot_num_samps;
extern uint16_t tot_num_regulators;

static unordered_map<string, uint16_t> compression_map;
static vector<string> decompression_map;

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
	/* NOTE** This map starts from values i+1 because we are only using it to make the compression step below faster.  The compression map is redundant for every uint16_t greater than the number of regulators (i.e., contents are emptied after readTransformedGexpMatrix
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
genemap readTransformedGexpMatrix(string filename) {
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
	
	// now, we can more efficiently load
	string gene;
	float expr_vals[tot_num_samps];
	while(getline(f, line, '\n')) {
		size_t prev = 0U, pos, i = 0U;
		pos = line.find_first_of("\t", prev);
		gene = line.substr(prev, pos-prev);
		prev = pos + 1;
		while ((pos = line.find_first_of("\t", prev)) != string::npos) {
			if (pos > prev) {
				expr_vals[i++] = stof(line.substr(prev, pos-prev));
			}
			prev = pos + 1;
		}
		expr_vals[tot_num_samps-1] = stof(line.substr(prev, string::npos));

		// copy the array to vector
		vector<float> expr_vec(&expr_vals[0], &expr_vals[tot_num_samps]);
		
		/*
		 This compression works as follows.  When you input a key (gene) not in the table, it is immediately value initialized to uint16_t = 0.  However, no values are 0 in the table, as we added 1 to the index (see NOTE** above).  Note that *as soon as* we try to check if tehre exists 'gene' as a KEY, it is instantaneously made into a "key" with its own bin.
		 */
		if (compression_map[gene] == 0) {
			// we must have a target
			decompression_map.push_back(gene);
			// the last index of decompression_vec is the new uint16_t
			gm[decompression_map.size()-1] = expr_vec;
		} else {
			/* we already mapped this regulator, so we must use the string map to find its compression value.  We do -1 because of NOTE** above */
			gm[compression_map[gene]-1] = expr_vec;
		}

	}
	return gm;
}

/*
 Function that prints the Regulator, Target, and MI to the given filename.  Does not print to the console.  The data structure input is a reg_web, which is defined in "ARACNe3.hpp".
 */
void printNetworkRegTarMI(const reg_web &network, const string &filename) {
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

#include "ARACNe3.hpp"

using namespace std;

static uint16_t samps = 0;

/*
 * Reads a newline-separated regulator list and outputs a vector of strings
 */
vector<string> readRegList(string filename = "regulators.txt") {
	fstream f {filename};
	vector<string> regs;
	if (!f.is_open()) {
        	cerr << "error: file open failed " << filename << ".\n";
		return regs;
	}
	string reg;
	while (getline(f, reg, '\n')) {
		regs.push_back(reg);
	}
	return regs;
}

/*
 * Reads a copula-transformed tab-separated (G+1)x(N+1) gene expression matrix and 
 * outputs an unordered hash table corresponding to the {gene, expression} 
 * values
 *
 */
genemap readTransformedGexpMatrix(string filename = "exp_mat.txt") {
	fstream f {filename};
	genemap gm;
	if (!f.is_open()) {
        	cerr << "error: file open failed " << filename << ".\n";
		return gm;
	}

	// for the first line, we simply want to count the number of samples
	string line;
	getline(f, line, '\n');
	for (size_t i = 0; (i = line.find('\t', i)) != string::npos; ++i) {
		++samps;
	}
	
	// now, we can more efficiently load
	string gene;
	float expr_vals[samps];
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
		expr_vals[samps-1] = stof(line.substr(prev, string::npos));

		// copy the array to vector
		vector<float> expr_vec(&expr_vals[0], &expr_vals[samps]);
		gm[gene] = expr_vec;
	}
	return gm;
}

/*
 * Reads a list of regulators and returns them as a string vector, along with a
 * hash map of string -> unsigned short, which is necessary to speed up the
 * compression of the gexp matrix
 */
const tuple<const vector<string>, const string_map> readRegList_compressed(string
		filename = "regulators.txt") {
	vector<string> regs_vec = readRegList(filename);
	string_map sm;
	// NOTE** This map starts from values i+1 because we are only using it
	// to make the compression step below faster.
	for (uint16_t i = 0; i < regs_vec.size(); ++i) {
		sm[regs_vec[i]] = i+1; //NOTE** it's "regulator" -> i+1
	}
	return  make_tuple(regs_vec, sm);
}

/*
 * Reads a copula-transformed tab-separated (G+1)x(N+1) gene expression matrix
 * and takes a tuple returned from readRegList_compressed, and outputs a
 * "compressed" unordered hash table corresponding to the {gene, expression}
 * values (unsigned short -> expression), as well as the compression mapping of
 * unsigned short -> gene.  This compression mapping is just a vector of
 * strings, as we can conveniently index the vector to obtain the mapping.  
 */ 
const tuple<genemap_compressed, vector<string>>
readTransformedGexpMatrix_compressed(const tuple<vector<string>, const
		string_map> regulator_mapping, string filename = "exp_mat.txt")
{ 
	fstream f {filename};
	genemap_compressed gm_c;
	vector<string> decompression_vec(get<0>(regulator_mapping));
	string_map regulator_sm = get<1>(regulator_mapping);
	// Assumes we will have less than 22,000 extra genes; it is not a
	// problem if we have more than 22,000.
	gm_c.reserve(22000);
	decompression_vec.reserve(22000);
	regulator_sm.reserve(22000);

	if (!f.is_open()) {
        	cerr << "error: file open failed " << filename << ".\n";
		return make_tuple(gm_c, decompression_vec);
	}

	// for the first line, we simply want to count the number of samples
	string line;
	getline(f, line, '\n');
	for (size_t i = 0; (i = line.find('\t', i)) != string::npos; ++i) {
		++samps;
	}
	
	// now, we can more efficiently load
	string gene;
	float expr_vals[samps];
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
		expr_vals[samps-1] = stof(line.substr(prev, string::npos));

		// copy the array to vector
		vector<float> expr_vec(&expr_vals[0], &expr_vals[samps]);

		/*
		 * The compression is quicker with the regulator_sm because as
		 * follows: Checking whether a string has already been mapped is
		 * an O(N) operation for each time we want to add a target to
		 * the decompression_vec.  However, we can make this O(1) if we
		 * use the string_map.  When you input a key not in the table,
		 * it is immediately value initialized (to uint16_t 0), and so
		 * this is an O(1) check if the gene has already been mapped.
		 * Note that *as soon as* we try to check if there exists
		 * regulator_sm[gene], 'gene' is now a KEY in the table, with
		 * its own bin.
		 */
		if (regulator_sm[gene] == 0) {
			// we must have a target
			decompression_vec.emplace_back(gene);
			// the last index of decompression_vec is the new uint16_t
			gm_c[decompression_vec.size()-1] = expr_vec;
		} else {
			// we already mapped this regulator, we must use the
			// string map to find its compression value.  We do -1
			// because of NOTE** above
			gm_c[regulator_sm[gene]-1] = expr_vec;
		}
	}
	/*
	 * The final compression result; everything needed to use compressed
	 * values and to decompress the values.  gm_c will map all unique
	 * integral values to their expression vector, and decompression_vec
	 * will map those integral values to their genes as strings by index.
	 * The first #reg indices (indices 0 to #regs) are the regulators, and
	 * the number of regulators is known by the size of the vector returned
	 * from readRegList_compressed.
	 */
	return make_tuple(gm_c, decompression_vec);
}


//int main() {
//     string regulator_file = "test/regulators.txt";
//     string matrix_file = "test/exp_mat.txt";
//
//     const vector<string> regulators = readRegList(regulator_file);
//     for (string reg : regulators) { cout << reg << endl; }
//     
//     genemap expression = readTransformedGexpMatrix(matrix_file);
//     for (string reg : regulators) { cout << reg << " " << 
//     	expression[reg][0] << " " <<
//     	expression[reg][samps - 1] << endl; }
//     cout << "SAMPS: " << samps << endl;
//
//     return 0;
//}

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

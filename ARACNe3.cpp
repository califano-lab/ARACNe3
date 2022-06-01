#include "ARACNe3.hpp"

using namespace std;
// inferred while reading txt files.  Rcpp will have to compensate
uint32_t size_of_network_unpruned = 0;
extern uint16_t tot_num_samps;

/*
 Macro for timing parts of ARACNe3.  Will set last.
 */
auto last = chrono::high_resolution_clock::now(), cur = chrono::high_resolution_clock::now();
void sinceLast() {
	cur = chrono::high_resolution_clock::now();
	cout << chrono::duration_cast<chrono::milliseconds>(cur-last).count() << "ms" << endl;
	last = cur;
}


/*
 * Assumes that simply the path to the regulator list and the path to the gene
 * expression matrix are commandline arguments
 *
 * e.g. ./ARACNe3 test/regfile.txt test/matrixfile.txt
 */
int main(int argc, char *argv[]) {
	bool prune_FDR = true;
	bool prune_DPI = false;
	
	vector<string> regs = readRegList(string(argv[1]));
	genemap matrix = readTransformedGexpMatrix(string(argv[2]));
	size_of_network_unpruned = static_cast<uint32_t>(regs.size()*matrix.size()-regs.size());
	
	//-------time module-------
	cout << "INIT NULL BEGIN" << endl;
	last = chrono::high_resolution_clock::now();
	//-------------------------
	
	initNullMIs(tot_num_samps);
	
	//-------time module-------
	cout << "INIT NULL DONE" << endl;
	sinceLast();
	//-------------------------
	
	//-------time module-------
	cout << "COMPUTING REGULATOR \"WEBS\" NO P BEGIN" << endl;
	last = chrono::high_resolution_clock::now();
	//-------------------------
	
	reg_web network;
	network.reserve(regs.size());
	for (auto &reg : regs) {
		network[reg] = genemapAPMI(matrix, reg, 7.815, 4);
	}
	
	//-------time module-------
	cout << "REGULATOR \"WEBS\" DONE" << endl;
	sinceLast();
	//-------------------------
	
	if (prune_FDR) {
		//-------time module-------
		cout << "FDR PRUNING BEGIN" << endl;
		last = chrono::high_resolution_clock::now();
		//-------------------------
		/*
		 We could prune in-network, but that would require many search operations.  It is better to extract edges and reform the entire network, then free memory, it seems.
		 */
		reg_web pruned = pruneFDR(network, regs, size_of_network_unpruned, 0.05f);
		
		// frees some memory as well
		reg_web().swap(network);
		
		//-------time module-------
		cout << "FDR PRUNING DONE" << endl;
		sinceLast();
		//-------------------------
		
		
		//-------time module-------
		cout << "PRINTING NETWORK REG-TAR-MI" << endl;
		last = chrono::high_resolution_clock::now();
		//-------------------------
		
		printNetworkRegTarMI(pruned, "output.txt");
		
		//-------time module-------
		cout << "PRINTING DONE" << endl;
		sinceLast();
		//-------------------------
		
		if (prune_DPI) {
			//-------time module-------
			cout << "DPI PRUNING BEGIN" << endl;
			last = chrono::high_resolution_clock::now();
			//-------------------------
			
			//-------time module-------
			cout << "DPI PRUNING DONE" << endl;
			last = chrono::high_resolution_clock::now();
			//-------------------------
		}
	}
}

/* timing funcs
#include <chrono>
auto start = chrono::high_resolution_clock::now();
auto end = chrono::high_resolution_clock::now();
cout << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms" << endl;
 */

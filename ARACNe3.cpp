#include "ARACNe3.hpp"

using namespace std;
// inferred while reading txt files.  Rcpp will have to compensate
uint32_t size_of_network = 0;
uint16_t tot_num_samps = 0;
uint16_t tot_num_regulators = 0;
bool prune_FDR = false;
bool prune_DPI = false;

/*
 Convenient function for timing parts of ARACNe3.  Will set last.
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
	prune_FDR = true;
	prune_DPI = true;
	
	readRegList(string(argv[1]));
	genemap matrix = readTransformedGexpMatrix(string(argv[2]));
	size_of_network = static_cast<uint32_t>(tot_num_regulators*matrix.size()-tot_num_regulators);
	
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
	cout << "COMPUTING UNPRUNED NETWORK BEGIN" << endl;
	last = chrono::high_resolution_clock::now();
	//-------------------------
	
	reg_web network;
	network.reserve(tot_num_regulators);
	for (uint16_t reg = 0; reg < tot_num_regulators; ++reg) {
		network[reg] = genemapAPMI(matrix, reg, 7.815, 4);
		
	}
	
	//-------time module-------
	cout << "UNPRUNED NETWORK DONE" << endl;
	sinceLast();
	cout << "SIZE OF NETWORK: " << size_of_network << " EDGES." << endl;
	//-------------------------
	
	if (prune_FDR) {
		//-------time module-------
		cout << "FDR PRUNING BEGIN" << endl;
		last = chrono::high_resolution_clock::now();
		//-------------------------
		/*
		 We could prune in-network, but that would require many search operations.  It is better to extract edges and reform the entire network, then free memory, it seems.
		 */
		network = pruneFDR(network, size_of_network, 0.05f);
		
		//-------time module-------
		cout << "FDR PRUNING DONE" << endl;
		sinceLast();
		cout << "SIZE OF NETWORK: " << size_of_network << " EDGES." << endl;
		//-------------------------
		
		if (prune_DPI) {
			//-------time module-------
			cout << "DPI PRUNING BEGIN" << endl;
			last = chrono::high_resolution_clock::now();
			//-------------------------
			
			network = pruneDPI(network);
			
			//-------time module-------
			cout << "DPI PRUNING DONE" << endl;
			last = chrono::high_resolution_clock::now();
			cout << "SIZE OF NETWORK: " << size_of_network << " EDGES." << endl;
			//-------------------------
		}
	}
	
	//-------time module-------
	cout << "PRINTING NETWORK REG-TAR-MI" << endl;
	last = chrono::high_resolution_clock::now();
	//-------------------------
	
	printNetworkRegTarMI(network, "output.txt");
	
	//-------time module-------
	cout << "PRINTING DONE" << endl;
	sinceLast();
	//-------------------------
}

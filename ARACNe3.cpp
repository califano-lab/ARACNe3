#include "ARACNe3.hpp"
#include <chrono>

using namespace std;
extern uint16_t tot_num_samps;

/*
 * Assumes that simply the path to the regulator list and the path to the gene
 * expression matrix are commandline arguments
 *
 * e.g. ./ARACNe3 test/regfile.txt test/matrixfile.txt
 */
int main(int argc, char *argv[]) {
	bool prune_1 = true;
	
	const vector<string> regs = readRegList(string(argv[1]));
	genemap matrix = readTransformedGexpMatrix(string(argv[2]));

	//ofstream ofs{"output.txt"};
	//auto cout_buff = cout.rdbuf();
	//cout.rdbuf(ofs.rdbuf());
	
	if (prune_1) {
		cout << "REGULATOR\tTARGET\tMI\tP-Value" << endl;
		
		
		cout << "INIT NULL" << endl;
		//-------time module-------
		auto last = chrono::high_resolution_clock::now();
		//-------------------------
		
		initNullMIs(tot_num_samps);
		
		cout << "NULL DONE" << endl;
		//-------time module-------
		auto cur = chrono::high_resolution_clock::now();
		cout << chrono::duration_cast<chrono::milliseconds>(cur - last).count() << "ms" << endl;
		last = cur;
		//-------------------------
		
		cout << "COMPUTING REGULATOR \"WEBS\" WITH P" << endl;
		reg_web_p reg_edge_tars_p;
		reg_edge_tars_p.reserve(regs.size());
		for (auto &reg : regs) {
			reg_edge_tars_p[reg] = genemapAPMI_p(matrix, reg, 7.815, 4);
		}
		cout << "REGULATOR \"WEBS\" DONE" << endl;
		//-------time module-------
		cur = chrono::high_resolution_clock::now();
		cout << chrono::duration_cast<chrono::milliseconds>(cur - last).count() << "ms" << endl;
		last = cur;
		//-------------------------
		
		for (auto it = reg_edge_tars_p.begin(); it != reg_edge_tars_p.end(); ++it) {
			for (auto &edge_tar_p : it->second) {
				cout << it->first << '\t' << edge_tar_p.target << '\t' << edge_tar_p.mi << '\t' << edge_tar_p.p_value << endl;
			}
		}
		
	} else {
		cout << "REGULATOR\tTARGET\tMI" << endl;
		reg_web reg_edge_tars;
		reg_edge_tars.reserve(regs.size());
		for (auto &reg : regs) {
			reg_edge_tars[reg] = genemapAPMI(matrix, reg, 7.815, 4);
		}
		for (auto it = reg_edge_tars.begin(); it != reg_edge_tars.end(); ++it) {
			for (auto &edge_tar : it->second) {
				cout << it->first << '\t' << edge_tar.target << '\t' << edge_tar.mi << endl;
			}
		}
	}

	//cout.rdbuf(cout_buff);
}

/* timing funcs
#include <chrono>
auto start = chrono::high_resolution_clock::now();
auto end = chrono::high_resolution_clock::now();
cout << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms" << endl;
 */

#include "ARACNe3.hpp"
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>

using namespace std;

/*
 * Assumes that simply the path to the regulator list and the path to the gene
 * expression matrix are commandline arguments
 *
 * e.g. ./ARACNe3 test/regfile.txt test/matrixfile.txt
 */
int main(int argc, char *argv[]) {
	bool multithread = false;
	bool compressed = false;

	vector<string> regs;
	genemap matrix;

	if (compressed == true) {
		regs = readRegList(string(argv[1]));
		matrix = readTransformedGexpMatrix(string(argv[2]));
	} else if (compressed == false) {
		regs = readRegList(string(argv[1]));
		matrix = readTransformedGexpMatrix(string(argv[2]));
	}

	/*
	 * Multithreading in this application is dynamic and will handle
	 * computation for each edge MI in parallel.  This is the recommended
	 * way to run ARACNe3.
	 */
	if (multithread == true) {
		if (mkdir("output", 0777) != 0) return 1;
		//cout << "REGULATOR\tTARGET\tMI\n";
		for (auto &reg : regs) {
			pid_t pid = fork();
			if (pid == 0) {
				// makes the regulator name the name of the file
				ofstream ofs{"output/" + reg + ".txt"};
				auto cout_buff = cout.rdbuf();
				cout.rdbuf(ofs.rdbuf());
				
				genemapAPMI(matrix, reg, 7.815, 4);

				cout.rdbuf(cout_buff);
				return 0;
			}
		}
	} else if (multithread == false) {
		ofstream ofs{"output.txt"};
		auto cout_buff = cout.rdbuf();
		cout.rdbuf(ofs.rdbuf());

		cout << "REGULATOR\tTARGET\tMI" << endl;
		for (auto &reg : regs) {
			genemapAPMI(matrix, reg, 7.815, 4);
		}
		cout.rdbuf(cout_buff);
	}
}

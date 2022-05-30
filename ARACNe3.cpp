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
	const vector<string> regs = readRegList(string(argv[1]));
	genemap matrix = readTransformedGexpMatrix(string(argv[2]));

	ofstream ofs{"output.txt"};
	auto cout_buff = cout.rdbuf();
	cout.rdbuf(ofs.rdbuf());

	cout << "REGULATOR\tTARGET\tMI" << endl;
	for (auto &reg : regs) {
		genemapAPMI(matrix, reg, 7.815, 4);
	}
	cout.rdbuf(cout_buff);
}

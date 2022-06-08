#include "ARACNe3.hpp"

using namespace std;

/*
 These variables are inferred while reading text files.
 */
uint32_t size_of_network = 0;
uint16_t tot_num_samps = 0;
uint16_t tot_num_regulators = 0;

/*
 These variables are tuned according to user preferences.  Some of these the user doesn't choose, such as the cached_dir, which is always the working directory of the ARACNe3 script.
 */
bool prune_FDR = true;
float FDR = 0.05f;
double subsampling_percent = 1 - std::exp(-1);
bool prune_MaxEnt = true;
bool verbose = true;
std::string cached_dir;
uint32_t global_seed = 0;
float EXPERIMENTAL_mi_cutoff = 0;

/*
 Convenient function for timing parts of ARACNe3.  It's only used to time from the pipeline function, so it's included in ARACNe3.cpp.
 */
static auto last = std::chrono::high_resolution_clock::now(), cur = std::chrono::high_resolution_clock::now();
static void sinceLast() {
	cur = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(cur-last).count() << "ms" << std::endl;
	last = cur;
}

/*
 This function is the ARACNe3 main pipeline, called from main().  The main function just parses command line arguments and options, and it sets global variables, before calling the ARACNe3 function here.
 */
void ARACNe3(string normalized_exp_mat_tsv_filename = "exp_mat.txt", string newline_separated_regulator_list_file = "regulators.txt", string output_dir = "output", double subsampling_percent = (1 - std::exp(-1)), bool prune_FDR = true, float FDR = 0.05f, bool prune_MaxEnt = true, bool verbose = true) {
	if (verbose) {
		// we don't print "... TIME:" here because we have prints in the function below
		//-------time module-------
		last = chrono::high_resolution_clock::now();
		//-------------------------
	}
	
	readRegList(newline_separated_regulator_list_file);
	genemap matrix = readExpMatrix(normalized_exp_mat_tsv_filename, subsampling_percent);
	
	if (verbose) {
		//-------time module-------
		std::cout << "\nMATRIX & REGULATORS READ TIME:" << std::endl;
		sinceLast();
		//-------------------------
		
		//-------time module-------
		cout << "\nNULL MI MODEL TIME:" << endl;
		last = chrono::high_resolution_clock::now();
		//-------------------------
	}
	
	initNullMIs(tot_num_samps);
	
	if (verbose) {
		//-------time module-------
		sinceLast();
		//-------------------------
	
		//-------time module-------
		cout << "\nRAW NETWORK COMPUTATION TIME:" << endl;
		last = chrono::high_resolution_clock::now();
		//-------------------------
	}
	
	reg_web network;
	network.reserve(tot_num_regulators);
	for (reg_id_t reg = 0; reg < tot_num_regulators; ++reg) {
		network[reg] = genemapAPMI(matrix, reg, 7.815, 4);
		
	}
	
	if (verbose) {
		//-------time module-------
		sinceLast();
		cout << "SIZE OF NETWORK: " << size_of_network << " EDGES." << endl;
		//-------------------------
	}
	
	if (1 /*pruneFDR, but we always pruneFDR*/) {
		if (!prune_FDR) FDR = 1.01f; // we must set to 1.01f to preserve all edges; rounding issue.
		if (verbose) {
			//-------time module-------
			cout << "\nFDR PRUNING TIME:" << endl;
			last = chrono::high_resolution_clock::now();
			//-------------------------
		}
		
		/*
		 We could prune in-network, but that would require many search operations.  It is better to extract edges and reform the entire network, then free memory, it seems.
		 */
		network = pruneFDR(network, size_of_network, FDR);
		
		if (verbose) {
			//-------time module-------
			sinceLast();
			cout << "SIZE OF NETWORK: " << size_of_network << " EDGES." << endl;
			//-------------------------
		}
		
		if (prune_MaxEnt) {
			// you _must_ prune FDR to do MaxEnt, but you can always FDR = 1.00
			if (verbose) {
				//-------time module-------
				std::cout << "\nMaxEnt PRUNING TIME:" << std::endl;
				last = std::chrono::high_resolution_clock::now();
				//-------------------------
			}

			
			network = pruneMaxEnt(network);
			
			if (verbose) {
				//-------time module-------
				sinceLast();
				std::cout << "SIZE OF NETWORK: " << size_of_network << " EDGES." << std::endl;
				//-------------------------
			}
		}
	}
	
	if (verbose) {
		//-------time module-------
		cout << "\nPRINTING NETWORK!" << endl;
		last = chrono::high_resolution_clock::now();
		//-------------------------
	}
	
	writeNetworkRegTarMI(network, output_dir, std::to_string(size_of_network));
	
	if (verbose) {
		//-------time module-------
		sinceLast();
		//-------------------------
		std::string success_A3 = "Success!";
		cout << success_A3 << endl;
	}
}



//--------------------cmd line parser------------------------

char* getCmdOption(char **begin, char **end, const std::string & option)
{
    char **itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
	return *itr;
    }
    return 0;
}

bool cmdOptionExists(char **begin, char **end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

//-----------------------------------------------------------


/*
 Main function is the command line executable; this primes the global variables and parses the command line.  It will also return usage notes if the user incorrectly calls ./ARACNe3.
 
 Example:
 ./ARACNe3 -e test/matrix.txt -r test/regulators.txt -o test/output --noFDR --FDR 0.05 --noMaxEnt --subsample 0.6321 --noverbose --seed 1
 */
int main(int argc, char *argv[]) {
	if (cmdOptionExists(argv, argv+argc, "-h") || cmdOptionExists(argv, argv+argc, "--help") || !cmdOptionExists(argv, argv+argc, "-e") || !cmdOptionExists(argv, argv+argc, "-r") || !cmdOptionExists(argv, argv+argc, "-o")) {
		cout << "usage: " + ((string) argv[0]) + " -e path/to/matrix.txt -r path/to/regulators.txt -o path/to/output/directory --FDR 0.05" << endl;
		return 1;
	}
	
	//--------------------cmd line parsing------------------------
	
	string matrix = (string) getCmdOption(argv, argv+argc, "-e");
	string regulators = (string) getCmdOption(argv, argv+argc, "-r");
	string output_dir = (string) getCmdOption(argv, argv+argc, "-o");
	
	// make sure output_dir has a trailing slash
	if (output_dir.back() != '/')
	    output_dir += '/';
	
	if (cmdOptionExists(argv, argv+argc, "--FDR"))
		FDR = stof(getCmdOption(argv, argv+argc, "--FDR"));
	if (FDR >= 1.00f || FDR <= 0) {
		std::cout << "FDR not on range [0,1], setting to 1.00" << std::endl;
		FDR = 1.01f;
	}
	
	if (cmdOptionExists(argv, argv+argc, "--seed"))
		global_seed = stoi(getCmdOption(argv, argv+argc, "--seed"));
	
	if (cmdOptionExists(argv, argv+argc, "--subsample"))
		subsampling_percent = stod(getCmdOption(argv, argv+argc, "--subsample"));
	
	if (subsampling_percent >= 1.00 || subsampling_percent <= 0)
		{
			std::cout << "Subsampling percent not on range [0,1]; setting to 1.00." << std::endl;
			subsampling_percent = 1.00;
	}

	if (cmdOptionExists(argv, argv+argc, "--noFDR"))
	    prune_FDR = false;
	if (cmdOptionExists(argv, argv+argc, "--noMaxEnt"))
	    prune_MaxEnt = false;
	if (cmdOptionExists(argv, argv+argc, "--noverbose"))
	    verbose = false;
	//----------------------EXPERIMENTAL--------------------------
	
	if (cmdOptionExists(argv, argv+argc, "--mithresh"))
		EXPERIMENTAL_mi_cutoff = stof(getCmdOption(argv, argv+argc, "--mithresh"));
	if (EXPERIMENTAL_mi_cutoff < 0)
		EXPERIMENTAL_mi_cutoff = 0.0f;
	
	//------------------------------------------------------------
	
	//------------------------------------------------------------
	

	cached_dir = "./"+ hiddenfpre + "ARACNe3_cached/";
	
	makeOutputDir(output_dir);
	makeCachedDir(cached_dir);
	
	ARACNe3(matrix, regulators, output_dir, subsampling_percent, prune_FDR, FDR, prune_MaxEnt, verbose);
	return 0;
}

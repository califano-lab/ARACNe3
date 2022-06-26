#include "ARACNe3.hpp"

using namespace std;


/*
 These variables are tuned according to user preferences.  Some of these the user doesn't choose, such as the cached_dir, which is always the working directory of the ARACNe3 script.
 */
bool prune_alpha = true;
bool adaptive = false;
float alpha = 0.05f;
double subsampling_percent = 1 - std::exp(-1);
bool prune_MaxEnt = true;
std::string cached_dir;
std::string output_dir;
std::string log_dir;
std::string subnets_dir;
std::string method = "FDR";
float DEVELOPER_mi_cutoff = 0;
uint16_t num_subnets = 1;

uint32_t global_seed = 0;

/*
 These variables represent the original data and do not change after matrix files are read.
 */
extern uint16_t tot_num_samps;
extern uint16_t tot_num_subsample;
extern uint16_t tot_num_regulators;
extern genemap global_gm;

extern uint32_t num_null_marginals;

/*
 Convenient function for timing parts of ARACNe3.  It's only used to time from the pipeline function, so it's included in ARACNe3.cpp.
 */
static void sinceLast(decltype(std::chrono::high_resolution_clock::now()) &last, std::ostream &ostream) {
	auto cur = std::chrono::high_resolution_clock::now();
	ostream << std::chrono::duration_cast<std::chrono::milliseconds>(cur-last).count() << "ms" << std::endl;
	last = cur;
}

/*
 This function is the ARACNe3 main pipeline, called from main().  The main function just parses command line arguments and options, and it sets global variables, before calling the ARACNe3 function here.
 */
reg_web ARACNe3_subnet(genemap& subnet_matrix, uint16_t subnet_idx) {
	auto last = std::chrono::high_resolution_clock::now();
	
	// set the individual subnet log file
	std::ofstream log_output(log_dir + "log_subnet" + std::to_string(subnet_idx) + ".txt");
	
	/*
	 Log file header
	 */
	std::time_t t = std::time(nullptr);
	log_output << "---------" << std::put_time(std::localtime(&t), "%c %Z") << "---------\n" << std::endl;
	log_output << "Subnetwork #: " + std::to_string(subnet_idx) << std::endl;
	log_output << "Total # regulators: " + std::to_string(tot_num_regulators) << std::endl;
	log_output << "Total # targets: " + std::to_string(subnet_matrix.size()) << std::endl;
	log_output << "Total # samples: " + std::to_string(tot_num_samps) << std::endl;
	log_output << "Subsampled quantity: " + std::to_string(tot_num_subsample) << std::endl;
	log_output << "Total possible edges: " + std::to_string(tot_num_regulators*subnet_matrix.size()-tot_num_regulators) << std::endl;
	log_output << "Method of first pruning step: " + method << std::endl;
	log_output << "Alpha: " + std::to_string(alpha) << std::endl;
	log_output << "\n-----------Begin Network Generation-----------\n" << std::endl;
	
	/*
	 Begin Network computation
	 */
	//-------time module-------
	log_output << "\nRAW NETWORK COMPUTATION TIME:" << std::endl;
	last = std::chrono::high_resolution_clock::now();
	//-------------------------
	
	uint32_t size_of_network = 0;
	reg_web network;
	network.reserve(tot_num_regulators);
	for (gene_id_t reg = 0; reg < tot_num_regulators; ++reg) {
		network[reg] = genemapAPMI(subnet_matrix, reg, 7.815, 4);
		size_of_network += network[reg].size();
	}
	
	//-------time module-------
	sinceLast(last, log_output);
	log_output << "SIZE OF NETWORK: " << size_of_network << " EDGES." << std::endl;
	//-------------------------
	
	if (!prune_alpha) alpha = 1.01f; // we must set to 1.01f to preserve all edges; rounding issue.
	
	//-------time module-------
	log_output << "\nALPHA PRUNING TIME (" + method + "): " << std::endl;
	last = std::chrono::high_resolution_clock::now();
	//-------------------------
	
	/*
	 We could prune in-network, but that would require many search operations.  It is better to extract edges and reform the entire network, then free memory, it seems.
	 */
	std::pair<reg_web, map_map> pair = pruneAlpha(network, size_of_network);
	network = pair.first;
	map_map& tftfNetwork = pair.second;
	
	//-------time module-------
	sinceLast(last, log_output);
	log_output << "SIZE OF NETWORK: " << size_of_network << " EDGES." << std::endl;
	//-------------------------
	
	if (prune_MaxEnt) {
		//-------time module-------
		log_output << "\nMaxEnt PRUNING TIME:" << std::endl;
		last = std::chrono::high_resolution_clock::now();
		//-------------------------

		
		network = pruneMaxEnt(network, tftfNetwork, size_of_network);
		
		//-------time module-------
		sinceLast(last, log_output);
		log_output << "SIZE OF NETWORK: " << size_of_network << " EDGES." << std::endl;
		//-------------------------
	}
	
	//-------time module-------
	log_output << "\nPRINTING NETWORK IN DIRECTORY \"" + output_dir + "\"....." << std::endl;
	last = chrono::high_resolution_clock::now();
	//-------------------------
	
	// writes the individual subnet output
	writeNetworkRegTarMI(network, subnets_dir, "subnet" + std::to_string(subnet_idx));
	
	//-------time module-------
	sinceLast(last, log_output);
	//-------------------------
	
	return network;
}


/* This is the consolidation step.  Takes the subnetworks generated, 
 */


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
 ./ARACNe3 -e test/matrix.txt -r test/regulators.txt -o test/output --noAlpha -a 0.05 --alpha 0.05 --noMaxEnt --subsample 0.6321 --seed 1 --mithresh 0.2 --numnulls 1000000
 */
int main(int argc, char *argv[]) {
	auto last = std::chrono::high_resolution_clock::now();
	
	if (cmdOptionExists(argv, argv+argc, "-h") || cmdOptionExists(argv, argv+argc, "--help") || !cmdOptionExists(argv, argv+argc, "-e") || !cmdOptionExists(argv, argv+argc, "-r") || !cmdOptionExists(argv, argv+argc, "-o")) {
		std::cout << "usage: " + ((std::string) argv[0]) + " -e path/to/matrix.txt -r path/to/regulators.txt -o path/to/output/directory --alpha 0.05" << std::endl;
		return 1;
	}
	
	//--------------------cmd line parsing------------------------
	
	std::string exp_file = (std::string) getCmdOption(argv, argv+argc, "-e");
	std::string reg_file = (std::string) getCmdOption(argv, argv+argc, "-r");
	
	output_dir = (std::string) getCmdOption(argv, argv+argc, "-o");
	
	// make sure output_dir has a trailing slash
	if (output_dir.back() != '/')
	    output_dir += '/';
	
	if (cmdOptionExists(argv, argv+argc, "--alpha"))
		alpha = stof(getCmdOption(argv, argv+argc, "--alpha"));
	if (alpha >= 1.00f || alpha <= 0) {
		std::cout << "alpha not on range [0,1], setting to 1.00" << std::endl;
		alpha = 1.01f;
	}
	
	if (cmdOptionExists(argv, argv+argc, "--seed"))
		global_seed = stoi(getCmdOption(argv, argv+argc, "--seed"));
	
	if (cmdOptionExists(argv, argv+argc, "--subsample"))
		subsampling_percent = stod(getCmdOption(argv, argv+argc, "--subsample"));
	
	if (subsampling_percent > 1.0000001 || subsampling_percent <= 0) {
			std::cout << "Subsampling percent not on range (0,1]; setting to 1.00." << std::endl;
			subsampling_percent = 1.00;
	}
	
	if (cmdOptionExists(argv, argv+argc, "-x"))
		num_subnets = stoi(getCmdOption(argv, argv+argc, "-x"));

	if (cmdOptionExists(argv, argv+argc, "--numNetworks"))
		num_subnets = stoi(getCmdOption(argv, argv+argc, "--numNetworks"));

	if (cmdOptionExists(argv, argv+argc, "--noAlpha"))
	    	prune_alpha = false;
	if (cmdOptionExists(argv, argv+argc, "--noMaxEnt"))
	    	prune_MaxEnt = false;
	if (cmdOptionExists(argv, argv+argc, "--FDR"))
		method = "FDR";
	if (cmdOptionExists(argv, argv+argc, "--FWER"))
		method = "FWER";
	if (cmdOptionExists(argv, argv+argc, "--adaptive"))
		adaptive = true;

	//----------------------DEVELOPER--------------------------
	
	if (cmdOptionExists(argv, argv+argc, "--mithresh"))
		DEVELOPER_mi_cutoff = stof(getCmdOption(argv, argv+argc, "--mithresh"));
	if (DEVELOPER_mi_cutoff < 0)
		DEVELOPER_mi_cutoff = 0.0f;
	
	if (cmdOptionExists(argv, argv+argc, "--numnulls"))
		num_null_marginals = stoi(getCmdOption(argv, argv+argc, "--numnulls"));
	if (num_null_marginals < 0) {
		std::cout << "Number of null marginals not on range (0,inf); setting to 1000000." << std::endl;
		num_null_marginals = 1000000;
	}
	
	//------------------------------------------------------------

	cached_dir = "./"+ hiddenfpre + "ARACNe3_cached/";
	
	makeDir(output_dir);
	makeDir(cached_dir);
	
	log_dir = output_dir + "log/";
	makeDir(log_dir);
	
	subnets_dir = output_dir + "subnets/";
	makeDir(subnets_dir);
	
	//-------time module-------
	last = chrono::high_resolution_clock::now();
	//-------------------------

	std::ofstream log_output(output_dir + "finalLog.txt");
	std::time_t t = std::time(nullptr);
	log_output << "---------" << std::put_time(std::localtime(&t), "%c %Z") << "---------\n" << std::endl;
	std::cout << "\nBeginning ARACNe3 network generation.  See logs and progress reports in \"" + output_dir + "finalLog.txt\"" << std::endl;
	log_output << "\nBeginning ARACNe3 network generation..." << std::endl;
	
	readRegList(reg_file);
	
	readExpMatrix(exp_file);
	
	//-------time module-------
	log_output << "\nMATRIX & REGULATORS READ TIME:" << std::endl;
	sinceLast(last, log_output);
	//-------------------------
	
	//-------time module-------
	log_output << "\nNULL MI MODEL TIME:" << endl;
	last = chrono::high_resolution_clock::now();
	//-------------------------
	
	initNullMIs(tot_num_subsample);
	
	//-------time module-------
	sinceLast(last, log_output);
	//-------------------------
	
	//-------time module-------
	log_output << "\nCREATING SUB-NETWORK(s) TIME: " << std::endl;
	//-------------------------
	
	std::vector<reg_web> subnets;
	if (adaptive) {
		subnets = std::vector<reg_web>(num_subnets);
		for (uint16_t i = 0; i < num_subnets; ++i) {
			genemap subnet_matrix = sampleFromGlobalGenemap();
			subnets.push_back(ARACNe3_subnet(subnet_matrix, i));
		}
	} else {
		subnets = std::vector<reg_web>(num_subnets);
		for (uint16_t i = 0; i < num_subnets; ++i) {
			genemap subnet_matrix = sampleFromGlobalGenemap(); 
			subnets[i] = ARACNe3_subnet(subnet_matrix, i);
		}
	}
	
	//-------time module-------
	sinceLast(last, log_output);
	//-------------------------
	
	//-------time module-------
	log_output << "\nCONSOLIDATING SUB-NETWORK(s) TIME: " << std::endl;
	//-------------------------
	
	std::vector<consolidated_df> final_df = consolidate(subnets);
	
	//-------time module-------
	sinceLast(last, log_output);
	//-------------------------
	
	//-------time module-------
	log_output << "\nWRITING FINAL NETWORK..." << std::endl;
	//-------------------------
	
	writeConsolidatedNetwork(final_df, output_dir);
	
		using namespace std::string_literals;
		const char* success_A3 =
R"(

                |
                |
                ;                            
                ;                            
                |                            
           ,  / | \   ,
         , ;_/ ,L-, `_;  ,
         \._/.ARACNe3.\_./
           \_./(::)\._/                      
                ''


SUCCESS!
)";
		std::cout << success_A3 << std::endl;
		log_output << success_A3 << std::endl;
	
	return 0;
}

#include "ARACNe3.hpp"
#include "cmdline_parser.hpp"
#include "stopwatch.hpp"
#include "subnet_operations.hpp"
#include "io.hpp"

/*
 These variables are tuned according to user preferences.  Some of these the user doesn't choose, such as the cached_dir, which is always the working directory of the ARACNe3 script.
 */
bool prune_alpha = true;
bool adaptive = false;
bool do_not_consolidate = false;
bool go_to_consolidate = false;
float alpha = 0.05f;
double subsampling_percent = 1 - std::exp(-1);
bool prune_MaxEnt = true;
std::string cached_dir;
std::string output_dir;
std::string log_dir;
std::string subnets_dir;
std::string method = "FDR";
float DEVELOPER_mi_cutoff = 0.0f;
uint16_t num_subnets = 1U;
uint16_t num_subnets_to_consolidate = 0U;
uint16_t targets_per_regulator = 30U;
uint16_t nthreads = 1U;

uint32_t seed = 0U;

/*
 These variables represent the original data and do not change after matrix files are read.
 */
extern std::vector<std::string> decompression_map;

extern uint32_t num_null_marginals;

extern float FPR_estimate;
extern std::vector<float> FPR_estimates;

/*
 Main function is the command line executable; this primes the global variables and parses the command line.  It will also return usage notes if the user incorrectly calls ./ARACNe3.
 
 Example:
 ./ARACNe3 -e test/matrix.txt -r test/regulators.txt -o test/output --noAlpha -a 0.05 --alpha 0.05 --noMaxEnt --subsample 0.6321 --seed 1 --mithresh 0.2 --numnulls 1000000
 */
int main(int argc, char *argv[]) {
        if (cmdOptionExists(argv, argv + argc, "-h") ||
            cmdOptionExists(argv, argv + argc, "--help") ||
            !cmdOptionExists(argv, argv + argc, "-e") ||
            !cmdOptionExists(argv, argv + argc, "-r") ||
            !cmdOptionExists(argv, argv + argc, "-o")) {
          std::cout
              << "usage: " + ((std::string)argv[0]) +
                     makeUnixDirectoryNameUniversal(
                         " -e path/to/matrix.txt -r path/to/regulators.txt -o "
                         "path/to/output/directory")
              << std::endl;
          return EXIT_FAILURE;
        }

        //--------------------cmd line parsing------------------------
	
	const std::string exp_mat_file = makeUnixDirectoryNameUniversal((std::string) getCmdOption(argv, argv+argc, "-e"));
	const std::string reg_list_file = makeUnixDirectoryNameUniversal((std::string) getCmdOption(argv, argv+argc, "-r"));
	
	output_dir = (std::string) getCmdOption(argv, argv+argc, "-o");
	
	// make sure output_dir has a trailing slash
	if (output_dir.back() != directory_slash)
	    output_dir += directory_slash;
	
	if (cmdOptionExists(argv, argv+argc, "--alpha"))
		alpha = std::stof(getCmdOption(argv, argv+argc, "--alpha"));
	if (alpha >= 1.00f || alpha <= 0) {
		std::cout << "alpha not on range [0,1], setting to 1.00" << std::endl;
		alpha = 1.01f;
	}
	
	if (cmdOptionExists(argv, argv+argc, "--seed"))
		seed = std::stoi(getCmdOption(argv, argv+argc, "--seed"));
	
	if (cmdOptionExists(argv, argv+argc, "--subsample"))
		subsampling_percent = std::stod(getCmdOption(argv, argv+argc, "--subsample"));
	
	if (subsampling_percent > 1.0000001 || subsampling_percent <= 0) {
			std::cout << "Subsampling percent not on range (0,1]; setting to 1.00." << std::endl;
			subsampling_percent = 1.00;
	}
	
	if (cmdOptionExists(argv, argv+argc, "-x")) 
		num_subnets = targets_per_regulator = num_subnets_to_consolidate = std::stoi(getCmdOption(argv, argv+argc, "-x"));

	if (cmdOptionExists(argv, argv+argc, "--threads"))
		nthreads = std::stoi(getCmdOption(argv, argv+argc, "--threads"));

	if (cmdOptionExists(argv, argv+argc, "--noAlpha"))
	    	prune_alpha = false;
        alpha = 1.f;
	
	if (cmdOptionExists(argv, argv+argc, "--noMaxEnt"))
	    	prune_MaxEnt = false;
	if (cmdOptionExists(argv, argv+argc, "--FDR"))
		method = "FDR";
	if (cmdOptionExists(argv, argv+argc, "--FWER"))
		method = "FWER";
	if (cmdOptionExists(argv, argv+argc, "--FPR"))
		method = "FPR";
	if (cmdOptionExists(argv, argv+argc, "--adaptive"))
		adaptive = true;
	if (cmdOptionExists(argv, argv+argc, "--noconsolidate"))
		do_not_consolidate = true;
	if (cmdOptionExists(argv, argv+argc, "--consolidate"))
		go_to_consolidate = true;

	//----------------------DEVELOPER--------------------------
	
	if (cmdOptionExists(argv, argv+argc, "--mithresh"))
		DEVELOPER_mi_cutoff = std::stof(getCmdOption(argv, argv+argc, "--mithresh"));
	if (DEVELOPER_mi_cutoff < 0)
		DEVELOPER_mi_cutoff = 0.0f;
	
	if (cmdOptionExists(argv, argv+argc, "--numnulls"))
		num_null_marginals = std::stoi(getCmdOption(argv, argv+argc, "--numnulls"));
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
	
	std::ofstream log_output(output_dir + "finalLog.txt");
	
	// print the initial command to the log output
	for (uint16_t i = 0; i < argc; ++i)
		log_output << std::string(argv[i]) << " ";
	log_output << std::endl;
	
	std::time_t t = std::time(nullptr);
	std::cout << "\n---------" << std::put_time(std::localtime(&t), "%c %Z") << "---------" << std::endl;
	log_output << "\n---------" << std::put_time(std::localtime(&t), "%c %Z") << "---------" << std::endl;
	
	std::cout << "Beginning ARACNe3 instance.  See logs and progress reports in \"" + makeUnixDirectoryNameUniversal(output_dir) + "finalLog.txt\"." << std::endl;
	log_output << "Beginning ARACNe3 instance..." << std::endl;

	std::mt19937 rand{seed};
  Watch watch1;
  watch1.reset();
	
  log_output << "\nGene expression matrix & regulators list read time: ";

  const auto [exp_mat, ranks_mat, regulators, tot_num_samps] = readExpMatrixAndCopulaTransform(exp_mat_file, subsampling_percent, rand);

	uint16_t tot_num_subsample = std::ceil(subsampling_percent * tot_num_samps);
	if (tot_num_subsample >= tot_num_samps || tot_num_subsample < 0) {
		std::cerr << "Warning: subsample quantity invalid. All samples will be used." << std::endl;
		tot_num_subsample = tot_num_samps;
	}

	std::cout << "\nTotal N Samples: " + std::to_string(tot_num_samps) << std::endl;
	std::cout << "Subsampled N Samples: " + std::to_string(tot_num_subsample) << std::endl;

  readRegList(reg_list_file);

  //-------time module-------
	log_output << watch1.getSeconds() << std::endl;
  log_output << "\nMutual Information null model calculation time: ";
  watch1.reset();
	//-------------------------
	
	initNullMIs(tot_num_subsample);
	
	//-------time module-------
	log_output << watch1.getSeconds() << std::endl;
	//-------------------------
	
	// Must exist regardless of whether we skip to consolidation
	std::vector<gene_to_edge_tars> subnets;
	if(!go_to_consolidate) {

		//-------time module-------
		log_output << "\nCreating subnetwork(s) time: ";
    watch1.reset();
		//-------------------------
		
		if (adaptive) {
			bool stoppingCriteriaMet = false;
			gene_to_geneset regulon_set;
			for (const uint16_t &reg : regulators) regulon_set[reg];
			uint16_t i = 0U;
			while (!stoppingCriteriaMet) {
				gene_to_floats subnet_matrix = sampleExpMatAndReCopulaTransform(exp_mat, rand);
				subnets.push_back(ARACNe3_subnet(subnet_matrix, i+1, prune_alpha, method, alpha, prune_MaxEnt, output_dir, subnets_dir, log_dir, nthreads));
				
				// add any new edges to the regulon_set
				for (const auto &[reg, edge_tars] : subnets[i])
					for (const auto &edge : edge_tars)
						regulon_set[reg].insert(edge.target);
				
				// check stoping criteria
				uint16_t min = 65535U;
				for (const auto &[reg, regulon] : regulon_set)
					if (exp_mat.find(reg) != exp_mat.end())
						if (regulon.size() < min)
							min = regulon.size();
				if (min >= targets_per_regulator) 
					stoppingCriteriaMet = true;
				
				++i;
			}
			num_subnets = i;
		} else {
			if (nthreads > 1 && num_subnets > 1)
				std::cout << "Note: Because more than one thread is used to compute a fixed number of subnetworks (--adaptive not specified), completion times may not be in order.  This is not an error.\n" << std::endl;
			subnets = std::vector<gene_to_edge_tars>(num_subnets);
			for (int i = 0; i < num_subnets; ++i) {
				gene_to_floats subnet_matrix = sampleExpMatAndReCopulaTransform(exp_mat, rand);
				subnets[i] = ARACNe3_subnet(subnet_matrix, i+1, prune_alpha, method, alpha, prune_MaxEnt, output_dir, subnets_dir, log_dir, nthreads);
			}
		}
		
		//-------time module-------
    log_output << watch1.getSeconds() << std::endl;
		//-------------------------
    
		log_output << "Total subnetworks generated: " + std::to_string(num_subnets) << std::endl;
	} else if (go_to_consolidate) {

		//-------time module-------
		log_output << "\nReading subnetwork(s) time: ";
    watch1.reset();
		//-------------------------
		
		for (uint16_t subnet_num = 1; subnet_num <= num_subnets_to_consolidate; ++subnet_num) {
			try {
				subnets.push_back(readSubNetAndUpdateFPRFromLog(output_dir, subnet_num));
				num_subnets = subnet_num;
			} catch (TooManySubnetsRequested e) {
				std::cout << "WARNING: " + std::string(e.what()) << std::endl;
				break;
			}
		}
		
		//-------time module-------
		log_output << watch1.getSeconds() << std::endl;
		//-------------------------

		log_output << "Total subnets read: " + std::to_string(num_subnets) << std::endl;
	}
	
	// set the FPR estimate
	FPR_estimate = std::accumulate(FPR_estimates.begin(), FPR_estimates.end(), 0.0f) / FPR_estimates.size();
	
	if (!do_not_consolidate) {	
		//-------time module-------
		log_output << "\nConsolidating subnetwork(s) time: ";
    watch1.reset();
		//-------------------------
		
		std::vector<consolidated_df_row> final_df = consolidate_subnets_vec(subnets);
		
		//-------time module-------
    log_output << watch1.getSeconds() << std::endl;
		log_output << "\nWriting final network..." << std::endl;
		//-------------------------
		
		writeConsolidatedNetwork(final_df, output_dir + "finalNet_" + std::to_string(num_subnets) + "subnets.txt");
		
		/* Now that we know how many subnets, we can rename finalLog to include that. */
		std::string final_log_newname = "finalLog_" + std::to_string(num_subnets) + "subnets-consolidate.txt";

		//-------time module-------
		log_output << std::endl << "Renaming \"finalLog.txt\" to \"" + final_log_newname + "\"..." << std::endl;
		std::cout << std::endl << "Renaming \"finalLog.txt\" to \"" + final_log_newname + "\"..." << std::endl;
		//-------------------------
    
		std::filesystem::rename(output_dir + "finalLog.txt", output_dir  + final_log_newname);
		
	} else if (do_not_consolidate) {

		//-------time module-------
		log_output << "\nNo consolidation requested." << std::endl;
		//-------------------------
		
		/* Now that we know how many subnets, we can rename finalLog to include that. */
		std::string final_log_newname = "finalLog_" + std::to_string(num_subnets) + "subnets-noconsolidate.txt";

		//-------time module-------
		log_output << std::endl << "Renaming \"finalLog.txt\" to \"" + final_log_newname + "\"..." << std::endl;
		std::cout << std::endl << "Renaming \"finalLog.txt\" to \"" + final_log_newname + "\"..." << std::endl;
		//-------------------------

		std::filesystem::rename(output_dir + "finalLog.txt", output_dir  + final_log_newname);
	}
	
	
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
	
	return EXIT_SUCCESS;
}

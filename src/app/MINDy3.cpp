#include "MINDy3.hpp"
#include "io.hpp"
#include "algorithms.hpp"
#include "apcmi_nullmodel.hpp"

extern uint16_t tot_num_samps;
extern std::set<gene_id> regulators, targets, modulators, genes;
extern std::vector<std::string> decompression_map;

//--------------------timer----------------------------------

static std::chrono::time_point<std::chrono::high_resolution_clock> watch1 = std::chrono::high_resolution_clock::now(), watch2 = watch1;

void setWatch(std::chrono::time_point<std::chrono::high_resolution_clock> &watch_to_reset = watch1) {
	watch_to_reset = std::chrono::high_resolution_clock::now();
}

void printWatch(std::ostream &ofs, std::string message = "", std::chrono::time_point<std::chrono::high_resolution_clock> &zero = watch1) {
	auto cur = std::chrono::high_resolution_clock::now();
	ofs << message;
	ofs << std::chrono::duration_cast<std::chrono::seconds>(cur-zero).count() << "s" << std::endl;
}

//-----------------------------------------------------------

//--------------------cmd line parser------------------------

char* getCmdOption(char **begin, char **end, const std::string &option) {
	char **itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
		return *itr;
	return 0;
}

bool cmdOptionExists(char **begin, char **end, const std::string &option) {
	return std::find(begin, end, option) != end;
}

//-----------------------------------------------------------

int main(int argc, char *argv[]) {
	std::string cached_dir;
	std::string output_dir;
	uint32_t seed = 0U;
	uint16_t nthreads = 1U;
	
	bool only_save_integrated_outputs = false;
	
	if (cmdOptionExists(argv, argv+argc, "-h") || cmdOptionExists(argv, argv+argc, "--help") || !cmdOptionExists(argv, argv+argc, "-e") || !cmdOptionExists(argv, argv+argc, "-n") || !cmdOptionExists(argv, argv+argc, "-o") || !cmdOptionExists(argv, argv+argc, "-m")) {
		std::cout << "usage: " + ((std::string) argv[0]) + makeUnixDirectoryNameUniversal(" -e path/to/matrix.txt -n path/to/network.txt -m path/to/modulators.txt -o path/to/output/directory") << std::endl;
		return EXIT_FAILURE;
	}
	
	
	//--------------------cmd line parsing------------------------
	
	std::string exp_mat_file = (std::string) getCmdOption(argv, argv+argc, "-e");
	std::string network_file = (std::string) getCmdOption(argv, argv+argc, "-n");
	std::string modulators_list_file = (std::string) getCmdOption(argv, argv+argc, "-m");
	
	output_dir = (std::string) getCmdOption(argv, argv+argc, "-o");
	
	// make sure output_dir has a trailing slash
	if (output_dir.back() != directory_slash)
		output_dir += directory_slash;
	
	if (cmdOptionExists(argv, argv+argc, "--seed"))
		seed = std::stoi(getCmdOption(argv, argv+argc, "--seed"));
	
	if (cmdOptionExists(argv, argv+argc, "--threads"))
		nthreads = std::stoi(getCmdOption(argv, argv+argc, "--threads"));
	
	if (cmdOptionExists(argv, argv+argc, "--only-integrate"))
		only_save_integrated_outputs = true;
	
	//--------------------developer options------------------------
	
	uint8_t n_bins = 10U;
	if (cmdOptionExists(argv, argv+argc, "--bins"))
		n_bins = std::stoi(getCmdOption(argv, argv+argc, "--bins"));
	
	uint32_t n_nulls = 1000000U;
	if (cmdOptionExists(argv, argv+argc, "--nulls"))
		n_nulls = std::stoi(getCmdOption(argv, argv+argc, "--nulls"));
	
	//-------------------------------------------------------------
	
	cached_dir = "./"+ hiddenfpre + "MINDy3_cached/";
	
	makeDir(output_dir);
	if (!only_save_integrated_outputs)
		makeDir(output_dir + "modulator_output/");
	makeDir(output_dir + "modulator_log/");
	makeDir(output_dir + "modulator_output-integrated/");
	makeDir(cached_dir);
	
	std::ofstream log_output(output_dir + "log.txt");
	
	// print the initial command to the log output
	for (uint16_t i = 0; i < argc; ++i)
		log_output << std::string(argv[i]) << " ";
	log_output << std::endl;
	
	std::time_t t = std::time(nullptr);
	std::cout << "\n---------" << std::put_time(std::localtime(&t), "%c %Z") << "---------" << std::endl;
	log_output << "\n---------" << std::put_time(std::localtime(&t), "%c %Z") << "---------" << std::endl;
	
	std::cout << "Beginning MINDy3 instance.  See logs and progress reports in \"" + makeUnixDirectoryNameUniversal(output_dir) + "log.txt\"." << std::endl;
	log_output << "Beginning MINDy3 instance..." << std::endl;
	
	// initialize random object
	std::mt19937 rand{seed};
	
	// Read the provided files
	setWatch();
	log_output << "Files read time: " << std::flush;
	std::pair<hash_float_vec, hash_short_vec> matrices = readExpMatrixAndCopulaTransform(exp_mat_file, rand);
	hash_float_vec &gm = matrices.first;
	hash_short_vec &gm_ranks = matrices.second;
	hash_network nw = readNetworkFile(network_file);
	readModList(modulators_list_file);
	printWatch(log_output);
	
	// Instantiate and cache null model
	setWatch();
	log_output << "Null model generation/read time: " << std::flush;
	APCMINullModel nullmodel(tot_num_samps, n_bins, n_nulls, rand, cached_dir);
	nullmodel.cacheNullModel(cached_dir);
	printWatch(log_output);

	// Copy sets to vectors
	const std::vector<gene_id> regulators_vec(regulators.begin(), regulators.end()), targets_vec(targets.begin(), targets.end()), modulators_vec(modulators.begin(), modulators.end()), genes_vec(genes.begin(), genes.end());
	
	// Same regulator-target MIs are used for each modulator test.  Calculate those MIs
	setWatch();
	log_output << "Calculate all reg-tar MIs and SCC time: " << std::flush;
	hash_edge_float reg_tar_mis, reg_tar_sccs; // TODO: Can we parallelize with hashmap if container sizes are reserved?
	for (int reg_idx = 0; reg_idx < regulators_vec.size(); ++reg_idx) {
		const gene_id &reg = regulators_vec[reg_idx];
		for (int tar_idx = 0; tar_idx < targets_vec.size(); ++tar_idx) {
			const gene_id &tar = targets_vec[tar_idx];
			if (nw[reg].find(tar) != nw[reg].end() && reg != tar) {
				reg_tar_mis[reg][tar] = calcAPMI(gm[reg], gm[tar]);
				reg_tar_sccs[reg][tar] = spearmanCorrelate(gm_ranks[reg], gm_ranks[tar]);
			}
		}
	}
	printWatch(log_output);
	
	setWatch(watch1);
	log_output << "Total modulator inference time: " << std::flush;
#pragma omp parallel for firstprivate(gm, gm_ranks, nw, reg_tar_mis, reg_tar_sccs) num_threads(nthreads)
	for (std::vector<gene_id>::size_type mod_idx = 0U; mod_idx < modulators_vec.size(); ++mod_idx) {
		const gene_id &mod = modulators_vec[mod_idx];
		
		std::ofstream mod_ofs;
		if (!only_save_integrated_outputs)
			mod_ofs.open(output_dir + "modulator_output/" + decompression_map[mod] + ".txt");
		
		std::ofstream mod_ofs_integrated{output_dir + "modulator_output-integrated/" + decompression_map[mod] + ".txt"};
		std::ofstream mod_log{output_dir + "modulator_log/" + decompression_map[mod] + "_log.txt"};
		
		/*
		 Log file header
		 */
		std::time_t t = std::time(nullptr);
		mod_log << "---------" << std::put_time(std::localtime(&t), "%c %Z") << "---------" << std::endl << std::endl;
		mod_log << "Modulator: " + decompression_map[mod] << std::endl;
		mod_log << std::endl << "-----------Begin Modulator Inference-----------" << std::endl;
		setWatch(watch2);
		mod_log << "Modulator inference time for " + decompression_map[mod] + ": " << std::flush;
		
		// output header row to mod_ofs
		if (!only_save_integrated_outputs) {
			mod_ofs << 
			"regulator.values" << '\t' <<
			"target.values" << '\t' <<
			"modulator.values" << '\t' <<
			"reg-tar-mi.values" << '\t' << 
			"reg-mod-mi.values" << '\t' <<
			"mod-tar-mi.values" << '\t' <<
			"reg-tar-scc.values" << '\t' <<
			"reg-mod-scc.values" << '\t' <<
			"mod-tar-scc.values" << '\t' <<
			"apcmi.stat.values" << '\t' <<
			"apcmi.p.values" << '\t' <<
			"moa.values" << '\t' <<
			"ba.values" << '\n';
		}
		
		// output header row to mod_ofs_integrated
		mod_ofs_integrated <<
		"modulator.values" << '\t' <<
		"regulator.values" << '\t' <<
		"integrated.p.values" << '\t' <<
		"effect.size" << '\t' <<
		"directionality" << '\t' <<
		"mean.ba" << '\n';
		
		for (size_t reg_idx = 0U; reg_idx < regulators_vec.size(); ++reg_idx) {
			const gene_id &reg = regulators_vec[reg_idx];
			
			// each regulator needs integrated statistics
			std::vector<float> apcmi_p_values, moa_values, ba_values;
			apcmi_p_values.reserve(regulators_vec.size());
			moa_values.reserve(regulators_vec.size());
			ba_values.reserve(regulators_vec.size());
			
			for (size_t tar_idx = 0U; tar_idx < targets_vec.size(); ++tar_idx) {
				const gene_id &tar = targets_vec[tar_idx];
				if (nw[reg].find(tar) != nw[reg].end() && reg != tar && reg != mod && mod != tar) {
					// calculate the MI of each edge based on APMI
					const float reg_tar_mi = reg_tar_mis[reg][tar], reg_mod_mi = calcAPMI(gm[mod], gm[reg]), mod_tar_mi = calcAPMI(gm[mod], gm[tar]), reg_tar_scc = reg_tar_sccs[reg][tar], reg_mod_scc = spearmanCorrelate(gm_ranks[reg], gm_ranks[mod]), mod_tar_scc = spearmanCorrelate(gm_ranks[mod], gm_ranks[tar]);
					
					// determine which edge is the weakest
					bool reg_tar_is_weakest = reg_tar_mi <= reg_mod_mi && reg_tar_mi <= mod_tar_mi, reg_mod_is_weakest = reg_mod_mi <= reg_tar_mi && reg_mod_mi <= mod_tar_mi;
					
					// random tie breaking if there are equal mi values
					if (reg_tar_mi == reg_mod_mi && reg_mod_mi == mod_tar_mi)
						while (reg_tar_is_weakest && reg_mod_is_weakest) {
							reg_tar_is_weakest = ((int) rand()) % 2;
							reg_mod_is_weakest = ((int) rand()) % 2;
						}
					else if (reg_tar_mi == reg_mod_mi && reg_tar_is_weakest && reg_mod_is_weakest)
						while (!(reg_tar_is_weakest ^ reg_mod_is_weakest)) {
							reg_tar_is_weakest = ((int) rand()) % 2;
							reg_mod_is_weakest = ((int) rand()) % 2;
						}
					else if (reg_tar_mi == mod_tar_mi && reg_tar_is_weakest)
							reg_tar_is_weakest = ((int) rand()) % 2;
					else if (reg_mod_mi == mod_tar_mi && reg_mod_is_weakest)
							reg_mod_is_weakest = ((int) rand()) % 2;
					
					// calculate APCMI given the weakest MI
					float apcmi;
					if (reg_tar_is_weakest)
						apcmi = APCMI(gm[reg], gm[tar], gm[mod], n_bins, rand);
					else if (reg_mod_is_weakest)
						apcmi = APCMI(gm[mod], gm[reg], gm[tar], n_bins, rand);
					else
						apcmi = APCMI(gm[mod], gm[tar], gm[reg], n_bins, rand);
					
					const float cmipval = nullmodel.getCMIPVal(apcmi, 0.001f);
					const float moa = calcModeOfAction(gm[reg], gm[tar], gm[mod], n_bins, rand);
					const float bio_activity = spearmanCorrelate(gm_ranks[reg], gm_ranks[tar]) * spearmanCorrelate(gm_ranks[mod], gm_ranks[tar]);
					
					apcmi_p_values.emplace_back(cmipval);
					moa_values.emplace_back(moa);
					ba_values.emplace_back(bio_activity);
					
					if (!only_save_integrated_outputs) {
						// output row to mod_ofs
						mod_ofs << 
						decompression_map[reg] << '\t' <<
						decompression_map[tar] << '\t' <<
						decompression_map[mod] << '\t' <<
						reg_tar_mi << '\t' <<
						reg_mod_mi << '\t' <<
						mod_tar_mi << '\t' <<
						reg_tar_scc << '\t' <<
						reg_mod_scc << '\t' <<
						mod_tar_scc << '\t' <<
						apcmi << '\t' <<
						cmipval << '\t' <<
						moa << '\t' <<
						bio_activity << '\n';
					}
				}
			}
			
			mod_ofs_integrated <<
			decompression_map[mod] << '\t' <<
			decompression_map[reg] << '\t' <<
      fishersMethodP(apcmi_p_values) << '\t' <<
			1 - std::accumulate(apcmi_p_values.begin(), apcmi_p_values.end(), 0.0f)/apcmi_p_values.size() << '\t' <<
			std::accumulate(moa_values.begin(), moa_values.end(), 0.0f)/moa_values.size() << '\t' <<
			std::accumulate(ba_values.begin(), ba_values.end(), 0.0f)/ba_values.size() << '\n';
		}
		
		printWatch(mod_log, "", watch2);
	}
	printWatch(log_output, "", watch1);
	
	const char* success_M3 =
R"(
SUCCESS!
)";
		std::cout << success_M3 << std::endl;
		log_output << success_M3 << std::endl;
	
	return EXIT_SUCCESS;
}

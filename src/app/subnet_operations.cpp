/*
 subnet_operations for ARACNe3.  Contains various functions only needed in the consolidation step, such as calculation of the p-value for an edge based on the number of subnetworks it appeared in, the calculation of the SCC for an edge, etc.
 */

#include "subnet_operations.hpp"
#include "stopwatch.hpp"
#include "io.hpp"

std::vector<float> FPR_estimates;
float FPR_estimate = 1.5E-4f;

extern uint16_t tot_num_samps;
extern uint16_t tot_num_subsample;
extern uint16_t tot_num_regulators, defined_regulators;
extern gene_to_floats global_gm;
extern gene_to_shorts global_gm_r; 
extern uint16_t num_subnets;

float consolidate_scc(const std::vector<uint16_t>& x_ranked, const std::vector<uint16_t>& y_ranked) {
	const auto &n = x_ranked.size();
	int sigma_dxy = 0;
	for (uint16_t i = 0; i < n; ++i)
		sigma_dxy += (x_ranked[i] - y_ranked[i]) * (x_ranked[i] - y_ranked[i]);
	return 1 - 6.0f * sigma_dxy / n / (n * n - 1);
}

double lchoose(const uint16_t &n, const uint16_t &k) {
	return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
}

double right_tail_binomial_p(const uint16_t& num_occurrences) {
	double p = 0.0;
	for (uint16_t i = num_subnets; i >= num_occurrences; --i)
		p += std::exp(lchoose(num_subnets, i) + i * std::log(FPR_estimate) + (num_subnets - i) * std::log(1-FPR_estimate));
	return p;
}
/*
 Generates an ARACNe3 subnet (called from main).   
*/
gene_to_edge_tars ARACNe3_subnet(gene_to_floats subnet_matrix, const uint16_t subnet_num, const bool prune_alpha, const std::string& method, const float alpha, const bool prune_MaxEnt, const std::string& output_dir, const std::string& subnets_dir, const std::string& log_dir, const uint16_t nthreads) {
	std::ofstream log_output(log_dir + "log_subnet" + std::to_string(subnet_num) + ".txt");
	std::time_t t = std::time(nullptr);
	log_output << "---------" << std::put_time(std::localtime(&t), "%c %Z") << "---------" << std::endl << std::endl;
	log_output << "Subnetwork #: " + std::to_string(subnet_num) << std::endl;
	log_output << "Total # regulators (with gexp profile defined): " + std::to_string(defined_regulators) << std::endl;
	log_output << "Total # targets: " + std::to_string(subnet_matrix.size()) << std::endl;
	log_output << "Total # samples: " + std::to_string(tot_num_samps) << std::endl;
	log_output << "Subsampled quantity: " + std::to_string(tot_num_subsample) << std::endl;
	log_output << "Total possible edges: " + std::to_string(defined_regulators*subnet_matrix.size()-defined_regulators) << std::endl;
	log_output << "Method of first pruning step: " + method << std::endl;
	log_output << "Alpha: " + std::to_string(alpha) << std::endl;
	log_output << "MaxEnt Pruning: " + std::to_string(prune_MaxEnt) << std::endl;
	log_output << std::endl << "-----------Begin Network Generation-----------" << std::endl;
	
  // begin subnet computation

	//-------time module-------
  Watch watch1;
  log_output << "\nRaw network computation time: ";
  watch1.reset();
	//-------------------------
	
	uint32_t size_of_network = 0;
	std::vector<std::vector<edge_tar>> network_vec(tot_num_regulators); 
#pragma omp parallel for firstprivate(subnet_matrix) num_threads(nthreads)
	for (int reg = 0; reg < tot_num_regulators; ++reg) {
		if (global_gm.find(reg) != global_gm.end()) {
			network_vec[reg] = gene_to_floatsAPMI(subnet_matrix, reg, 7.815, 4);
			size_of_network += network_vec[reg].size();
		}
	}
	gene_to_edge_tars network;
	network.reserve(tot_num_regulators);
	for (gene_id reg = 0; reg < tot_num_regulators; ++reg)
		if (global_gm.find(reg) != global_gm.end())
			network[reg] = network_vec[reg];
	std::vector<std::vector<edge_tar>>().swap(network_vec);
	
	//-------time module-------
	log_output << watch1.getSeconds() << std::endl;
	log_output << "Size of network: " << size_of_network << " edges." << std::endl;
	//-------------------------
	
	//-------time module-------
  log_output << "\nThreshold pruning time (" + method + "): ";
  watch1.reset();
	//-------------------------
	
	auto size_prev = size_of_network;
	
	/*
	 We could prune in-network, but that would require many search operations.  It is better to extract edges and reform the entire network, then free memory, it seems.
	 */
	
	std::pair<gene_to_edge_tars, gene_to_gene_to_float> pair = pruneAlpha(network, size_of_network);
	network = pair.first;
	gene_to_gene_to_float& tftfNetwork = pair.second;
	
	//-------time module-------
	log_output << watch1.getSeconds() << std::endl;
	log_output << "Edges removed: " << size_prev - size_of_network << " edges." << std::endl;
	log_output << "Size of network: " << size_of_network << " edges." << std::endl;
	//-------------------------
	
  // save for binomial distribution parameter (theta)
	uint32_t num_edges_after_threshold_pruning = size_of_network; 
	
	if (prune_MaxEnt) {
		//-------time module-------
    log_output << "\nMaxEnt pruning time: ";
    watch1.reset();
		//-------------------------

		size_prev = size_of_network;
		network = pruneMaxEnt(network, tftfNetwork, size_of_network);
		
		//-------time module-------
		log_output << watch1.getSeconds() << std::endl;
		log_output << "Edges removed: " << size_prev - size_of_network << " edges." << std::endl;
		log_output << "Size of network: " << size_of_network << " edges." << std::endl;
		//-------------------------
		
		uint32_t num_edges_after_MaxEnt_pruning = size_of_network;
		if (method == "FDR")
			FPR_estimates.emplace_back((alpha*num_edges_after_MaxEnt_pruning)/(defined_regulators*global_gm.size()-(1-alpha)*num_edges_after_threshold_pruning));
		else if (method == "FWER")
			FPR_estimates.emplace_back((alpha/(defined_regulators*(global_gm.size()-1)))*(num_edges_after_MaxEnt_pruning)/(num_edges_after_threshold_pruning));
		else if (method == "FPR")
			FPR_estimates.emplace_back(alpha*num_edges_after_MaxEnt_pruning/num_edges_after_threshold_pruning);
	} else {
		if (method == "FDR")
			FPR_estimates.emplace_back((alpha*num_edges_after_threshold_pruning)/(defined_regulators*global_gm.size()-(1-alpha)*num_edges_after_threshold_pruning));
		else if (method == "FWER")
			FPR_estimates.emplace_back(alpha/(defined_regulators*(global_gm.size()-1)));
		else if (method == "FPR")
			FPR_estimates.emplace_back(alpha);
	}
	
	//-------time module-------
  log_output << "\nPrinting network in directory \"" + output_dir + "\".....";
  watch1.reset();
	//-------------------------
	
	// writes the individual subnet output
	writeNetworkRegTarMI(network, subnets_dir, "subnet" + std::to_string(subnet_num));
	
	//-------time module-------
	log_output << watch1.getSeconds() << std::endl;
	//-------------------------
	
	std::cout << "... subnetwork " + std::to_string(subnet_num) + " completed = " + std::to_string(size_of_network) + " edges returned ..." << std::endl;
	
	return network;
}


std::vector<consolidated_df_row> consolidate_subnets_vec(std::vector<gene_to_edge_tars> &subnets) {
	std::vector<consolidated_df_row> final_df;
	const auto tot_poss_edgs = defined_regulators*(global_gm.size()-1);
	final_df.reserve(tot_poss_edgs);
	
	std::vector<gene_to_gene_to_float> subnets_mpmp;
	for (uint16_t i = 0; i < subnets.size(); ++i)
		subnets_mpmp.emplace_back(regweb_to_mapmap(subnets[i]));
	
	for (uint16_t reg = 0; reg < tot_num_regulators; ++reg) {
		for (const auto &[tar, tar_vec] : global_gm) {
			uint16_t num_occurrences = 0;
			for (uint16_t sn = 0; sn < subnets.size(); ++sn) {
				if (subnets_mpmp[sn][reg].find(tar) != subnets_mpmp[sn][reg].end())
					++num_occurrences;
			}
			if (num_occurrences > 0) {
				const float final_mi = APMI(global_gm[reg], global_gm[tar]);
				const float final_scc = consolidate_scc(global_gm_r[reg], global_gm_r[tar]);
				const double final_p = right_tail_binomial_p(num_occurrences);
				final_df.emplace_back(reg, tar, final_mi, final_scc, num_occurrences, final_p);
			}
		}
	}
	
	return final_df;
}

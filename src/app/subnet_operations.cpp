/*
 subnet_operations for ARACNe3.  Contains various functions only needed in the consolidation step, such as calculation of the p-value for an edge based on the number of subnetworks it appeared in, the calculation of the SCC for an edge, etc.
 */

#include "subnet_operations.hpp"
#include "stopwatch.hpp"
#include "io.hpp"
#include "algorithms.hpp"

std::vector<float> FPR_estimates;
float FPR_estimate = 1.5E-4f;

extern uint16_t tot_num_samps;
extern uint16_t tot_num_subsample;
extern std::unordered_set<gene_id> regulators, genes;
extern uint32_t num_null_marginals;

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

double right_tail_binomial_p(const uint16_t &num_occurrences, const uint16_t &num_subnets) {
	double p = 0.0;
	for (uint16_t i = num_subnets; i >= num_occurrences; --i)
		p += std::exp(lchoose(num_subnets, i) + i * std::log(FPR_estimate) + (num_subnets - i) * std::log(1-FPR_estimate));
	return p;
}

/*
 Prunes a network by control of alpha using the Benjamini-Hochberg Procedure if method = FDR, or FWER if method = FWER.
 */
std::tuple<gene_to_gene_to_float, uint32_t, gene_to_gene_to_float>
pruneAlpha(const gene_to_gene_to_float &network, uint32_t size_of_network,
           const std::string &method, const float &alpha) {

  /* A vector that describes each regulator-mi-target interaction must be
   * initialized for sorting-based pruning */
  std::vector<const std::tuple<gene_id, gene_id, float>> reg_tar_mi;
  reg_tar_mi.reserve(size_of_network);

  for (gene_id reg : regulators)
    for (const auto &[tar, mi] : network.at(reg))
      reg_tar_mi.emplace_back(reg, tar, mi);

  // sort descending
  std::sort(reg_tar_mi.begin(), reg_tar_mi.end(),
            [](const std::tuple<gene_id, gene_id, float> &rtm1,
               const std::tuple<gene_id, gene_id, float> &rtm2) -> bool {
              return std::get<2>(rtm1) > std::get<2>(rtm2);
            });

  uint32_t argmax_k = 0;
  uint32_t m = size_of_network;
  if (method == "FDR") {
    // Benjamini-Hochberg
    for (auto it = reg_tar_mi.begin(); it != reg_tar_mi.end(); ++it) {
      const auto k = it - reg_tar_mi.begin();
      const float p_k = getMIPVal(std::get<2>(*it), 100.0 / num_null_marginals);
      if (p_k < k * alpha / m)
        argmax_k = static_cast<uint32_t>(k) + 1;
    }
  } else if (method == "FWER") {
    for (auto it = reg_tar_mi.begin(); it != reg_tar_mi.end(); ++it) {
      const auto k = it - reg_tar_mi.begin();
      const float p_k = getMIPVal(std::get<2>(*it), 100.0 / num_null_marginals);
      if (p_k < alpha / m)
        argmax_k = static_cast<uint32_t>(k) + 1;
    }
  } else if (method == "FPR") {
    // Only for benchmarking
    for (auto it = reg_tar_mi.begin(); it != reg_tar_mi.end(); ++it) {
      const auto k = it - reg_tar_mi.begin();
      const float p_k = getMIPVal(std::get<2>(*it), 100.0 / num_null_marginals);
      if (p_k < alpha)
        argmax_k = static_cast<uint32_t>(k) + 1;
    }
  }

  // create the new vector that is a pruned version of original
  std::vector<std::tuple<gene_id, gene_id, float>> pruned_vec(
      &reg_tar_mi[0], &reg_tar_mi[argmax_k]);

  // rebuild network
  size_of_network = pruned_vec.size();
  gene_to_gene_to_float pruned_net;
  gene_to_gene_to_float pruned_net_reg_reg_only;

  pruned_net.reserve(regulators.size());
  pruned_net_reg_reg_only.reserve(regulators.size());

  for (const auto &[reg, tar, mi] : pruned_vec) {
    pruned_net[reg][tar] = mi;
    if (regulators.find(tar) != regulators.end())
      pruned_net_reg_reg_only[reg][tar] = mi;
  }

  return std::make_tuple(pruned_net, size_of_network, pruned_net_reg_reg_only);
}

/*
 Prune the network according to the MaxEnt weakest-edge reduction.
 */
gene_to_gene_to_float pruneMaxEnt(gene_to_gene_to_float network, uint32_t size_of_network, gene_to_gene_to_float &network_reg_reg_only) {

  gene_to_geneset edges_to_remove;
  edges_to_remove.reserve(regulators.size());

  // triangular matrix (reg and reg+1)
	for (const auto &[reg1, reg2_mi] : network_reg_reg_only) {
    std::unordered_map<gene_id, float> &reg1_regulon = network[reg1];
    std::unordered_set<gene_id> &remove_from_reg1 = edges_to_remove[reg1];
    for (const auto &[reg2, mi_regs] : reg2_mi) {
      network_reg_reg_only[reg2].erase(reg1);
      std::unordered_map<gene_id, float> &reg2_regulon = network[reg2];
      std::unordered_set<gene_id> &remove_from_reg2 = edges_to_remove[reg2];
      for (const auto &[tar, mi_reg1_tar] : reg1_regulon) {
        if (reg2_regulon.find(tar) != reg2_regulon.end()) {
          const float &mi_reg2_tar = reg2_regulon[tar];
          if (mi_reg1_tar < mi_regs && mi_reg1_tar < mi_reg2_tar)
            remove_from_reg1.insert(tar);
          else if (mi_reg2_tar < mi_regs && mi_reg2_tar < mi_reg1_tar)
            remove_from_reg2.insert(tar);
          else {
            remove_from_reg1.insert(reg2);
            remove_from_reg2.insert(reg1);
          }
        }
      }
    }
  }
	
	for (const auto &[reg, remove_from_reg] : edges_to_remove) {
    for (const auto &tar : remove_from_reg)
      network[reg].erase(tar);
		 size_of_network -= remove_from_reg.size();
  }
	
	return network;
}

/*
 Generates an ARACNe3 subnet (called from main).
*/
gene_to_gene_to_float ARACNe3_subnet(const gene_to_floats &subsample_exp_mat, const uint16_t subnet_num, const bool prune_alpha, const std::string& method, const float alpha, const bool prune_MaxEnt, const std::string& output_dir, const std::string& subnets_dir, const std::string& subnet_log_dir, const uint16_t nthreads) {
	std::ofstream log_output(subnet_log_dir + "log_subnet" + std::to_string(subnet_num) + ".txt");
	std::time_t t = std::time(nullptr);
	log_output << "---------" << std::put_time(std::localtime(&t), "%c %Z") << "---------" << std::endl << std::endl;
	log_output << "Subnetwork #: " + std::to_string(subnet_num) << std::endl;
	log_output << "Total # regulators (defined in gexp mat): " + std::to_string(regulators.size()) << std::endl;
	log_output << "Total # targets: " + std::to_string(genes.size()) << std::endl;
	log_output << "Total # samples: " + std::to_string(tot_num_samps) << std::endl;
	log_output << "Subsampled quantity: " + std::to_string(tot_num_subsample) << std::endl;
	log_output << "Total possible edges: " + std::to_string(regulators.size()*(genes.size()-1)) << std::endl;
	log_output << "Method of first pruning step: " + method << std::endl;
	log_output << "Alpha: " + std::to_string(alpha) << std::endl;
	log_output << "MaxEnt Pruning: " + std::to_string(prune_MaxEnt) << std::endl;
	log_output << "\n-----------Begin Network Generation-----------" << std::endl;
	
  // begin subnet computation

	//-------time module-------
  Watch watch1;
  log_output << "\nRaw network computation time: ";
  watch1.reset();
	//-------------------------
	
	uint32_t size_of_network = 0U;

  gene_to_gene_to_float network;
  network.reserve(regulators.size());
  for (uint16_t reg : regulators) {
    network[reg].reserve(regulators.size()*(genes.size()-1));
    for (uint16_t tar : genes)
      if (reg != tar) {
        network[reg][tar] = calcAPMI(subsample_exp_mat.at(reg), subsample_exp_mat.at(tar));
        size_of_network += 1;
      }
  }

	//-------time module-------
	log_output << watch1.getSeconds() << std::endl;
	log_output << "Size of network: " << size_of_network << " edges." << std::endl;
	//-------------------------
	
	//-------time module-------
  log_output << "\nThreshold pruning time (" + method + "): ";
  watch1.reset();
	//-------------------------
	
	uint32_t size_prev = size_of_network;
	
  // unpack tuple into objects
	gene_to_gene_to_float network_reg_reg_only;
  std::tie(network, size_of_network, network_reg_reg_only) = pruneAlpha(network, size_of_network, method, alpha);
	
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
		network = pruneMaxEnt(network, size_of_network, network_reg_reg_only);
		
		//-------time module-------
		log_output << watch1.getSeconds() << std::endl;
		log_output << "Edges removed: " << size_prev - size_of_network << " edges." << std::endl;
		log_output << "Size of network: " << size_of_network << " edges." << std::endl;
		//-------------------------
		
		uint32_t num_edges_after_MaxEnt_pruning = size_of_network;
		if (method == "FDR")
			FPR_estimates.emplace_back((alpha*num_edges_after_MaxEnt_pruning)/(regulators.size()*genes.size()-(1-alpha)*num_edges_after_threshold_pruning));
		else if (method == "FWER")
			FPR_estimates.emplace_back((alpha/(regulators.size()*(genes.size()-1)))*(num_edges_after_MaxEnt_pruning)/(num_edges_after_threshold_pruning));
		else if (method == "FPR")
			FPR_estimates.emplace_back(alpha*num_edges_after_MaxEnt_pruning/num_edges_after_threshold_pruning);
	} else {
		if (method == "FDR")
			FPR_estimates.emplace_back((alpha*num_edges_after_threshold_pruning)/(regulators.size()*genes.size()-(1-alpha)*num_edges_after_threshold_pruning));
		else if (method == "FWER")
			FPR_estimates.emplace_back(alpha/(regulators.size()*(genes.size()-1)));
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


const std::vector<consolidated_df_row> consolidate_subnets_vec(const std::vector<gene_to_gene_to_float> &subnets, const gene_to_floats &exp_mat, const gene_to_shorts &ranks_mat) {
	std::vector<consolidated_df_row> final_df;
	const uint32_t tot_poss_edgs = regulators.size()*(genes.size()-1);
	
	for (const gene_id &reg : regulators) {
		for (const gene_id &tar : genes) {
			uint16_t num_occurrences = 0;
			for (uint16_t sn = 0U; sn < subnets.size(); ++sn) {
				if (subnets[sn].at(reg).find(tar) != subnets[sn].at(reg).end())
					++num_occurrences;
			}
			if (num_occurrences > 0) {
				const float final_mi = calcAPMI(exp_mat.at(reg), exp_mat.at(tar));
				const float final_scc = consolidate_scc(ranks_mat.at(reg), ranks_mat.at(tar));
				const double final_p = right_tail_binomial_p(num_occurrences, subnets.size());
				final_df.emplace_back(reg, tar, final_mi, final_scc, num_occurrences, final_p);
			}
		}
	}
	
	return final_df;
}

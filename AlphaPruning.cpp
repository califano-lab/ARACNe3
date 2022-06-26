#include "ARACNe3.hpp"

extern uint16_t tot_num_regulators;
extern bool prune_MaxEnt;
extern float alpha;
extern std::string method;


/*
 Prunes a network by control of alpha using the Benjamini-Hochberg Procedure if method = FDR, or FWER if method = FWER..
 */
std::pair<reg_web, map_map> pruneAlpha(reg_web &network, uint32_t& size_of_network) {

	std::vector<std::pair<gene_id_t, edge_tar>> reg_edge_tar;
	reg_edge_tar.reserve(size_of_network);
	
	for (gene_id_t reg = 0; reg < tot_num_regulators; ++reg) {
		for (auto &et : network[reg]) {
			reg_edge_tar.emplace_back(reg, et);
		}
	}
	
	/*
	 Sort using the Comparator as a Lambda function.  std::sort sends < to sort in ascending order, so we use > for descending order.
	 */
	std::sort(reg_edge_tar.begin(), reg_edge_tar.end(), [](const std::pair<gene_id_t, edge_tar> &ret1, const std::pair<gene_id_t, edge_tar> &ret2) -> bool {return ret1.second.mi > ret2.second.mi;});
	
	uint32_t argmax_k = 0;
	uint32_t m = size_of_network;
	if (method == "FDR") {
		/*
		Apply the Benjamini-Hochberg principle; find argmax_k, the index to which we will prune.  We will include everything up to reg_tar_mi[argmax_k].
		*/
		for (auto it = reg_edge_tar.begin(); it != reg_edge_tar.end(); ++it) {
			auto k = it - reg_edge_tar.begin();
			float p_k = getMIPVal(it->second.mi);
			if (p_k <= k*alpha/m)
				argmax_k = static_cast<uint32_t>(k);
		}
	} else if (method == "FWER") {
	/*
	Control for FWER
	*/
		for (auto it = reg_edge_tar.begin(); it != reg_edge_tar.end(); ++it) {
			auto k = it - reg_edge_tar.begin();
			float p_k = getMIPVal(it->second.mi);
			if (p_k <= alpha/m)
				argmax_k = static_cast<uint32_t>(k);
		}
	}
	
	// create the new vector that is a pruned version of original
	std::vector<std::pair<gene_id_t, edge_tar>> pruned_vec(&reg_edge_tar[0], &reg_edge_tar[argmax_k+1]);
	
	// submit new network size to global variable defined in ARACNe3.cpp
	size_of_network = static_cast<uint32_t>(pruned_vec.size());
	
	/*
	 Rebuild the network.  We also fill out a data structure for the regulator-regulator interactions, so we can quickly identify regulators that regulate each other.  Rebuilding the network is inefficient due to the push_back.
	 */
	reg_web pruned_net;
	map_map tftfNetwork;
	pruned_net.reserve(tot_num_regulators);
	for (auto &pair : pruned_vec) {
		pruned_net[pair.first].push_back(pair.second);
		if (prune_MaxEnt) {
			/* Now, for this one, due to compression we know that if the identifier is below tot_num_regulators, we have a regulator!  Blazing fast! */
			if (pair.first < tot_num_regulators && pair.second.target < tot_num_regulators)
				tftfNetwork[pair.first][pair.second.target] = pair.second.mi;
		}
	}
	
	return std::make_pair(pruned_net, tftfNetwork);
}

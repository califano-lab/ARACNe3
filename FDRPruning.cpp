#include "ARACNe3.hpp"

extern uint32_t size_of_network;
extern uint16_t tot_num_regulators;
extern bool prune_DPI;

// adjacency matrix for DPI pruning
extern std::vector<std::vector<float>> reg_reg_edges;

/*
 Prunes a network by control of FDR using the Benjamini-Hochberg Procedure.
 */
reg_web pruneFDR(reg_web &network, uint32_t network_size, float FDR) {

	std::vector<std::pair<uint16_t, edge_tar>> reg_edge_tar;
	reg_edge_tar.reserve(network_size);
	
	for (uint16_t reg = 0; reg < tot_num_regulators; ++reg) {
		for (auto &et : network[reg]) {
			reg_edge_tar.emplace_back(reg, et);
		}
	}
	
	/*
	 Sort using the Comparator as a Lambda function.  std::sort sends < to sort in ascending order, so we use > for descending order.
	 */
	std::sort(reg_edge_tar.begin(), reg_edge_tar.end(), [](const std::pair<uint16_t, edge_tar> &ret1, const std::pair<uint16_t, edge_tar> &ret2) -> bool {return ret1.second.mi > ret2.second.mi;});
	
	/*
	 Apply the Benjamini-Hochberg principle; find argmax_k, the index to which we will prune.  We will include everything up to reg_tar_mi[argmax_k].
	 */
	uint32_t argmax_k = 0;
	uint32_t &m = network_size;
	for (auto it = reg_edge_tar.begin(); it != reg_edge_tar.end(); ++it) {
		auto k = it - reg_edge_tar.begin();
		float p_k = getMIPVal(it->second.mi);
		if (p_k <= k*FDR/m)
			argmax_k = static_cast<uint32_t>(k);
	}
	
	// create the new vector that is a pruned version of original
	std::vector<std::pair<uint16_t, edge_tar>> pruned_vec(&reg_edge_tar[0], &reg_edge_tar[argmax_k+1]);
	
	// submit new network size to global variable defined in ARACNe3.cpp
	size_of_network = static_cast<uint32_t>(pruned_vec.size());
	
	/*
	 Rebuild the network.  We also fill out an adjacency matrix for the regulator-regulator interactions, so we can quickly identify regulators that regulate each other.  Rebuilding the network is inefficient due to the push_back.
	 */
	if (prune_DPI)
		reg_reg_edges = std::vector<std::vector<float>>(tot_num_regulators, std::vector<float>(tot_num_regulators, 0.0f));
		
	reg_web pruned_net;
	pruned_net.reserve(tot_num_regulators);
	for (auto &pair : pruned_vec) {
		pruned_net[pair.first].push_back(pair.second);
		if (prune_DPI) {
			/* Now, for this one, due to compression we know that if the identifier is below tot_num_regulators, we have a regulator!  Blazing fast! */
			if (pair.first < tot_num_regulators && pair.second.target < tot_num_regulators)
				reg_reg_edges[pair.first][pair.second.target] = pair.second.mi;
		}
	}
	
	return pruned_net;
}

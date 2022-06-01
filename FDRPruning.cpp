#include "ARACNe3.hpp"

extern uint32_t size_of_network_unpruned;

/*
 Prunes a network by control of FDR using the Benjamini-Hochberg Procedure.
 */
reg_web pruneFDR(reg_web &network, std::vector<std::string> &regs, uint32_t network_size, float FDR) {

	std::vector<std::pair<std::string, edge_tar>> reg_edge_tar;
	reg_edge_tar.reserve(network_size);
	
	for (auto &reg : regs) {
		for (auto &et : network[reg]) {
			reg_edge_tar.emplace_back(reg, et);
		}
	}
	
	/*
	 Sort using the Comparator as a Lambda function.  std::sort sends < to sort in ascending order, so we use > for descending order.
	 */
	std::sort(reg_edge_tar.begin(), reg_edge_tar.end(), [](const std::pair<std::string, edge_tar> &ret1, const std::pair<std::string, edge_tar> &ret2) -> bool {return ret1.second.mi > ret2.second.mi;});
	
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
	std::vector<std::pair<std::string, edge_tar>> pruned_vec(&reg_edge_tar[0], &reg_edge_tar[argmax_k+1]);
	
	// submit new network size to global variable defined in ARACNe3.cpp
	size_of_network_unpruned = static_cast<uint32_t>(pruned_vec.size());
	
	/*
	 Rebuild the network.  This is where we should ideally also fill out an adjacency matrix for the regulators, so we can quickly identify regulators that regulate each other.  It's convenient but highly inefficient due to the push_back.
	 */
	reg_web pruned_net;
	pruned_net.reserve(regs.size());
	for (auto &pair : pruned_vec) {
		pruned_net[pair.first].push_back(pair.second);
	}
	
	return pruned_net;
}

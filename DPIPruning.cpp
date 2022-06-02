#include "ARACNe3.hpp"

// for DPI pruning, regulator -> regulator -> mi
std::unordered_map<reg_id_t, std::unordered_map<reg_id_t, float>> tftfNetwork;

extern uint32_t size_of_network;
extern uint16_t tot_num_regulators;

reg_web pruneDPI(reg_web &finalNet) {
	// primed will store the edges that are weakest.  we use a set to eliminate redundancy if the same edge is identified twice; same as hash set?
	std::unordered_map<reg_id_t, std::unordered_set<reg_id_t>> removedEdges;
	
	// must sort the network edge_tars based on target identifier (least->greatest) for below
	for (reg_id_t reg1 = 0; reg1 < tot_num_regulators; ++reg1) {
		if (tftfNetwork.contains(reg1)) {
			
			
		}
	}
	
/*
	// if regulator 1 -> regulator 2.  However then implicitly regulator 1->regulator 2. So this could be halved in terms of runtime
	const float r1_r2_mi = reg_reg_edges[reg1][reg2], r1_t_mi = r1_t[i].mi, r2_t_mi = r2_t[i].mi;
	if (r1_r2_mi >= r1_t_mi && r2_t_mi >= r1_t_mi) // if r1->t is weakest
		weakest.insert(std::make_pair(reg1, r1_t[i].target));
	else if (r1_r2_mi >= r2_t_mi && r1_t_mi >= r2_t_mi) // if r2->t is weakest
		weakest.insert(std::make_pair(reg2, r1_t[i].target));
	else // r1->r2 must be weakest
		weakest.insert(std::make_pair(reg1, reg2));*/
	
	/*
	 Now, we use 'weakest' to remove weakest edges.
	 
	 Consider using erase(remove_if()) idiom
	 */
	/*
	for (auto &p : weakest) {
		for (auto it = network[p.first].begin(); it != network[p.first].end(); ++it) {
			if (it->target == p.second)
			{
				network[p.first].erase(it--);
			}
		}
	}
	size_of_network -= weakest.size();
	 */
	return network;
}

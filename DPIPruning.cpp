#include "ARACNe3.hpp"

// adjacency matrix is a 2D array of float, for each regulator, of the r->r MI
std::vector<std::vector<float>> reg_reg_edges;

extern uint32_t size_of_network;
extern uint16_t tot_num_regulators;

reg_web pruneDPI(reg_web &network) {
	// primed will store the edges that are weakest.  we use a set to eliminate redundancy if the same edge is identified twice
	std::set<std::pair<uint16_t, uint16_t>> weakest;
	
	// if regulator 1 -> regulator 2.  However then implicitly regulator 1->regulator 2. So this could be halved in terms of runitme
	for (uint16_t reg1 = 0; reg1 < tot_num_regulators; ++reg1){
		for (uint16_t reg2 = 0; reg2 < tot_num_regulators; ++reg2) {
			if(reg_reg_edges[reg1][reg2]) {
			// iterators to the edge_tar structs must check if any target is equal.  O(N^2)
				for (auto it1 = network[reg1].begin(); it1 != network[reg1].end(); ++it1) {
					for (auto it2 = network[reg2].begin(); it2 != network[reg2].end(); ++it2) {
						if (it1->target == it2->target) {
							/*
							 If they share a target, we identify the weakest edge and prime it for deletion.  Remember that reg1 and reg2 are the gene identifiers for the regulators, and the it->target's are the gene identifiers for the targets.  Therefore, we extract mutual information, compare, and append the identifiers as a pair
							 */
							const float r1_r2_mi = reg_reg_edges[reg1][reg2], r1_t_mi = it1->mi, r2_t_mi = it2->mi;
							if (r1_r2_mi >= r1_t_mi && r2_t_mi >= r1_t_mi) // if r1->t is weakest
								weakest.insert(std::make_pair(reg1, it1->target));
							else if (r1_r2_mi >= r2_t_mi && r1_t_mi >= r2_t_mi) // if r2->t is weakest
								weakest.insert(std::make_pair(reg2, it2->target));
							else // r1->r2 must be weakest
								weakest.insert(std::make_pair(reg1, reg2));
						}
					}
				}
			}
		}
	}
	
	/*
	 Now, we use 'weakest' to remove weakest edges.
	 
	 Consider using erase(remove_if()) idiom
	 */
	for (auto &p : weakest) {
		for (auto it = network[p.first].begin(); it != network[p.first].end(); ++it) {
			if (it->target == p.second)
			{
				network[p.first].erase(it);
			}
		}
	}
	size_of_network -= weakest.size();
	return network;
}

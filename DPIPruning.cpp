#include "ARACNe3.hpp"

// adjacency matrix is a 2D array of float, for each regulator, of the r->r MI
std::vector<std::vector<float>> reg_reg_edges;

extern uint32_t size_of_network;
extern uint16_t tot_num_regulators;

reg_web pruneDPI(reg_web &network) {
	// primed will store the edges that are weakest.  we use a set to eliminate redundancy if the same edge is identified twice; same as hash set?
	std::set<std::pair<uint16_t, uint16_t>> weakest;
	
	// must sort the network edge_tars based on target identifier (least->greatest) for below
	std::cout << "SORTING EDGE TARS" << std::endl;
	sort_edge_tars(network);
	std::cout << "SORTING DONE" << std::endl;
	
	
	// if regulator 1 -> regulator 2.  However then implicitly regulator 1->regulator 2. So this could be halved in terms of runtime
	for (uint16_t reg1 = 0; reg1 < tot_num_regulators; ++reg1){
		std::cout << "reg1 " << reg1 << std::endl;
		for (uint16_t reg2 = 0; reg2 < tot_num_regulators; ++reg2) {
			std::cout << "reg2 " << reg2 << std::endl;
			if(reg_reg_edges[reg1][reg2]) {
				// identify shared targets
				std::vector<edge_tar> r1_t; //r1_edge_tars from intersection
				std::vector<edge_tar> r2_t; //r2_edge_tars from intersection
				
				/*
				 std::set_intersection REQUIRES that the vector is sorted and unique.  After the FDRPruning, each edge_tar vector is sorted largest -> smallest wrt MI, and does not contain duplicates, but we must sort by gene identifier, as this is what the intersection is on.
				 */
				std::cout << "TAKING INTERSECTION 1" << std::endl;
				std::set_intersection(network[reg1].begin(), network[reg1].end(), network[reg2].begin(), network[reg2].end(), std::back_inserter(r1_t));
				std::cout << "TAKING INTERSECTION 2" << std::endl;
				std::set_intersection(network[reg2].begin(), network[reg2].end(), network[reg1].begin(), network[reg1].end(), std::back_inserter(r2_t));
				std::cout << "INTERSECTIONS DONE" << std::endl;
				
				std::cout << "PUSHING WEAKEST EDGES" << std::endl;
				for (uint32_t i = 0; i < r1_t.size(); ++i) {
						/*
						 If they share a target, we identify the weakest edge and prime it for deletion.  Remember that reg1 and reg2 are the gene identifiers for the regulators, and the it->target's are the gene identifiers for the targets.  Therefore, we extract mutual information, compare, and append the identifiers as a pair
						 */
						const float r1_r2_mi = reg_reg_edges[reg1][reg2], r1_t_mi = r1_t[i].mi, r2_t_mi = r2_t[i].mi;
						if (r1_r2_mi >= r1_t_mi && r2_t_mi >= r1_t_mi) // if r1->t is weakest
							weakest.insert(std::make_pair(reg1, r1_t[i].target));
						else if (r1_r2_mi >= r2_t_mi && r1_t_mi >= r2_t_mi) // if r2->t is weakest
							weakest.insert(std::make_pair(reg2, r1_t[i].target));
						else // r1->r2 must be weakest
							weakest.insert(std::make_pair(reg1, reg2));
				}
				std::cout << "PUSHING WEAKEST EDGES DONE" << std::endl;
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
				network[p.first].erase(it--);
			}
		}
	}
	size_of_network -= weakest.size();
	return network;
}

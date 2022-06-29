#include "ARACNe3.hpp"
#include "omp.h"

// for MaxEnt pruning, regulator -> regulator -> mi
extern uint16_t tot_num_regulators;
extern uint16_t nthreads;


/*
 Prune the network according to the MaxEnt weakest-edge reduction.
 */
reg_web pruneMaxEnt(reg_web& network, map_map& tftfNetwork, uint32_t &size_of_network) {
	// primed will store the edges that are weakest.  we use a set to eliminate redundancy if the same edge is identified twice; same as hash set?
	std::vector<std::vector<std::set<gene_id_t>>> removedEdgesForThread(nthreads, std::vector<std::set<gene_id_t>>(tot_num_regulators));
	
	/*
	 Inefficient conversion operation here.  Makes searching whether a target is contained an O(1) operation due to the hash map of hash maps, as opposed to a hash map of edge_tar vectors, which would make checking for a particular edge_tar.target an O(n) operation.
	 */
			
	map_map finalNet = regweb_to_mapmap(network); 
	
	// must sort the network edge_tars based on target identifier (least->greatest) for below
#pragma omp parallel for
	for (gene_id_t reg1 = 0; reg1 < tot_num_regulators; ++reg1) {
		if (tftfNetwork.contains(reg1)) {
			auto &fin1 = finalNet[reg1];
			auto &tft1 = tftfNetwork[reg1]; 
			auto &rem1 = removedEdgesForThread[omp_get_thread_num()][reg1];
			for (gene_id_t reg2 = reg1 + 1; reg2 < tot_num_regulators; ++reg2) {
				if (tft1.contains(reg2)) {
					auto &fin2 = finalNet[reg2];
					auto &rem2 = removedEdgesForThread[omp_get_thread_num()][reg2];
					const float tftfMI = tft1[reg2];
					for(const auto &[target, v2] : fin2) {
						if (fin1.contains(target)) {
							const float v1 = fin1[target];
							if (v1 < tftfMI && v1 < v2)
								rem1.insert(target);
							else if (v2 < tftfMI && v2 < v1)
								rem2.insert(target);
							else {
								// is ARACNe-AP removing the r1-r2 edge?
								rem1.insert(reg2);
								rem2.insert(reg1);
							}
						}
					}
				}
			}
		}
	}
	
	std::vector<std::set<gene_id_t>> removedEdges(tot_num_regulators);
	for (gene_id_t reg = 0; reg < tot_num_regulators; ++reg) {
		for (uint16_t th = 0; th < nthreads; ++th) {
			for (const auto &tar : removedEdgesForThread[th][reg])
				removedEdges[reg].insert(tar);
		}
	}
	
	for (const auto &removedSet : removedEdges) 
		 size_of_network -= removedSet.size();
	
	reg_web pruned_net;
	pruned_net.reserve(tot_num_regulators);
	for (const auto &[reg, tarmap] : finalNet) {
		pruned_net[reg].reserve(network[reg].size());
		auto &rem = removedEdges[reg];
		for (const auto &[tar, mi] : tarmap) {
			if (!rem.contains(tar)) {
				pruned_net[reg].emplace_back(tar, mi);
			}
		}
	}

	return pruned_net;
}

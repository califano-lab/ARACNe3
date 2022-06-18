/*
 Consolidator for ARACNe3.  Contains various functions only needed in the consolidation step, such as calculation of the p-value for an edge based on the number of subnetworks it appeared in, the calculation of the SCC for an edge, etc.
 */

#include "ARACNe3.hpp"

extern uint16_t tot_num_samps_pre_subsample;
extern uint16_t tot_num_samps;
extern uint16_t tot_num_regulators;
extern genemap global_gm;
extern uint16_t num_subnets;

float consolidate_scc(const std::vector<float>& vec_x, const std::vector<float>& vec_y) {
	/*
	 A ranking is formed in the following way.  Indices index = [0,subsample_quant) are sorted based on the ranking of expr_vec_sampled[index], so that we get some new sorted set of indexes (5, 2, 9, ... ) that is the rank of each element in expr_vec_sampled
	 
	 For a lambda function, brackets indicate the scope of the function.
	 */
	std::vector<uint16_t> x_ranked = rank_vals(vec_x), y_ranked = rank_vals(vec_y);
	return 0.0f;
}

float consolidate_p(const uint16_t& num_occurrences) {
	return 0.0f;
}

std::vector<consolidated_df> consolidate(std::vector<reg_web> &subnets) {
	std::vector<consolidated_df> final_df;
	const auto tot_poss_edgs = tot_num_regulators*global_gm.size()-tot_num_regulators;
	final_df.reserve(tot_poss_edgs);
	
	std::vector<map_map> subnets_mpmp;
	for (uint16_t i = 0; i < subnets.size(); ++i) {
		subnets_mpmp.emplace_back(regweb_to_mapmap(subnets[i]));
	}
	
	
		for (uint16_t reg = 0; reg < tot_num_regulators; ++reg) {
			for (uint16_t tar = 0; tar < global_gm.size(); ++tar) {
				uint16_t num_occurrences = 0;
				for (uint16_t sn = 0; sn < subnets.size(); ++sn) {
					if (subnets_mpmp[sn][reg].contains(tar))
						++num_occurrences;
				}
				if (num_occurrences > 0) {
					const float final_mi = APMI(global_gm[reg], global_gm[tar]);
					const float final_scc = consolidate_scc(global_gm[reg], global_gm[tar]);
					const float final_p = consolidate_p(num_occurrences);
					
					final_df.emplace_back(reg, tar, final_mi, final_scc, num_occurrences, final_p);
				}
		}
	}
	
	return final_df;
}

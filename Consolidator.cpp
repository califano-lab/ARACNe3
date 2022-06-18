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
	std::vector<uint16_t> x_ranked = rank_vals(vec_x), y_ranked = rank_vals(vec_y);
	float sigma = 0, sigmaxy = 0, sigmasq = 0;
	for (uint16_t i = 0; i < vec_x.size(); ++i) {
		sigma += x_ranked[i]; // same for x and y
		sigmaxy += x_ranked[i]*y_ranked[i];
		sigmasq += x_ranked[i]*x_ranked[i]; // same for x and y
	}
	return (vec_x.size() * sigmaxy - sigma*sigma)/
			(float) (vec_x.size() * sigmasq - sigma*sigma);
}


/*
 We can do factorial on at most an 8-bit integer.
 */
uint32_t factorial(const uint8_t& n) {
	return n == 1U ? n : n * factorial(n - 1);
}

uint32_t n_choose_r(uint8_t n, uint8_t r) {
	return factorial(n) / (factorial(n - r) * factorial(r));
}

double right_tail_binomial_p(const uint16_t& num_occurrences) {
	float theta = 1.5E-4f;
	double p = 0.0f;
	for (uint16_t i = num_subnets; i >= num_occurrences; --i) 
		p += n_choose_r(num_subnets, i) * std::pow(theta, i) * std::pow(1-theta,num_subnets-i);
	return p;
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
					const double final_p = right_tail_binomial_p(num_occurrences);
					
					final_df.emplace_back(reg, tar, final_mi, final_scc, num_occurrences, final_p);
				}
		}
	}
	
	return final_df;
}

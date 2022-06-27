/*
 Consolidator for ARACNe3.  Contains various functions only needed in the consolidation step, such as calculation of the p-value for an edge based on the number of subnetworks it appeared in, the calculation of the SCC for an edge, etc.
 */

#include "ARACNe3.hpp"

float FPR_estimate = 1.5E-4f;

extern uint16_t tot_num_samps;
extern uint16_t tot_num_subsample;
extern uint16_t tot_num_regulators;
extern genemap global_gm;
extern genemap_r global_gm_r; 
extern uint16_t num_subnets;

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

double right_tail_binomial_p(const uint16_t& num_occurrences) {
	double p = 0.0;
	if (num_subnets == 1)
		return std::numeric_limits<double>::quiet_NaN(); // cannot have a p-value for 1 subnet (1 network)
	for (uint16_t i = num_subnets; i >= num_occurrences; --i)
		p += std::exp(lchoose(num_subnets, i) + i * std::log(FPR_estimate) + (num_subnets - i) * std::log(1-FPR_estimate));
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
					const float final_scc = consolidate_scc(global_gm_r[reg], global_gm_r[tar]);
					const double final_p = right_tail_binomial_p(num_occurrences);
					final_df.emplace_back(reg, tar, final_mi, final_scc, num_occurrences, final_p);
				}
		}
	}
	
	return final_df;
}

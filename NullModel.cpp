#include "ARACNe3.hpp"

/*
 * This file is the null-model module of ARACNe3.  It is a separate file because
 * for the first pruning step, we will need to calculate a p-value based on the
 * eCDF for the null mutual information of gene expression marginals, calculated
 * via APMI from 1,000,000 random vectors.  The vector is sorted and stored on
 * the heap, and it will be accessible for p-value calculation from the global
 * variable pointer null_mis.  null_mis is to be free'd after the first pruning
 */

static const uint32_t num_null_marginals = 1000000;
static std::vector<float> null_mis;

extern bool verbose;
extern std::string cached_dir;

/*
 Will cache the null_mis vector in the cached directory.  The vector will be a BLOB named 'Null_tot_num_samps_###'.
 */
void cacheNullModel(uint16_t tot_num_samps, const std::vector<float> &mi_vec) {
	std::string filename = cached_dir + "Null_tot_num_samps_" + std::to_string(tot_num_samps);
	std::ofstream cached(filename, std::ios::out | std::ios::binary);

	//tab-delinate the data
	std::ostream_iterator<float> cached_iterator(cached, "\t");
	std::copy(mi_vec.begin(), mi_vec.end(), cached_iterator);
	return;
}

/*
 Computes 1 million null mutual information values for the sample size.  Checks whether there already exists a null_mi vector (filename) in the cached directory.
 */
const std::vector<float> initNullMIs(uint16_t tot_num_samps) {
	std::string filename = cached_dir + "Null_tot_num_samps_" + std::to_string(tot_num_samps);
	if (std::filesystem::exists(filename)) {
		std::ifstream cached(filename, std::ios::in | std::ios::binary);
		
		//read tab-delineated data
		std::vector<float> mi_vec;
		mi_vec.reserve(num_null_marginals);
		std::istream_iterator<float> cached_iterator(cached);
		while(cached_iterator != std::istream_iterator<float>()) /*default is off-end iterator */
			mi_vec.emplace_back(*cached_iterator++);
		null_mis = mi_vec;
		return mi_vec;
	} else {
		// make the permute vector, the ref vector, send to permuteAPMI
		std::vector<float> ref_vec;
		ref_vec.reserve(tot_num_samps);
		for (uint16_t i = 1; i <= tot_num_samps; ++i) {
			ref_vec.push_back(((float) i)/(tot_num_samps+1));
		}

		// vector of vectors, 1mil rows
		std::vector<std::vector<float>> target_vec(num_null_marginals, ref_vec);

		auto rng = std::default_random_engine {};
		for (unsigned int i = 0; i < num_null_marginals; ++i) {
			std::shuffle(std::begin(target_vec[i]), std::end(target_vec[i]),
					rng);
		}

		// get the 1 million MI values
		std::vector<float> mi_vec = permuteAPMI(ref_vec, target_vec, 7.815, 4);
		
		// standard sorting, smallest to largest
		std::sort(mi_vec.begin(), mi_vec.end(), std::less<float>());
		
		// cache the null model in memory for future calculations
		cacheNullModel(tot_num_samps, mi_vec);
		
		
		// finally, set file static variable and/or return mi_vec to caller
		null_mis = mi_vec;
		return mi_vec;
	}
	

}

const float getMIPVal(const float &mi) {
	// returns an iterator to the highest index mi is less than
	auto it = std::upper_bound(null_mis.begin(), null_mis.end(), mi);
	
	// p-value is 1 minus the percentile
	return (1 - (it-null_mis.begin())/(float)num_null_marginals);
}

const std::vector<float> getMIPVals(const std::vector<float> &mis) {
	std::vector<float> ps;
	ps.reserve(mis.size());
	for (auto &mi : mis)
		ps.emplace_back(getMIPVal(mi));
	return ps;
}

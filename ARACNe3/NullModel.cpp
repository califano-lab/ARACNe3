#include "ARACNe3.hpp"

/*
 * This file is the null-model module of ARACNe3.  It is a separate file because
 * for the first pruning step, we will need to calculate a p-value based on the
 * eCDF for the null mutual information of gene expression marginals, calculated
 * via APMI from 1,000,000 random vectors.  The vector is sorted and stored on
 * the heap, and it will be accessible for p-value calculation from the global
 * variable pointer null_mis.  null_mis is to be free'd after the first pruning
 * step.
 *
 */

static std::vector<float> *null_mis;

/*
 * Computes 1 million null mutual information values from the number of samples.
 * To reduce runtime, we emulate rowAPMI but 
 */
const std::vector<const float> initNullMIs(int tot_num_samps) {
	// make the permute vector, the ref vector, send to permuteAPMI
	std::vector<float> ref_vec;
	ref_vec.reserve(tot_num_samps);
	for (uint16_t i = 1; i <= tot_num_samps; ++i) {
		ref_vec.push_back(((float) i)/(tot_num_samps+1));
	}

	// an array of vectors, pointer on stack; array on heap
	std::vector<std::vector<float>> target_vec;
	target_vec.reserve(1000000);

	auto rng = std::default_random_engine {};
	for (unsigned int i = 0; i < 1000000; ++i) {
		target_vec.emplace_back(std::vector<float>(ref_vec));
		std::shuffle(std::begin(target_vec[i]), std::end(target_vec[i]),
				rng);
	}

	// vector is now on heap
	const std::vector<const float> mi_vec = permuteAPMI(ref_vec, target_vec,
			7.815, 4);

	// NEXT -> Sort vector largest to smallest and output.  

	return mi_vec;
}

int main() {
	const std::vector<const float> mis = initNullMIs(100);
	//for (auto &num : mis) { std::cerr << num << std::endl; }
	return 0;
}

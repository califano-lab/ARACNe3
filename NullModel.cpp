#include "ARACNe3.hpp"

/*
 This file contains the null-model module of ARACNe3.  We store the 1,000,000 null MI values in a file static variable vector (which is really a file static pointer to the vector on the heap).  Storing this as a global variable makes it efficient to get p-values because we don't need to pass a pointer each time we calculate p, but we can refer to the file static variable.  We could get around using file static variables by making this entire file into more of a class system, but what we have here is fine.
 */
uint32_t num_null_marginals = 1000000; // number of marginals (1 mil.)

static std::vector<float> null_mis; // null MI vec
static float m, b; // parameters for OLS regression p-val estimation

/*
 Global variables are passed from ARACNe3.cpp, which are the user-defined parameters.  
 */
extern bool verbose;
extern std::string cached_dir;
extern uint32_t global_seed;

/*
 Will cache the null_mis vector in the cached directory.  The vector will be a BLOB named 'Null_tot_num_samps_###'.
 */
void cacheNullModel(const std::vector<float>& mi_vec, const float& m, const double& b, const std::string& filename) {
	std::ofstream cached(filename, std::ios::out | std::ios::binary);

	//comma-delinate the data
	std::ostream_iterator<float> cached_iterator(cached, " ");
	std::copy(mi_vec.begin(), mi_vec.end(), cached_iterator);
	
	cached << "\n" << m << " " << b;
	
	return;
}


inline static float sqr(float x) { return x*x; }

std::pair<float, float> linreg(uint32_t N, const std::vector<float>& x, const std::vector<float>& y){
	
	float sumx = 0.0f, sumx2 = 0.0f, sumxy = 0.0f, sumy = 0.0f, sumy2 = 0.0f;
	
	for (uint32_t i = 0; i < N; ++i) { 
		sumx  += x[i];       
		sumx2 += sqr(x[i]);  
		sumxy += x[i] * y[i];
		sumy  += y[i];      
		sumy2 += sqr(y[i]); 
	}

	float denom = (N * sumx2 - sqr(sumx));
	float m = 0.0f, b = 0.0f;
	if (denom == 0) {
		std::cout << "SINGULAR MATRIX. COULD NOT FIT PIECEWISE NULL MODEL." << std::endl;
		return std::make_pair(m, b);
	}
	
	m = (N * sumxy  -  sumx * sumy) / denom;
	b = (sumy * sumx2  -  sumx * sumxy) / denom;
	
	return std::make_pair(m, b); 
}


/*
 Computes 1 million null mutual information values for the sample size.  Checks whether there already exists a null_mi vector (nulls_filename) in the cached directory.
 */
const std::vector<float> initNullMIs(const uint16_t& tot_num_samps) {
	std::string nulls_filename = cached_dir + "numnulls-" + std::to_string(num_null_marginals) + "_numsamps-" + std::to_string(tot_num_samps);
	/*
	 If there already is a null model for this number of samples cached in the cached_dir, then we just pull values from that.  Also pull parameters from regression.
	 */
	if (std::filesystem::exists(nulls_filename)) {
		std::ifstream cached(nulls_filename, std::ios::in | std::ios::binary);
		
		//read tab-delineated data
		std::vector<float> mi_vec;
		mi_vec.reserve(num_null_marginals);
		std::istream_iterator<float> cached_iterator(cached);
		uint32_t i = 0;
		for (; i < num_null_marginals; ++i) {
			mi_vec.emplace_back(*cached_iterator++);
		}
		m = *cached_iterator++;
		b = *cached_iterator;
		
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

		auto rng = std::default_random_engine{global_seed};
		for (unsigned int i = 0; i < num_null_marginals; ++i) {
			std::shuffle(std::begin(target_vec[i]), std::end(target_vec[i]),
					rng);
		}

		// get the 1 million MI values
		std::vector<float> mi_vec = permuteAPMI(ref_vec, target_vec, 7.815, 4);
		
		// standard sorting, largest to smallest
		std::sort(mi_vec.begin(), mi_vec.end(), std::greater<float>());
		
		/* Fit an ordinary least squares regression log(p) vs. MI values for MI values with with eCDF p < 0.01
		 */
		uint32_t p_99th_percentile_idx = std::ceil(num_null_marginals/100.0);
		std::vector<float> indices(p_99th_percentile_idx), log_ps(p_99th_percentile_idx);
		
		std::iota(indices.begin(), indices.end(), 0.0f);
		std::transform(mi_vec.begin(), mi_vec.begin()+p_99th_percentile_idx, log_ps.begin(), [](const auto &mi) -> float { return std::log(mi); });
		
		
		std::pair<float, float> sol = linreg(p_99th_percentile_idx, indices, log_ps);
		m = sol.first;
		b = sol.second;
		
		// cache the null model in memory for future calculations
		cacheNullModel(mi_vec, m, b, nulls_filename);
		
		// finally, set file static variable and/or return mi_vec to caller
		null_mis = mi_vec;
		return mi_vec;
	}
	

}

/*
 This function uses the sorted null distribution vector to calculate the p-value of a MI value.  If the associated p-value is below p_precise, the p-value will be calculated from an OLS regression of log-p values from the eCDF (p < 0.01)
 */
const float getMIPVal(const float& mi, const float& p_precise) {
	// . 
	/* Returns an iterator that points to the first index for which mi is greater than the rest.  
	 */
	auto it = std::upper_bound(null_mis.begin(), null_mis.end(), mi, std::greater<float>());
	
	/* p-value is the percentile.  Since it's an index, it underestimates p-value, so we add 1 later.
	 */
	const float p = (it-null_mis.begin()+1)/(float)num_null_marginals;
	
	if (p < p_precise)
		return std::exp(m*std::logf(p)+b); // invert log
	else
		return p;
}

/*
 This is a wrapper for getMIPVal which vectorizes the getMIPVal calculation.  This will mainly be used for Rcpp.
 */
const std::vector<float> getMIPVals(const std::vector<float>& mis, const float& p_precise) {
	std::vector<float> ps;
	ps.reserve(mis.size());
	for (auto &mi : mis)
		ps.emplace_back(getMIPVal(mi, p_precise));
	return ps;
}

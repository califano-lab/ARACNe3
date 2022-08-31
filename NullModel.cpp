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
extern std::string output_dir;
extern uint32_t global_seed;
extern uint16_t nthreads;

/*
 Will cache the null_mis vector as the specified filename.
 */
void cacheNullModel(const std::vector<float>& mi_vec, const float& m, const double& b, const std::string& nulls_filename) {
	std::string OLS_coef_filename = nulls_filename + "_OLS";
	
	std::ofstream nulls_file(nulls_filename + ".txt", std::ios::out | std::ios::binary);
	std::ofstream OLS_coef_file(OLS_coef_filename + ".txt", std::ios::out | std::ios::binary);
	
	for (auto it = mi_vec.cbegin(); it != mi_vec.cend(); ++it)
		nulls_file << *it << '\n';
	
	OLS_coef_file << m << '\n' << b << '\n';
	
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
const std::vector<float> initNullMIs(const uint16_t& tot_num_subsample) {
	std::string nulls_filename = "numnulls-" + std::to_string(num_null_marginals) + "_numsamps-" + std::to_string(tot_num_subsample);
	std::string OLS_coef_filename = nulls_filename + "_OLS";
	/*
	 If there already is a null model for this number of samples cached in the cached_dir, then we just pull values from that.  Also pull parameters from regression.
	 */
	if ( /*never use cache for debugging;*/ std::filesystem::exists(cached_dir + nulls_filename + ".txt") && std::filesystem::exists(cached_dir + OLS_coef_filename + ".txt")) {

		std::ifstream nulls_file(cached_dir + nulls_filename + ".txt", std::ios::in | std::ios::binary);
		std::ifstream OLS_coef_file(cached_dir + OLS_coef_filename + ".txt", std::ios::in | std::ios::binary);
		
		std::vector<float> mi_vec;
		mi_vec.reserve(num_null_marginals);
		std::istream_iterator<float> nulls_iterator(nulls_file);
		std::istream_iterator<float> OLS_iterator(OLS_coef_file);
		uint32_t i = 0;
		for (; i < num_null_marginals; ++i)
			mi_vec.emplace_back(*nulls_iterator++);
		m = *OLS_iterator++;
		b = *OLS_iterator;
		
		global_seed++; //global_seed is incremented when calculating null
		null_mis = mi_vec;
	} else {
		// make the permute vector, the ref vector, send to permuteAPMI
		std::vector<float> ref_vec;
		ref_vec.reserve(tot_num_subsample);
		for (uint16_t i = 1; i <= tot_num_subsample; ++i)
			ref_vec.emplace_back(((float) i)/(tot_num_subsample+1));

		// vector of vectors, 1mil rows
		std::vector<std::vector<float>> target_vec(num_null_marginals, ref_vec);

		static std::mt19937 rand{global_seed++};
		
		// Cannot be parallelized because access to random generator
		for (unsigned int i = 0; i < num_null_marginals; ++i)
			std::shuffle(target_vec[i].begin(), target_vec[i].end(), rand);

		// get the 1 million MI values
		std::vector<float> mi_vec = permuteAPMI(ref_vec, target_vec, 7.815, 4);
		
		// standard sorting, largest to smallest
		std::sort(mi_vec.begin(), mi_vec.end(), std::greater<float>());
		
		/* Fit an ordinary least squares regression log(p) vs. MI values for MI values with with eCDF p < 0.01
		 */
		uint32_t p_99th_percentile_idx = std::ceil(num_null_marginals/100.0);

		// we will be doing a linear regression of the log(mi_ps)
		std::vector<float> mis(mi_vec.begin(), mi_vec.begin() + p_99th_percentile_idx), mi_ps(p_99th_percentile_idx);

#pragma omp parallel for num_threads(nthreads)
		for (int i = 0; i < p_99th_percentile_idx; ++i)
			mi_ps[i] = ((i+1)/(float)num_null_marginals);
		
		// mi_ps turned into log(mi_ps)
		std::transform(mi_ps.begin(), mi_ps.begin()+p_99th_percentile_idx, mi_ps.begin(), [](const auto &p) -> float { return std::log(p); });
		
		
		std::pair<float, float> sol = linreg(p_99th_percentile_idx, mis, mi_ps);
		m = sol.first;
		b = sol.second;
		
		// cache the null model in memory for future calculations
		cacheNullModel(mi_vec, m, b, cached_dir + nulls_filename);
		
		// finally, set file static variable and/or return mi_vec to caller
		null_mis = mi_vec;
	}
	cacheNullModel(null_mis, m, b, output_dir + nulls_filename);
	return null_mis;
}

/*
 This function uses the sorted null distribution vector to calculate the p-value of a MI value.  If the associated p-value is below p_precise, the p-value will be calculated from an OLS regression of log-p values from the eCDF (p < 0.01)
 */
const float getMIPVal(const float& mi, const float& p_precise) {
	/* Returns an iterator that points to the first index for which mi > the rest.  
	 */
	auto it = std::upper_bound(null_mis.begin(), null_mis.end(), mi, std::greater<float>());
	
	/* p-value is the percentile.  Since it's an index, it underestimates p-value, so we add 1 later.
	 */
	const float p = (it-null_mis.begin()+1)/(float)num_null_marginals;
	
	if (p < p_precise)
		return std::min(p, std::exp(m*mi+b)); // invert log
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

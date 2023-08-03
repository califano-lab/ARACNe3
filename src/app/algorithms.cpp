#include "MINDy3.hpp"
#include "algorithms.hpp"
#include <boost/math/distributions/chi_squared.hpp>

static float q_thresh = 7.815f;
static uint16_t size_thresh = 4U;

/** @brief Ranks indices based on the values in vec. 
 * 
 * This function sorts the indices in the range [1, size) based on the values
 * of vec[index-1]. The returned vector represents the rank of each element in
 * vec with the smallest element ranked as 1 and the largest as size. If two
 * elements in vec have the same value, their corresponding indices in the
 * ranking are randomly shuffled.
 *
 * @param vec The input vector for which the ranking should be formed.
 * @param rand A Mersenne Twister pseudo-random generator of 32-bit numbers
 * with a state size of 19937 bits. Used to shuffle indices corresponding to
 * equal values in vec.
 * 
 * @return A vector representing the rank of indices in the input vector vec.
 *
 * @example vec = {9.2, 3.5, 7.4, 3.5} The function returns {4, 1, 3, 2}. Note
 * that the ranks for the elements with the same value 3.5 (indices 1 and 3)
 * may be shuffled differently in different runs.
 */
std::vector<uint16_t> rankIndices(const std::vector<float>& vec, std::mt19937 &rand) {
	std::vector<uint16_t> idx_ranks(vec.size());
	std::iota(idx_ranks.begin(), idx_ranks.end(), 0U); /* 0, 1, ..., size-1 */
	std::sort(idx_ranks.begin(), idx_ranks.end(), [&vec](const uint16_t &num1, const uint16_t &num2) -> bool { return vec[num1] < vec[num2]; }); /* sort ascending */
	for (uint16_t r = 0U; r < idx_ranks.size();) {
		uint16_t same_range = 1U;
		while (r + same_range < idx_ranks.size() && vec[idx_ranks[r]] == vec[idx_ranks[r+same_range]])
			++same_range; // same_range is off-end index
		if (same_range > 1U) {
			std::shuffle(idx_ranks.begin()+r, idx_ranks.begin()+r+same_range, rand);
			r = r + same_range;
		} else {
			++r;
		}
	}
	return idx_ranks;
}

/**
 * @brief Calculates the Spearman's Rank Correlation Coefficient (SCC) between
 * two ranked vectors.
 * 
 * This function computes the Spearman's Rank Correlation Coefficient, a
 * statistical measure of the monotonic relationship between two ranked
 * variables. Both x_ranked and y_ranked are expected to have equal size and to
 * be ranked in ascending order.
 *
 * @param x_ranked Ranks of variable X.
 * @param y_ranked Ranks of variable Y.
 * 
 * @return The SCC. This value lies between -1 and 1 inclusive, where 1
 * signifies a perfect positive correlation, -1 a perfect negative correlation,
 * and 0 indicates no correlation.
 */
float spearmanCorrelate(const std::vector<uint16_t>& x_ranked, const std::vector<uint16_t>& y_ranked) {
	const auto &n = x_ranked.size();
	int sigma_dxy = 0;
	for (uint16_t i = 0; i < n; ++i)
		sigma_dxy += (x_ranked[i] - y_ranked[i]) * (x_ranked[i] - y_ranked[i]);
	return 1 - 6.0f * sigma_dxy / n / (n * n - 1);
}

/**
 * @brief Combines multiple p-values into one, using Fisher's method.
 * 
 * This function uses Fisher's method for combining independent p-values from
 * several statistical tests bearing upon the same overall hypothesis. This
 * method transforms the p-values into a chi-squared distribution.
 * 
 * @param p_vals A vector of p-values obtained from independent tests.
 * 
 * @return A single p-value obtained by combining the input p-values using
 * Fisher's method. 
 */
float fishersMethodP(const std::vector<float>& p_vals) {
  float log_sum = -2*std::accumulate(p_vals.begin(), p_vals.end(), 0.0f,
      [](float acc, float value) { return acc + std::log(value); });
  uint32_t degrees_of_freedom = 2 * p_vals.size();
  boost::math::chi_squared_distribution dist(degrees_of_freedom);
  return float(1.0 - boost::math::cdf(dist, log_sum));
}

/**
 * @brief Calculate linear regression to find the slope and y-intercept.
 * 
 * This function takes two float vector parameters, x and y, representing data
 * points on a two-dimensional plane. It returns a pair of floats where the
 * first float is the slope (m) and the second float is the y-intercept (b)
 * from the line equation y = mx + b.
 * 
 * @param x A vector of x floats.
 * @param y A vector of y floats.
 *
 * @return A pair of floats (m,b).
 */
std::pair<float, float> linearRegress(const std::vector<float>& x, const std::vector<float>& y){
	std::vector<float>::size_type N = x.size();
	
	float sumx = 0.0f, sumx2 = 0.0f, sumxy = 0.0f, sumy = 0.0f, sumy2 = 0.0f;
	
	for (uint32_t i = 0; i < N; ++i) { 
		sumx  += x[i];       
		sumx2 += x[i]*x[i];  
		sumxy += x[i] * y[i];
		sumy  += y[i];      
		sumy2 += y[i]*y[i]; 
	}

	float denom = (N * sumx2 - sumx*sumx);
	float m = 0.0f, b = 0.0f;
	if (denom == 0) {
		std::cout << "SINGULAR MATRIX. COULD NOT FIT PIECEWISE NULL MODEL." << std::endl;
		return std::make_pair(m, b);
	}
	
	m = (N * sumxy  -  sumx * sumy) / denom;
	b = (sumy * sumx2  -  sumx * sumxy) / denom;
	
	return std::make_pair(m, b); 
}

/**
 * @brief Calculate the Mutual Information (MI) for a square struct.
 *
 * @param s The square structure for which MI is calculated.
 * 
 * @return A float representing the MI of the input square struct.
 */
float calcMI(const square &s) {
	const float pxy = s.num_pts/(float)s.tot_num_pts, marginal = s.width, mi = pxy*std::log(pxy/(marginal*marginal));
	return std::isfinite(mi) ? mi : 0.0f;
}

/**
 * @brief Perform recursive tessellation of XY plane and MI calculation at
 * dead-ends.
 *
 * This function recursively tessellates the XY plane and calculates the Mutual
 * Information (MI) at dead-ends. The function divides the plane into four
 * quadrants and computes the chi-square statistic for the distribution of
 * points across the quadrants. If the chi-square statistic exceeds a
 * predefined threshold or if the function is operating on the initial square,
 * the function continues to subdivide the plane. Otherwise, it calculates the
 * MI for the current square.
 *
 * Note: The use of alloca is not safe from stack overflow, strictly speaking,
 * but stack overflow is only conceivable in a contrived case, and any generic
 * usage of MINDy3 should never create a stack allocation too large.
 *
 * @param x_ptr Pointer to the x-coordinate data.
 * @param y_ptr Pointer to the y-coordinate data.
 * @param s The square struct on which to perform a tessellation.
 *
 * @return A float value representing the result of MI calculations, or
 * performs 
 * recursion, depending on the condition.
 */
const float calcAPMISplit(const float *const x_ptr, const float *const y_ptr, const square &s) {
	const float &x_bound1=s.x_bound1, &y_bound1=s.y_bound1, &width=s.width;

	// if we have less points in the square than size_thresh, calc MI
	if (s.num_pts < size_thresh) { return calcMI(s);}

	// thresholds for potential partition of XY plane
	const float x_thresh = x_bound1 + width*0.5f,
		  y_thresh = y_bound1 + width*0.5f;

	// indices for quadrants, to test chi-square, with num_pts for each
	uint16_t *tr_pts, *br_pts, *bl_pts, *tl_pts, tr_num_pts=0U, br_num_pts=0U, bl_num_pts=0U, tl_num_pts=0U;
	tr_pts = (uint16_t*)alloca(s.num_pts * sizeof(uint16_t));
	br_pts = (uint16_t*)alloca(s.num_pts * sizeof(uint16_t));
	bl_pts = (uint16_t*)alloca(s.num_pts * sizeof(uint16_t));
	tl_pts = (uint16_t*)alloca(s.num_pts * sizeof(uint16_t));

	// points that belong to each quadrant are discovered and sorted
	// outer for loop will iterate through the pts array
	for (uint16_t i = 0U; i < s.num_pts; ++i) {
		// we must pull the actual point index from the pts array
		const uint16_t p = s.pts[i];
		const bool top = y_ptr[p] >= y_thresh,
			  right = x_ptr[p] >= x_thresh;
		if (top && right) { tr_pts[tr_num_pts++] = p; } 
		else if (right) { br_pts[br_num_pts++] = p; } 
		else if (top) { tl_pts[tl_num_pts++] = p; } 
		else { bl_pts[bl_num_pts++] = p; }
	}

	// compute chi-square, more efficient not to use pow()
	const float E = s.num_pts*0.25f, chisq = ((tr_num_pts-E)*(tr_num_pts-E) +
		(br_num_pts-E)*(br_num_pts-E) +
		(bl_num_pts-E)*(bl_num_pts-E) +
		(tl_num_pts-E)*(tl_num_pts-E))/E;

	// partition if chi-square or if initial square
	if (chisq > q_thresh || s.num_pts == s.tot_num_pts) {
		const square tr{x_thresh, y_thresh, width*0.5f, tr_pts, tr_num_pts, s.tot_num_pts},
			   br{x_thresh, y_bound1, width*0.5f, br_pts, br_num_pts, s.tot_num_pts},
			   bl{x_bound1, y_bound1, width*0.5f, bl_pts, bl_num_pts, s.tot_num_pts},
			   tl{x_bound1, y_thresh, width*0.5f, tl_pts, tl_num_pts, s.tot_num_pts};

		return calcAPMISplit(x_ptr, y_ptr, tr) + calcAPMISplit(x_ptr, y_ptr, br) + calcAPMISplit(x_ptr, y_ptr, bl) + calcAPMISplit(x_ptr, y_ptr, tl);
	} else {
		// if we don't partition, then we calc MI
		return calcMI(s);
	}
}

/**
 * @brief Calculates the Adaptive Partitioning Mutual Information (APMI)
 * between two vectors. 
 *
 * @param x_vec The first vector.
 * @param y_vec The second vector.
 * @param q_thresh A threshold for chi-square.
 * @param size_thresh A threshold for minimum partition size.
 * @return float The APMI value between the two input vectors.
 */
float calcAPMI(const std::vector<float>& x_vec, const std::vector<float>& y_vec, const float& q_thresh, const uint16_t& size_thresh) {
	// Set file static variables
	::size_thresh = size_thresh;
	::q_thresh = q_thresh;
	
	uint16_t tot_num_pts = x_vec.size();
	
	uint16_t *all_pts = (uint16_t*)alloca(tot_num_pts * sizeof(uint16_t));
	std::iota(all_pts, &all_pts[tot_num_pts], 0U);
	
	// Initialize plane and calc all MIs
	const square init{0.0f, 0.0f, 1.0f, &all_pts[0U], tot_num_pts, tot_num_pts};
	float *x_ptr = (float*)alloca(x_vec.size() * sizeof(float)), *y_ptr = (float*)alloca(y_vec.size() * sizeof(float));
	
	std::copy(x_vec.begin(), x_vec.end(), x_ptr);
	std::copy(y_vec.begin(), y_vec.end(), y_ptr);
	
	return calcAPMISplit(x_ptr, y_ptr, init);
}

/**
 * @brief Calculates the Adaptive Partitioning Conditional Mutual Information
 * (APCMI) between two vectors conditioned on a third vector.
 *
 * This function calculates the Adjusted Partition Conditional Mutual
 * Information (APCMI) by dividing the conditional vector into a specified
 * number of bins. APCMI is a measure of the dependency between two variables
 * given a third one. 
 *
 * @param x_vec The first vector.
 * @param y_vec The second vector.
 * @param cond_vec The conditioned vector.
 * @param n_bins The number of bins to divide the conditioned vector into.
 * @param rand A Mersenne Twister pseudo-random number generator.
 * @param q_thresh Chi-square threshold.
 * @param size_thresh A threshold for minimum partition size.
 * @return float The APCMI value between the two input vectors conditioned on
 * the third vector.
 */
float APCMI(const std::vector<float> &x_vec, std::vector<float> &y_vec, const std::vector<float>& cond_vec, const uint16_t &n_bins, std::mt19937 &rand, const float& q_thresh, const uint16_t& size_thresh) {
	// sort x and y vectors based on binned expression of conditional vector and calculate APCMI for bin
	float APCMI_total = 0.0f;
	for (uint16_t bin = 0U; bin < n_bins; ++bin) {
		std::vector<float> x_vec_bin, y_vec_bin;
		x_vec_bin.reserve(std::ceil(cond_vec.size()/(float) n_bins));
		y_vec_bin.reserve(std::ceil(cond_vec.size()/(float) n_bins));
		
		for (size_t idx = 0U; idx < cond_vec.size(); ++idx)
			if (bin/(float)n_bins < cond_vec[idx] && cond_vec[idx] < (bin+1)/(float)n_bins) {
				x_vec_bin.emplace_back(x_vec[idx]);
				y_vec_bin.emplace_back(y_vec[idx]);
			}
		
		// redo copula-transform for conditional x- and y-vecs
		std::vector<uint16_t> x_idx_ranks = rankIndices(x_vec_bin, rand), y_idx_ranks = rankIndices(y_vec_bin, rand);
		for (uint16_t r = 0; r < x_vec_bin.size(); ++r) {
			x_vec_bin[x_idx_ranks[r]] = (r + 1)/((float)x_vec_bin.size() + 1);
			y_vec_bin[y_idx_ranks[r]] = (r + 1)/((float)y_vec_bin.size() + 1);
		}
		
		// add weighted bin APMI to total APCMI
		APCMI_total += 1/(float)n_bins * calcAPMI(x_vec_bin, y_vec_bin);
	}
	return APCMI_total;
}

/**
 * @brief Calculates the mode of action by assessing the rank correlation
 * between APCMI and the modulator expression.
 * 
 * The function first sorts x and y vectors according to the values of mod_vec.
 * It then performs a copula-transform for the binned conditional x and y
 * vectors and calculates their APMCI. Lastly, a Spearman rank correlation is
 * performed between the APMCI and the modulator to determine the mode of
 * action.
 * 
 * @param x_vec Input x vector of floats
 * @param y_vec Input y vector of floats
 * @param mod_vec Input mode vector of floats
 * @param n_bins Number of bins to divide the vectors
 * @param rand Mersenne twister pseudo-random number generator
 * @param q_thresh
 * @param size_thresh
 * 
 * @return float Mode of Action statistic
 */
float calcModeOfAction(const std::vector<float> &x_vec, const std::vector<float> &y_vec, const std::vector<float>& mod_vec, const uint16_t &n_bins, std::mt19937 &rand, const float& q_thresh, const uint16_t& size_thresh) {
	
	// add x and y APCMI to vector based on sorted values of mod_vec
	std::vector<float> APCMIs;
	APCMIs.reserve(n_bins);
	for (uint16_t bin = 0U; bin < n_bins; ++bin) {
		std::vector<float> x_vec_bin, y_vec_bin;
		x_vec_bin.reserve(std::ceil(mod_vec.size()/(float) n_bins));
		y_vec_bin.reserve(std::ceil(mod_vec.size()/(float) n_bins));
		
		for (size_t idx = 0U; idx < mod_vec.size(); ++idx)
			if (bin/(float)n_bins < mod_vec[idx] && mod_vec[idx] < (bin+1)/(float)n_bins) {
				x_vec_bin.emplace_back(x_vec[idx]);
				y_vec_bin.emplace_back(y_vec[idx]);
			}
		
		// redo copula-transform for conditional x- and y-vecs
		std::vector<uint16_t> x_idx_ranks = rankIndices(x_vec_bin, rand), y_idx_ranks = rankIndices(y_vec_bin, rand);
		for (uint16_t r = 0; r < x_vec_bin.size(); ++r) {
			x_vec_bin[x_idx_ranks[r]] = (r + 1)/((float)x_vec_bin.size() + 1);
			y_vec_bin[y_idx_ranks[r]] = (r + 1)/((float)y_vec_bin.size() + 1);
		}
		
		// add weighted bin APMI to total APCMI
		APCMIs.emplace_back(calcAPMI(x_vec_bin, y_vec_bin));
	}
	
	// create rank vectors
	std::vector<uint16_t> mod_bin_ranks(n_bins, 0U), APCMIs_ranks(n_bins, 0U);
	std::iota(mod_bin_ranks.begin(), mod_bin_ranks.end(), 1U);
	std::vector<uint16_t> APCMIs_idx_ranks = rankIndices(APCMIs, rand);
	for (uint16_t r = 0U; r < APCMIs_idx_ranks.size(); ++r)
		APCMIs_ranks[APCMIs_idx_ranks[r]] = r + 1;
	
	return spearmanCorrelate(mod_bin_ranks, APCMIs_ranks);
}

#include "ARACNe3.hpp"
#include "algorithms.hpp"

extern float DEVELOPER_mi_cutoff;
static float q_thresh;
static uint16_t size_thresh;

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

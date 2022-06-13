/*
 This file is intended to contain functions callable from R that will allow for Rcpp usage and a custom analysis.  Right now, this file is not ready for release.  However, one might find these functions useful during ARACNe3 development.
 */

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp17")]]

using namespace Rcpp;

static float q_thresh;
static vector<float> *vec_x, *vec_y;
static vector<float> mis;
static uint16_t tot_num_pts, size_thresh;

/*
 * Calculate the MI for a square struct
 */
float calcMI(const square &s) {
	const float pxy = s.num_pts/(float)tot_num_pts, marginal = s.width;
	float mi;
	return isfinite(mi = pxy*log(pxy/marginal/marginal)) ? mi : 0.0;
}

/*
 * Recursive tessellation of XY plane and mi calculation at dead-ends
 *
 * s a pointer to the square we are considering to partition
 * mis a pointer to a vector of mis
 *
 * returns nothing; values computed from pointers to original
 */
void APMI_split(const square &s) {
	// extract values; memory disadvantage but runtime advantage
	const float x_bound1=s.x_bound1, y_bound1=s.y_bound1, width=s.width;
	const uint16_t *pts=s.pts, num_pts=s.num_pts;

	// if we have less points in the square than size_thresh, calc MI
	if (num_pts < size_thresh) {mis.push_back(calcMI(s)); return;}

	// thresholds for potential partition of XY plane
	const float x_thresh = x_bound1 + width*0.5,
		  y_thresh = y_bound1 + width*0.5;

	// indices for quadrants, to test chi-square, with num_pts for each
	uint16_t tr_pts[num_pts], br_pts[num_pts], bl_pts[num_pts],
		  tl_pts[num_pts], tr_num_pts=0, br_num_pts=0, bl_num_pts=0,
		  tl_num_pts=0;

	// points that belong to each quadrant are discovered and sorted
	// outer for loop will iterate through the pts array
	for (uint16_t i = 0; i < num_pts; ++i) {
		// we must pull the actual point index from the pts array
		const uint16_t p = pts[i];
		const bool top = (*vec_y)[p] >= y_thresh,
			  right = (*vec_x)[p] >= x_thresh;
		if (top && right) { tr_pts[tr_num_pts++] = p; }
		else if (right) { br_pts[br_num_pts++] = p; }
		else if (top) { tl_pts[tl_num_pts++] = p; }
		else { bl_pts[bl_num_pts++] = p; }
	}

	// compute chi-square, more efficient not to use pow()
	const float E = num_pts*0.25, chisq = ((tr_num_pts-E)*(tr_num_pts-E) +
		(br_num_pts-E)*(br_num_pts-E) +
		(bl_num_pts-E)*(bl_num_pts-E) +
		(tl_num_pts-E)*(tl_num_pts-E))/E;

	// partition if chi-square or if initial square
	if (chisq > q_thresh || num_pts == tot_num_pts) {
		const square tr{x_thresh, y_thresh, width*0.5, tr_pts, tr_num_pts},
			   br{x_thresh, y_bound1, width*0.5, br_pts, br_num_pts},
			   bl{x_bound1, y_bound1, width*0.5, bl_pts, bl_num_pts},
			   tl{x_bound1, y_thresh, width*0.5, tl_pts, tl_num_pts};

		APMI_split(tr);
		APMI_split(br);
		APMI_split(bl);
		APMI_split(tl);
	} else {
		// if we don't partition, then we calc MI
		mis.push_back(calcMI(s));
	}

	return;
}

/* Takes in two expression vectors (regulator-target) and computes APMI for
 * each partition of the observation space
 *
 * vec_x an X expression vector
 * vec_y a Y expression vector
 * q_thresh q-value for chi-square bin independence
 * size_thresh minimum points in a tile to consider chi-square and partition
 *
 * returns the APMI
 */
// [[Rcpp::export]]
float APMI(vector<float> &vec_x, vector<float> &vec_y,
		const float q_thresh,
		const uint16_t size_thresh) {
	// Set file static variables
	::size_thresh = size_thresh;
	::q_thresh = q_thresh;
	::vec_x = &vec_x;
	::vec_y = &vec_y;
	::tot_num_pts = (*::vec_x).size();


	// Make an array of all indices, to be partitioned later
	uint16_t all_pts[(*::vec_x).size()];
	for (uint16_t i = 0; i < tot_num_pts; i++) { all_pts[i] = i; }
	
	// Initialize plane and calc all MIs
	const square init{0.0, 0.0, 1.0,  all_pts, tot_num_pts};
	APMI_split(init);
	
	float mi = std::accumulate(mis.begin(), mis.end(), static_cast<float>(0.0));
	mis.clear();
	
	return mi;
}

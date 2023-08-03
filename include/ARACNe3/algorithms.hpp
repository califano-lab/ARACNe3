#ifndef algorithms_hpp
#define algorithms_hpp

#include <vector>
#include <random>

typedef struct {const float &x_bound1, &y_bound1, &width; uint16_t *const pts, &num_pts, &tot_num_pts; } square;

float spearmanCorrelate(const std::vector<uint16_t>& x_ranked, const std::vector<uint16_t>& y_ranked);
float fishersMethodP(const std::vector<float>& p_vals);
float calcAPMI(const std::vector<float>& x_vec, const std::vector<float>& y_vec, const float& q_thresh = 7.815, const uint16_t& size_thresh = 4);
std::vector<uint16_t> rankIndices(const std::vector<float>& vec, std::mt19937 &rand);
float APCMI(const std::vector<float> &x_vec, std::vector<float> &y_vec, const std::vector<float>& cond_vec, const uint16_t &n_bins, std::mt19937 &rand, const float &q_thresh = 7.815f, const uint16_t &size_thresh = 4U);
std::pair<float, float> linearRegress(const std::vector<float>& x, const std::vector<float>& y);
float calcModeOfAction(const std::vector<float> &x_vec, const std::vector<float> &y_vec, const std::vector<float>& mod_vec, const uint16_t &n_bins, std::mt19937 &rand, const float& q_thresh = 7.815f, const uint16_t& size_thresh = 4U);

#endif /* algorithms_h */

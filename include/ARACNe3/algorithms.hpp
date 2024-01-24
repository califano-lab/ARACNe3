#pragma once

#include <random>
#include <vector>

typedef struct {
  const float x_bound1, y_bound1, width;
  uint16_t *const pts, num_pts, tot_num_pts;
} square;

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
std::vector<uint16_t> rankIndices(const std::vector<float> &vec,
                                  std::mt19937 &rand);

/**
 * @brief Performs a copula transform on a given vector of floats.
 *
 * The copula transform ranks the elements of the input vector, and then
 * normalizes these ranks by dividing each by the total number of elements in
 * the vector + 1. Ties are broken by random shuffling
 *
 * @param data A const reference to a vector of floats to be transformed.
 * @param rand A reference to a random device for random tie breaking.
 * @return std::vector<float> A vector containing the copula transform of the
 * input.
 */
std::vector<float> copulaTransform(const std::vector<float> &data,
                                   std::mt19937 &rand);

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
float calcAPMI(const std::vector<float> &x_vec, const std::vector<float> &y_vec,
               const float q_thresh = 7.815, const uint16_t size_thresh = 4);

float calcSCC(const std::vector<uint16_t> &x_ranked,
              const std::vector<uint16_t> &y_ranked);

std::pair<float, float> OLS(const std::vector<float> &x,
                            const std::vector<float> &y);

double lchoose(const uint16_t &n, const uint16_t &k);

double rightTailBinomialP(const uint16_t n, const uint16_t k,
                          const float theta);
/**
 * Calculates the right-tail binomial probability density while remaining in log
 * space. Prevents underflow
 */
double lRightTailBinomialP(const uint16_t n, const uint16_t k,
                           const float theta);

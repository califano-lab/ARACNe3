#include "algorithms.hpp"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <stdexcept>

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
  const float pxy = s.num_pts / (float)s.tot_num_pts, marginal = s.width,
              mi = pxy * std::log(pxy / (marginal * marginal));
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
const float calcAPMISplit(const float *const x_ptr, const float *const y_ptr,
                          const square s) {
  // if we have less points in the square than size_thresh, calc MI
  if (s.num_pts < size_thresh) {
    return calcMI(s);
  }

  // thresholds for potential partition of XY plane
  const float x_thresh = s.x_bound1 + s.width * 0.5f,
              y_thresh = s.y_bound1 + s.width * 0.5f;

  // indices for quadrants, to test chi-square, with num_pts for each
  uint16_t *tr_pts, *br_pts, *bl_pts, *tl_pts, tr_num_pts = 0U, br_num_pts = 0U,
                                               bl_num_pts = 0U, tl_num_pts = 0U;
  tr_pts = (uint16_t *)alloca(s.num_pts * sizeof(uint16_t));
  br_pts = (uint16_t *)alloca(s.num_pts * sizeof(uint16_t));
  bl_pts = (uint16_t *)alloca(s.num_pts * sizeof(uint16_t));
  tl_pts = (uint16_t *)alloca(s.num_pts * sizeof(uint16_t));

  // points that belong to each quadrant are discovered and sorted
  // outer for loop will iterate through the pts array
  for (uint16_t i = 0U; i < s.num_pts; ++i) {
    // we must pull the actual point index from the pts array
    const uint16_t p = s.pts[i];
    const bool top = y_ptr[p] >= y_thresh, right = x_ptr[p] >= x_thresh;
    if (top && right) {
      tr_pts[tr_num_pts++] = p;
    } else if (right) {
      br_pts[br_num_pts++] = p;
    } else if (top) {
      tl_pts[tl_num_pts++] = p;
    } else {
      bl_pts[bl_num_pts++] = p;
    }
  }

  // compute chi-square, more efficient not to use pow()
  const float E = s.num_pts * 0.25f,
              chisq = ((tr_num_pts - E) * (tr_num_pts - E) +
                       (br_num_pts - E) * (br_num_pts - E) +
                       (bl_num_pts - E) * (bl_num_pts - E) +
                       (tl_num_pts - E) * (tl_num_pts - E)) /
                      E;

  // partition if chi-square or if initial square
  if (chisq > q_thresh || s.num_pts == s.tot_num_pts) {
    const square tr{x_thresh, y_thresh,   s.width * 0.5f,
                    tr_pts,   tr_num_pts, s.tot_num_pts},
        br{x_thresh, s.y_bound1, s.width * 0.5f,
           br_pts,   br_num_pts, s.tot_num_pts},
        bl{s.x_bound1, s.y_bound1, s.width * 0.5f,
           bl_pts,     bl_num_pts, s.tot_num_pts},
        tl{s.x_bound1, y_thresh,   s.width * 0.5f,
           tl_pts,     tl_num_pts, s.tot_num_pts};

    return calcAPMISplit(x_ptr, y_ptr, tr) + calcAPMISplit(x_ptr, y_ptr, br) +
           calcAPMISplit(x_ptr, y_ptr, bl) + calcAPMISplit(x_ptr, y_ptr, tl);
  } else {
    // if we don't partition, then we calc MI
    return calcMI(s);
  }
}

float calcAPMI(const std::vector<float> &x_vec, const std::vector<float> &y_vec,
               const float q_thresh, const uint16_t size_thresh) {
  // Set file static variables
  ::size_thresh = size_thresh;
  ::q_thresh = q_thresh;

  uint16_t tot_num_pts = x_vec.size();

  uint16_t *all_pts = (uint16_t *)alloca(tot_num_pts * sizeof(uint16_t));
  std::iota(all_pts, &all_pts[tot_num_pts], 0U);

  // Initialize plane and calc all MIs
  const square init{0.0f, 0.0f, 1.0f, &all_pts[0U], tot_num_pts, tot_num_pts};
  float *x_ptr = (float *)alloca(x_vec.size() * sizeof(float)),
        *y_ptr = (float *)alloca(y_vec.size() * sizeof(float));

  std::copy(x_vec.begin(), x_vec.end(), x_ptr);
  std::copy(y_vec.begin(), y_vec.end(), y_ptr);

  return calcAPMISplit(x_ptr, y_ptr, init);
}

std::vector<uint32_t> rankWithRandomTiebreak(const std::vector<float> &vec,
                                             std::mt19937 &rnd) {
  const std::size_t n = vec.size();

  // We need to sort an indexes vector, but it should be paired with another
  // random indexes vector for breaking ties
  std::vector<std::size_t> idx_ranks(n), shuffled_priorities(n);
  std::iota(idx_ranks.begin(), idx_ranks.end(), 0u);
  std::iota(shuffled_priorities.begin(), shuffled_priorities.end(), 0u);
  std::shuffle(shuffled_priorities.begin(), shuffled_priorities.end(), rnd);

  std::vector<std::pair<std::size_t, std::size_t>> idx_ranks_p(n);
  for (uint32_t i = 0u; i < n; ++i)
    idx_ranks_p[i] = {idx_ranks[i], shuffled_priorities[i]};

  std::sort(idx_ranks_p.begin(), idx_ranks_p.end(), [&](auto &p1, auto &p2) {
    return vec[p1.first] == vec[p2.first] ? p1.second < p2.second
                                          : vec[p1.first] < vec[p2.first];
  });

  std::vector<uint32_t> ranks(n);
  for (uint32_t r = 0u; r < n; ++r)
    ranks[idx_ranks_p[r].first] = r + 1u;

  return ranks;
}

std::vector<float> rankWithAverageTiebreak(const std::vector<float> &v) {
  std::vector<uint32_t> idx_rnks(v.size());
  std::iota(idx_rnks.begin(), idx_rnks.end(), 0u);
  std::sort(idx_rnks.begin(), idx_rnks.end(),
            [&](uint32_t i1, uint32_t i2) { return v[i1] < v[i2]; });

  // This will assign ranks like 1 2 2 4 5 6 6 6 9 (truncated ranks for ties)
  std::vector<float> ranks(v.size());
  for (uint32_t i = 0u; i < idx_rnks.size(); ++i) {
    if (i > 0u && v[idx_rnks[i]] == v[idx_rnks[i - 1]])
      ranks[idx_rnks[i]] = ranks[idx_rnks[i - 1]];
    else
      ranks[idx_rnks[i]] = i + 1u;
  }

  // Now we need to count how many are the same, add, and average
  auto cur_beg = idx_rnks.begin();
  auto cur_end = cur_beg + 1u;

  // Sum of 0 + 1 + 2 (the truncated portion of ranks r r r) is n(n - 1)/2
  while (cur_end != idx_rnks.end()) {
    // base case is rank goes 1 2 3, in which sum is 0 and ranks don't change
    if (ranks[*cur_end] > ranks[*cur_beg]) {
      double r = static_cast<double>(ranks[*cur_beg]);
      double n = static_cast<double>(cur_end - cur_beg);
      double sum = n * (n - 1.) / 2.;
      float avg_rank = (n * r + sum) / n;

      for (auto it = cur_beg; it != cur_end; ++it)
        ranks[*it] = avg_rank;

      // set the new cur_beg and cur_end
      cur_beg = cur_end;
      ++cur_end;
    } else { // when ranks is the same, we increment cur_end
      ++cur_end;
    }
  }
  // final transformation for cur_beg -> actual end
  double r = static_cast<double>(ranks[*cur_beg]);
  double n = static_cast<double>(cur_end - cur_beg);
  double sum = n * (n - 1.) / 2.;
  float avg_rank = (n * r + sum) / n;

  for (auto it = cur_beg; it != cur_end; ++it)
    ranks[*it] = avg_rank;

  return ranks;
}

float pearsonsR(const std::vector<float> &x_vec,
                const std::vector<float> &y_vec) {
  const std::size_t n = x_vec.size();

  if (x_vec.size() != y_vec.size())
    throw std::runtime_error(
        "Cannot perform regression on vectors of unequal size");

  float x_mean =
      std::reduce(x_vec.cbegin(), x_vec.cend()) / static_cast<float>(n);
  float y_mean =
      std::reduce(y_vec.cbegin(), y_vec.cend()) / static_cast<float>(n);

  float ssr_x = std::transform_reduce(
      x_vec.cbegin(), x_vec.cend(), 0.f, std::plus<>(),
      [x_mean](float x) { return (x - x_mean) * (x - x_mean); });
  float ssr_y = std::transform_reduce(
      y_vec.cbegin(), y_vec.cend(), 0.f, std::plus<>(),
      [y_mean](float y) { return (y - y_mean) * (y - y_mean); });

  float sum_prod = 0.f;
  for (std::size_t i = 0U; i < n; ++i)
    sum_prod += (x_vec[i] - x_mean) * (y_vec[i] - y_mean);

  return sum_prod / std::sqrt(ssr_x * ssr_y);
}

float spearmansRho(const std::vector<float> &x_vec,
                   const std::vector<float> &y_vec) {
  const double n = x_vec.size();

  std::vector<float> x_ranks = rankWithAverageTiebreak(x_vec);
  std::vector<float> y_ranks = rankWithAverageTiebreak(y_vec);

  return pearsonsR(x_ranks, y_ranks);
}

double lChoose(const uint16_t &n, const uint16_t &k) {
  return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
}

double rightTailBinomialP(uint16_t n, uint16_t k, float theta) {
  // If k is 0, the right-tail probability includes all possibilities, which sum
  // to 1.
  if (k == 0)
    return 1.0;

  double p = 0.0;
  // Start from k and go up to n to avoid underflow.
  for (uint16_t i = k; i <= n; ++i)
    p += std::exp(lChoose(n, i) + i * std::log(theta) +
                  (n - i) * std::log(1 - theta));

  return p;
}

double lRightTailBinomialP(uint16_t n, uint16_t k, float theta) {
  // If k is 0, the right-tail probability includes all possibilities, which is
  // p = 1.
  if (k == 0)
    return std::log(1.);

  double max_log_p = -std::numeric_limits<double>::infinity();
  std::vector<double> log_ps;

  // Calculate log probabilities and find the maximum log probability
  for (uint16_t i = k; i <= n; ++i) {
    double log_p =
        lChoose(n, i) + i * std::log(theta) + (n - i) * std::log(1. - theta);
    max_log_p = std::max(max_log_p, log_p);
    log_ps.push_back(log_p);
  }

  // Apply log-sum-exp trick to sum in log domain
  double log_sum_exp = 0.;
  for (double log_p : log_ps) {
    log_sum_exp += std::exp(log_p - max_log_p);
  }

  return max_log_p + std::log(log_sum_exp);
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
std::pair<float, float> OLS(const std::vector<float> &x_vec,
                            const std::vector<float> &y_vec) {
  const std::size_t n = x_vec.size();

  if (x_vec.size() != y_vec.size())
    throw std::runtime_error(
        "Cannot perform regression on vectors of unequal size");

  float x_mean =
      std::reduce(x_vec.cbegin(), x_vec.cend()) / static_cast<float>(n);
  float y_mean =
      std::reduce(y_vec.cbegin(), y_vec.cend()) / static_cast<float>(n);

  float ssr_x = std::transform_reduce(
      x_vec.cbegin(), x_vec.cend(), 0.f, std::plus<>(),
      [x_mean](float x) { return (x - x_mean) * (x - x_mean); });

  float sum_prod = 0.f;
  for (std::size_t i = 0U; i < n; ++i)
    sum_prod += (x_vec[i] - x_mean) * (y_vec[i] - y_mean);

  float slope = sum_prod / ssr_x;
  float intercept = y_mean - slope * x_mean;

  return {slope, intercept};
}

std::vector<float> copulaTransform(const std::vector<float> &data,
                                   std::mt19937 &rnd) {
  std::size_t n = data.size();

  std::vector<uint32_t> ranks = rankWithRandomTiebreak(data, rnd);
  std::vector<float> copulas(n);

  std::transform(ranks.cbegin(), ranks.cend(), copulas.begin(),
                 [=](auto &r) { return r / static_cast<float>(n); });

  return copulas;
}

float estimateFPRNoMaxEnt(const float alpha, const std::string &method,
                          const uint32_t num_edges_after_threshold_pruning,
                          const uint32_t tot_possible_edges) {
  if (method == "FDR")
    return alpha * num_edges_after_threshold_pruning /
           (tot_possible_edges -
            (1 - alpha) * num_edges_after_threshold_pruning);
  else if (method == "FWER")
    return alpha / tot_possible_edges;
  else // method == "FPR"
    return alpha;
}

float estimateFPRWithMaxEnt(const float alpha, const std::string &method,
                            const uint32_t num_edges_after_threshold_pruning,
                            const uint32_t num_edges_after_MaxEnt_pruning,
                            const uint32_t tot_possible_edges) {
  if (method == "FDR")
    return alpha * num_edges_after_MaxEnt_pruning /
           (tot_possible_edges -
            (1 - alpha) * num_edges_after_threshold_pruning);
  else if (method == "FWER")
    return alpha / tot_possible_edges * num_edges_after_MaxEnt_pruning /
           num_edges_after_threshold_pruning;
  else // method == "FPR"
    return alpha * num_edges_after_MaxEnt_pruning /
           num_edges_after_threshold_pruning;
}

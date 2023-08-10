#pragma once

#include "ARACNe3.hpp"
#include <random>
#include <vector>

typedef struct {
  const float &x_bound1, &y_bound1, &width;
  uint16_t *const pts, &num_pts, &tot_num_pts;
} square;

std::vector<uint16_t> rankIndices(const std::vector<float> &vec,
                                   std::mt19937 &rand);

float calcAPMI(const std::vector<float> &x_vec, const std::vector<float> &y_vec,
               const float q_thresh = 7.815, const uint16_t size_thresh = 4);

float calcSCC(const std::vector<uint16_t> &x_ranked,
              const std::vector<uint16_t> &y_ranked);

std::pair<float, float> linearRegress(const std::vector<float> &x,
                                      const std::vector<float> &y);

double lchoose(const uint16_t &n, const uint16_t &k);

double rightTailBinomialP(const uint16_t n, const uint16_t k,
                          const float theta);

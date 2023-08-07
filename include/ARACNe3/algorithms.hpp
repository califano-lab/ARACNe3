#pragma once

#include <vector>
#include "ARACNe3.hpp"

typedef struct {
  const float &x_bound1, &y_bound1, &width;
  uint16_t *const pts, &num_pts, &tot_num_pts;
} square;

float calcAPMI(const std::vector<float> &x_vec, const std::vector<float> &y_vec,
               const float &q_thresh = 7.815, const uint16_t &size_thresh = 4);

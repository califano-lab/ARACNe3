#pragma once

#include <vector>
#include <tuple>
#include <cstddef>

class APMINullModel {
public:
  APMINullModel();
  APMINullModel(const std::vector<float> &nm, float om, float ob, std::size_t ns,
                std::size_t nn, uint32_t s)
      : null_mis(nm), ols_m(om), ols_b(ob), n_samps(ns), n_nulls(nn), seed(s){};
  APMINullModel(const std::size_t n_samps, const std::size_t n_nulls,
                const uint32_t seed);

  float getMIPVal(const float mi, const float p_precise = 0.001f) const;

  std::tuple<std::vector<float>, float, float, std::size_t, std::size_t, uint32_t>
  getModel() {
    return {null_mis, ols_m, ols_b, n_samps, n_nulls, seed};
  }

private:
  std::vector<float> null_mis;
  float ols_m, ols_b;
  std::size_t n_samps, n_nulls;
  uint32_t seed;
};

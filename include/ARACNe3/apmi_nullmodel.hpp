#pragma once

#include <random>
#include <string>
#include <vector>

class APMINullModel {
private:
  std::vector<float> null_mis;
  float m, b;
  std::string nulls_filename_no_extension, OLS_coefs_filename_no_extension;

public:
  APMINullModel(const APMINullModel &copied); // copy ctor
  // rand should be passed from main based on seed for predictable behavior.
  APMINullModel(const uint32_t n_nulls, const uint16_t tot_num_subsample,
                const std::string &cached_dir, std::mt19937 &rand);
  ~APMINullModel();
  void cacheNullModel(const std::string cached_dir); // cache vec, m, and b
  const float
  getMIPVal(const float &mi,
            const float &p_precise = 0.001f) const; // return p value
};

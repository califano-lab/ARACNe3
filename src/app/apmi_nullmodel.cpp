#include "apmi_nullmodel.hpp"
#include "ARACNe3.hpp"
#include "algorithms.hpp"
#include <filesystem>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <omp.h>

extern uint16_t nthreads;

APMINullModel::APMINullModel(const APMINullModel &copied) {
  null_mis = copied.null_mis;
  m = copied.m;
  b = copied.b;
  nulls_filename_no_extension = copied.nulls_filename_no_extension;
  OLS_coefs_filename_no_extension = copied.OLS_coefs_filename_no_extension;
}

/*
 Computes 1 million null mutual information values for the sample size.  Checks
 whether there already exists a null_mi vector (nulls_filename) in the cached
 directory.
 */
APMINullModel::APMINullModel(const uint32_t n_nulls,
                             const uint16_t tot_num_subsample,
                             const std::string &cached_dir,
                             std::mt19937 &rand) {
  this->nulls_filename_no_extension = "Nssamp-" +
                                      std::to_string(tot_num_subsample) +
                                      "_Nnull-" + std::to_string(n_nulls);
  this->OLS_coefs_filename_no_extension = nulls_filename_no_extension + "_OLS";

#ifdef _DEBUG // If debug, we must generate a new null model each time
  constexpr bool debug = true;
#else
  constexpr bool debug = false;
#endif

  if (std::filesystem::exists(cached_dir + nulls_filename_no_extension +
                              ".txt") &&
      std::filesystem::exists(cached_dir + OLS_coefs_filename_no_extension +
                              ".txt") &&
      !debug) {
    this->null_mis.reserve(n_nulls);

    std::ifstream nulls_file(cached_dir + nulls_filename_no_extension + ".txt",
                             std::ios::in | std::ios::binary);
    std::ifstream OLS_coef_file(cached_dir + OLS_coefs_filename_no_extension +
                                    ".txt",
                                std::ios::in | std::ios::binary);
    std::istream_iterator<float> nulls_iterator(nulls_file);
    std::istream_iterator<float> OLS_iterator(OLS_coef_file);

    for (uint32_t i = 0; i < n_nulls; ++i)
      null_mis.emplace_back(*nulls_iterator++);
    this->m = *OLS_iterator++;
    this->b = *OLS_iterator;
  } else {
    // make the ref vector for null APMI against shuffled version
    std::vector<float> ref_vec;
    ref_vec.reserve(tot_num_subsample);

    for (uint16_t i = 1U; i <= tot_num_subsample; ++i)
      ref_vec.emplace_back(((float)i) / (tot_num_subsample + 1));

    std::vector<float> shuffle_vec = ref_vec;

    this->null_mis = std::vector<float>(n_nulls);

#pragma omp parallel for num_threads(nthreads)
    for (uint32_t i = 0U; i < n_nulls; ++i) {
#pragma omp critical
      { std::shuffle(shuffle_vec.begin(), shuffle_vec.end(), rand); }
      null_mis[i] = calcAPMI(ref_vec, shuffle_vec);
    }

    // sort largest to smallest
    std::sort(null_mis.begin(), null_mis.end(), std::greater<float>());

    // OLS regress log(p) vs MI for eCDF p < 0.01
    uint32_t significant_thresh_idx =
        std::ceil(n_nulls * 1.0f / 100.0f); // index of p = 0.01
    std::vector<float> significant_mis(
        null_mis.begin(), null_mis.begin() + significant_thresh_idx);
    std::vector<float> significant_mi_ps(significant_thresh_idx);
    for (uint32_t i = 0; i < significant_thresh_idx; ++i)
      significant_mi_ps[i] = ((i + 1) / (float)n_nulls); // fill p-vals
    std::transform(significant_mi_ps.begin(),
                   significant_mi_ps.begin() + significant_thresh_idx,
                   significant_mi_ps.begin(), [](const auto &p) -> float {
                     return std::log(p);
                   }); // log-transform

    std::pair<float, float> sol =
        linearRegress(significant_mis, significant_mi_ps);
    this->m = sol.first;
    this->b = sol.second;
  }
}

APMINullModel::~APMINullModel() {}

void APMINullModel::cacheNullModel(const std::string cached_dir) {
  std::string nulls_filename =
      cached_dir + nulls_filename_no_extension + ".txt";
  std::string ols_filename =
      cached_dir + OLS_coefs_filename_no_extension + ".txt";

  if (!std::filesystem::exists(nulls_filename) ||
      !std::filesystem::exists(ols_filename)) {
    std::ofstream nulls_file(nulls_filename, std::ios::out | std::ios::binary);
    std::ofstream OLS_coefs_file(ols_filename,
                                 std::ios::out | std::ios::binary);
    for (auto it = null_mis.cbegin(); it != null_mis.cend(); ++it)
      nulls_file << *it << '\n';
    OLS_coefs_file << m << '\n' << b << '\n';
  }
  return;
}

const float APMINullModel::getMIPVal(const float &mi,
                                     const float &p_precise) const {
  // points to first index for which mi > the rest.
  auto it = std::upper_bound(null_mis.cbegin(), null_mis.cend(), mi,
                             std::greater<float>());

  // p-value as a percentile.  We add 1 because it is an index
  const float p = (it - null_mis.cbegin() + 1) / (float)null_mis.size();

  if (p < p_precise)
    return std::min(p, std::exp(m * mi + b)); // invert log
  else
    return p;
}

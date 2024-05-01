#include <algorithm>
#include <random>
#include <fstream>

#include "apmi_nullmodel.hpp"
#include "algorithm.hpp"

APMINullModel::APMINullModel() {};

APMINullModel::APMINullModel(const std::size_t n_samps, const std::size_t n_nulls,
                             const uint32_t seed)
    : n_samps(n_samps), n_nulls(n_nulls), seed(seed) {
  std::mt19937 rnd(seed);

  // All null APMIs are based on a shuffled copula of ref_vec with n_samps
  std::vector<float> ref_vec;
  ref_vec.reserve(n_samps);
  for (std::size_t i = 0u; i < n_samps; ++i)
    ref_vec.emplace_back((i + 1.f) / (n_samps + 1.f));

  std::vector<std::vector<float>> shuffle_vecs(n_nulls, ref_vec);
  for (std::size_t i = 0u; i < n_nulls; ++i)
    std::shuffle(shuffle_vecs[i].begin(), shuffle_vecs[i].end(), rnd);

  this->null_mis = std::vector<float>(n_nulls);
  for (std::size_t i = 0u; i < n_nulls; ++i)
    this->null_mis[i] = calcAPMI(ref_vec, shuffle_vecs[i]);

  // sort largest to smallest
  std::sort(null_mis.begin(), null_mis.end(), std::greater<float>());

  // OLS regress log(p) vs MI for eCDF p < 0.01 (but use for p < 1e-4)
  float sig_p = 0.01f;
  std::size_t sig_p_offendidx = std::floor(n_nulls * sig_p);
  std::vector<float> sig_mis(null_mis.begin(),
                             null_mis.begin() + sig_p_offendidx);
  std::vector<float> sig_mi_ps(sig_mis.size());
  for (std::size_t i = 0u; i < sig_mi_ps.size(); ++i)
    sig_mi_ps[i] = (i + 1.f) / (n_nulls + 1.f); // fill p-vals
  std::transform(sig_mi_ps.begin(), sig_mi_ps.end(), sig_mi_ps.begin(),
                 [](const float p) { return std::log(p); }); // log-transform

  const auto [m, b] = OLS(sig_mis, sig_mi_ps);
  this->ols_m = m;
  this->ols_b = b;
}

float APMINullModel::getMIPVal(const float mi, const float p_precise) const {
  // points to the index where comp is true; i.e. off end of [elements] > mi.
  std::size_t n_nulls_gt_mi = std::upper_bound(null_mis.cbegin(), null_mis.cend(),
                                          mi, std::greater<float>()) -
                         null_mis.cbegin();

  // p-value as a percentile.
  const float p = (n_nulls_gt_mi + 1.f) / (null_mis.size() + 1.f);

  return p < p_precise ? std::exp(ols_m * mi + ols_b) : p;
}

template<class Archive>
void APMINullModel::serialize(Archive& ar, const unsigned int version)
{
    ar & null_mis;
    ar & ols_m;
    ar & ols_b;
    ar & n_samps;
    ar & n_nulls;
    ar & seed;
}


APMINullModel
APMINullModel::getCachedModel(const std::string &cached_blob_name) {
  APMINullModel apmi_null_model;
  std::ifstream ifs(cached_blob_name, std::ios::binary);
  if (ifs) {
    try {
      boost::archive::binary_iarchive ia(ifs);
      ia >> apmi_null_model;
    } catch (const std::exception &e) {
      throw std::runtime_error("Failed to get cached APMINullModel: " +
                               std::string(e.what()));
    }
  } else {
    throw std::runtime_error("File not found: " + cached_blob_name);
  }
  return apmi_null_model;
}

void APMINullModel::cacheModel(const std::string &cached_blob_name) const {
  std::ofstream ofs(cached_blob_name, std::ios::binary);
  if (ofs) {
    try {
      boost::archive::binary_oarchive oa(ofs);
      oa << *this;
    } catch (const std::exception &e) {
      throw std::runtime_error("Failed to serialize APMINullModel: " +
                               std::string(e.what()));
    }
  } else {
    throw std::runtime_error("Unable to open file for writing: " +
                             cached_blob_name);
  }
}

#pragma once

#include <vector>
#include <tuple>
#include <cstddef>
#include <cstdint>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>  // for vector serialize member

class APMINullModel {
public:
  APMINullModel();
  APMINullModel(const std::vector<float> &nm, float om, float ob, std::size_t ns,
                std::size_t nn, uint32_t s)
      : null_mis(nm), ols_m(om), ols_b(ob), n_samps(ns), n_nulls(nn), seed(s){};
  APMINullModel(const std::size_t n_samps, const std::size_t n_nulls,
                const uint32_t seed);

  float getMIPVal(const float mi, const float p_precise = 1e-4f) const;

  std::tuple<std::vector<float>, float, float, std::size_t, std::size_t, uint32_t>
  getModel() {
    return {null_mis, ols_m, ols_b, n_samps, n_nulls, seed};
  }

  static APMINullModel getCachedModel(const std::string &cached_blob_name);
  void cacheModel(const std::string& cached_blob_name) const;

private:
  std::vector<float> null_mis;
  float ols_m, ols_b;
  std::size_t n_samps, n_nulls;
  uint32_t seed;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);

  friend class boost::serialization::access;
};

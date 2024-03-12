#include <gtest/gtest.h> // resolved by CMake
#include "apmi_nullmodel.hpp"

TEST(APMINullModelTest, Initialization) {
  const std::size_t n_samps = 1000u;
  const std::size_t n_nulls = 10'000u;
  const uint32_t seed = 123u; // Use a fixed seed for reproducibility
  APMINullModel model(n_samps, n_nulls, seed);

  auto [null_mis, ols_m, ols_b, n_samps_nm, n_nulls_nm, seed_nm] =
        model.getModel();

  for (auto it = null_mis.cbegin()+1u; it != null_mis.cend(); ++it)
    EXPECT_GE(*(it-1u), *it);

  EXPECT_EQ(n_samps, n_samps_nm);
  EXPECT_EQ(n_nulls, n_nulls_nm);
  EXPECT_EQ(seed, seed_nm);
}

TEST(APMINullModelTest, PValueForExtremeMIScores) {
  const std::size_t n_samps = 1000u;
  const std::size_t n_nulls = 10'000u;
  const uint32_t seed = 123u;
  APMINullModel model(n_samps, n_nulls, seed);

  auto [null_mis, ols_m, ols_b, n_samps_nm, n_nulls_nm, seed_nm] =
      model.getModel();

  // Assuming getMIPVal is public or can be tested indirectly
  // Low MI score, should have p-value close to 1
  float low_mi = *std::min_element(null_mis.begin(), null_mis.end()) - 0.1f;
  float p_precise = 1e-4f;
  EXPECT_NEAR(1.0f, model.getMIPVal(low_mi, p_precise), 0.05f);

  // High MI score, should have p-value close to 0, both with p_precise and not
  float high_mi = *std::max_element(null_mis.begin(), null_mis.end()) + 0.1f;
  EXPECT_NEAR(0.0f, model.getMIPVal(high_mi, p_precise), 0.05f);

  EXPECT_NEAR(0.0f, model.getMIPVal(high_mi, p_precise = 0.f), 0.05f);
}

#include <filesystem>
#include <gtest/gtest.h>

#include "apmi_nullmodel.hpp"

class APMINullModelTest : public ::testing::Test {
protected:
  static constexpr std::size_t n_samps = 1000u;
  static constexpr std::size_t n_nulls = 10'000u;
  static constexpr uint32_t seed = 0u;

  std::filesystem::path temp_dir =
      std::filesystem::temp_directory_path() /
      ("apmi_nullmodel_test_" + std::to_string(std::rand()));

  APMINullModel model;

  void SetUp() override {
    // make the test dir root

    std::filesystem::remove_all(temp_dir); // sometimes artifacts from testing
    std::filesystem::create_directories(temp_dir);

    // ---- Set up model ----

    model = APMINullModel(n_samps, n_nulls, seed);
  }

  void TearDown() override { std::filesystem::remove_all(temp_dir); }
};

TEST_F(APMINullModelTest, Initialization) {
  auto [null_mis, ols_m, ols_b, n_samps_nm, n_nulls_nm, seed_nm] =
      model.getModel();

  // Check basic initialization
  ASSERT_EQ(n_samps_nm, n_samps);
  ASSERT_EQ(n_nulls_nm, n_nulls);
  ASSERT_EQ(seed_nm, seed);
  ASSERT_EQ(null_mis.size(), n_nulls);



TEST_F(APMINullModelTest, ModelCaching) {
  const std::string cached_blob_path = temp_dir / "apminullmodel.blob";

  model.cacheModel(cached_blob_path);

  // Retrieve the model from the cache
  APMINullModel retrieved_model =
      APMINullModel::getCachedModel(cached_blob_path);

  // Access internal state for comparison
  auto [original_null_mis, original_ols_m, original_ols_b, original_n_samps,
        original_n_nulls, original_seed] = model.getModel();
  auto [retrieved_null_mis, retrieved_ols_m, retrieved_ols_b, retrieved_n_samps,
        retrieved_n_nulls, retrieved_seed] = retrieved_model.getModel();

  // Compare the original and retrieved models
  EXPECT_EQ(original_null_mis, retrieved_null_mis);
  EXPECT_FLOAT_EQ(original_ols_m, retrieved_ols_m);
  EXPECT_FLOAT_EQ(original_ols_b, retrieved_ols_b);
  EXPECT_EQ(original_n_samps, retrieved_n_samps);
  EXPECT_EQ(original_n_nulls, retrieved_n_nulls);
  EXPECT_EQ(original_seed, retrieved_seed);
}

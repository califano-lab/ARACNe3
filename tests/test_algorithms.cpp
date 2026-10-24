#include <gtest/gtest.h>
#include "algorithms.hpp"
#include <algorithm>

TEST(AlgorithmsTest, RankIndicesTest) {
  // Test the rankIndices function

  // For example: ASSERT_EQ(f(x), y);
  ASSERT_EQ(0, 0);
}

TEST(AlgorithmsTest, CalcSCCPerfectPositive) {
    std::vector<uint16_t> x_ranked{1, 2, 3, 4, 5};
    std::vector<uint16_t> y_ranked{1, 2, 3, 4, 5};
    float expected = 1.0f; // Perfect correlation
    EXPECT_NEAR(expected, calcSCC(x_ranked, y_ranked), 0.0001f);
}

TEST(AlgorithmsTest, CalcSCCPerfectNegative) {
    std::vector<uint16_t> x_ranked{1, 2, 3, 4, 5};
    std::vector<uint16_t> y_ranked{5, 4, 3, 2, 1};
    float expected = -1.0f; // Perfect negative correlation
    EXPECT_NEAR(expected, calcSCC(x_ranked, y_ranked), 0.0001f);
}

TEST(AlgorithmsTest, CalcSCCNoCorrelation) {
    // Generate large vectors
    const size_t large_size = 10000;
    std::vector<uint16_t> x_ranked(large_size);
    std::vector<uint16_t> y_ranked(large_size);

    // Fill them with a sequence of numbers
    std::iota(x_ranked.begin(), x_ranked.end(), 0);
    std::iota(y_ranked.begin(), y_ranked.end(), 0);

    // Shuffle one of the vectors to simulate no correlation
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 eng(rd()); // Seed the generator
    std::shuffle(y_ranked.begin(), y_ranked.end(), eng);

    // We expect the correlation to be close to 0, but due to the random nature, it might not be exactly 0
    float expected = 0.0f;
    float result = calcSCC(x_ranked, y_ranked);

    // Since this is random, we allow a bit more tolerance here
    EXPECT_NEAR(expected, result, 0.1f);
}

// Regression test: concurrent calcAPMI() calls with different thresholds must
// each match their own sequential result (thresholds used to be shared statics).
TEST(AlgorithmsTest, CalcAPMIThreadSafeAcrossDifferingThresholds) {
  const int n = 200;
  std::mt19937 rng(42);
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  std::vector<float> x(n), y(n);
  for (int i = 0; i < n; ++i) {
    x[i] = dist(rng);
    y[i] = 0.7f * x[i] + 0.3f * dist(rng);
  }

  // Dependent x, y so that each threshold pair gives a different MI value
  const std::vector<std::pair<float, uint16_t>> thresholds = {
      {3.0f, 2}, {7.815f, 4}, {15.0f, 8}, {30.0f, 16}};

  std::vector<float> expected(thresholds.size());
  for (size_t t = 0; t < thresholds.size(); ++t) {
    expected[t] = calcAPMI(x, y, thresholds[t].first, thresholds[t].second);
  }

  const int reps = 200;
  const int total = reps * (int)thresholds.size();
  std::vector<float> got(total);

#pragma omp parallel for num_threads(8)
  for (int i = 0; i < total; ++i) {
    const size_t t = i % thresholds.size();
    got[i] = calcAPMI(x, y, thresholds[t].first, thresholds[t].second);
  }

  for (int i = 0; i < total; ++i) {
    const size_t t = i % thresholds.size();
    EXPECT_FLOAT_EQ(expected[t], got[i])
        << "mismatch at rep " << i << " for threshold pair " << t
        << " (q_thresh=" << thresholds[t].first
        << ", size_thresh=" << thresholds[t].second << ")";
  }
}

// lchoose
TEST(AlgorithmsTest, LchooseBasic) {
  EXPECT_NEAR(0.0, lchoose(5, 0), 1e-9); // n choose 0 should always be 1, log(1) is 0
  EXPECT_NEAR(0.0, lchoose(5, 5), 1e-9); // n choose n should always be 1, log(1) is 0
  EXPECT_NEAR(std::log(10.0), lchoose(5, 2), 1e-9); // 5 choose 2 is 10, log(10) should be the result
  EXPECT_NEAR(std::log(9.77449461715677e103), lchoose(350, 175), 1e-9);  // google'd result
}

// rightTailBinomialP 
TEST(AlgorithmsTest, RightTailBinomialPKnownValues) {
    EXPECT_NEAR(0.5, rightTailBinomialP(9, 5, 0.5), 1e-5);
    EXPECT_NEAR(0.0107421875, rightTailBinomialP(10, 9, 0.5), 1e-5);
    EXPECT_NEAR(0.00, rightTailBinomialP(350, 349, 0.01), 1e-5);
}

TEST(AlgorithmsTest, RightTailBinomialPExtremeValues) {
    EXPECT_NEAR(1.0, rightTailBinomialP(10, 0, 0.5), 1e-5); 
    EXPECT_NEAR(0.0009765625, rightTailBinomialP(10, 10, 0.5), 1e-5); 
}

// lRightTailBinomialP 
TEST(AlgorithmsTest, LogRightTailBinomialPKnownValues) {
    EXPECT_NEAR(std::log(0.5), lRightTailBinomialP(9, 5, 0.5), 1e-5);
    EXPECT_NEAR(std::log(0.0107421875), lRightTailBinomialP(10, 9, 0.5), 1e-5);
    EXPECT_NEAR(-200, lRightTailBinomialP(200, 200, 1./std::exp(1)), 1e-5);
    EXPECT_NEAR(-65534, lRightTailBinomialP(65534U, 65534U, 1./std::exp(1)), 1e-5);
}

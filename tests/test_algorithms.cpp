#include <gtest/gtest.h>
#include "algorithms.hpp"

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

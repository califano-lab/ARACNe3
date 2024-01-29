#include <gtest/gtest.h>
#include "algorithms.hpp"

TEST(AlgorithmsTest, RankIndicesTest) {
  // Test the rankIndices function

  // For example: ASSERT_EQ(f(x), y);
  ASSERT_EQ(0, 0);
}

//TEST(AlgorithmsTest, CalcSCCPerfectPositive) {
//    std::vector<uint16_t> x_ranked{1, 2, 3, 4, 5};
//    std::vector<uint16_t> y_ranked{1, 2, 3, 4, 5};
//    float expected = 1.0f; // Perfect correlation
//    EXPECT_NEAR(expected, calcSCC(x_ranked, y_ranked), 0.0001f);
//}
//
//TEST(AlgorithmsTest, CalcSCCPerfectNegative) {
//    std::vector<uint16_t> x_ranked{1, 2, 3, 4, 5};
//    std::vector<uint16_t> y_ranked{5, 4, 3, 2, 1};
//    float expected = -1.0f; // Perfect negative correlation
//    EXPECT_NEAR(expected, calcSCC(x_ranked, y_ranked), 0.0001f);
//}
//
//TEST(AlgorithmsTest, CalcSCCNoCorrelation) {
//    // Generate large vectors
//    const size_t large_size = 10000;
//    std::vector<uint16_t> x_ranked(large_size);
//    std::vector<uint16_t> y_ranked(large_size);
//
//    // Fill them with a sequence of numbers
//    std::iota(x_ranked.begin(), x_ranked.end(), 0);
//    std::iota(y_ranked.begin(), y_ranked.end(), 0);
//
//    // Shuffle one of the vectors to simulate no correlation
//    std::random_device rd;  // Obtain a random number from hardware
//    std::mt19937 eng(rd()); // Seed the generator
//    std::shuffle(y_ranked.begin(), y_ranked.end(), eng);
//
//    // We expect the correlation to be close to 0, but due to the random nature, it might not be exactly 0
//    float expected = 0.0f;
//    float result = calcSCC(x_ranked, y_ranked);
//
//    // Since this is random, we allow a bit more tolerance here
//    EXPECT_NEAR(expected, result, 0.1f);
//}

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

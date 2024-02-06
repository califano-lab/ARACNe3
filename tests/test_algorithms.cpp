#include <gtest/gtest.h>
#include "algorithms.hpp"
#include <unordered_set>

int seed = 0;

// ---- rankWithRandomTiebreak ----

TEST(RankWithRandomTiebreakTest, DistinctElements) {
    std::mt19937 rnd(seed);
    std::vector<float> vec{5.0, 1.0, 3.0, 4.0, 2.0};
    std::vector<uint32_t> expected{5, 1, 3, 4, 2};
    auto ranks = rankWithRandomTiebreak(vec, rnd);
    EXPECT_EQ(expected, ranks);
}

TEST(RankWithRandomTiebreakTest, AllElementsEqual) {
    std::mt19937 rnd(seed);
    std::vector<float> vec{1.0, 1.0, 1.0, 1.0, 1.0};
    auto ranks = rankWithRandomTiebreak(vec, rnd);
    // Check if ranks are a permutation of {1, 2, 3, 4, 5}
    std::vector<uint32_t> sortedRanks = ranks;
    std::sort(sortedRanks.begin(), sortedRanks.end());
    for (uint32_t i = 0; i < sortedRanks.size(); ++i) {
        EXPECT_EQ(i + 1, sortedRanks[i]);
    }
}

TEST(RankWithRandomTiebreakTest, SomeTiedElements) {
    std::mt19937 rnd(seed);
    std::vector<float> vec{2.0, 3.0, 2.0, 3.0, 1.0};
    auto ranks = rankWithRandomTiebreak(vec, rnd);
    // Check that ranks are unique and within the correct range
    std::unordered_set<uint32_t> uniqueRanks;
    for (auto r : ranks) {
        uniqueRanks.insert(r);
        EXPECT_GE(r, 1u);
        EXPECT_LE(r, vec.size());
    }
    EXPECT_EQ(uniqueRanks.size(), vec.size());
}

TEST(RankWithRandomTiebreakTest, LargeDataSet) {
    std::mt19937 rnd(seed);
    const size_t large_size = 100'000;
    std::vector<float> vec(large_size, 0.0);
    auto ranks = rankWithRandomTiebreak(vec, rnd);
    std::vector<uint32_t> sortedRanks = ranks;
    std::sort(sortedRanks.begin(), sortedRanks.end());
    for (uint32_t i = 0; i < large_size; ++i) {
        EXPECT_EQ(i + 1, sortedRanks[i]);
    }
}

TEST(RankWithRandomTiebreakTest, Randomness) {
    std::mt19937 rnd(seed);

    const size_t size = 100;
    std::vector<float> vec(size, 0.f); // All elements are equal
                                       //
    bool different = false;
    for (size_t i = 0; i < 10; ++i) {
        auto ranks1 = rankWithRandomTiebreak(vec, rnd);
        auto ranks2 = rankWithRandomTiebreak(vec, rnd);
        if (ranks1 != ranks2) {
            different = true;
            break;
        }
    }
    EXPECT_TRUE(different);
}

// ---- rankWithAverageTiebreak ----

TEST(RankWithAverageTiebreakTest, DistinctElements) {
    std::vector<float> vec{5.0, 1.0, 3.0, 4.0, 2.0};
    std::vector<float> expected{5, 1, 3, 4, 2};
    auto ranks = rankWithAverageTiebreak(vec);
    EXPECT_EQ(expected, ranks);
}

TEST(RankWithAverageTiebreakTest, AllElementsEqual) {
    std::vector<float> vec{1.0, 1.0, 1.0, 1.0, 1.0};
    std::vector<float> expected{3, 3, 3, 3, 3}; // Average rank for all elements
    auto ranks = rankWithAverageTiebreak(vec);
    EXPECT_EQ(expected, ranks);
}

TEST(RankWithAverageTiebreakTest, SomeTiedElements) {
    std::vector<float> vec{2.0, 3.0, 2.0, 3.0, 1.0};
    std::vector<float> expected{2.5, 4.5, 2.5, 4.5, 1}; // Average ranks for ties
    auto ranks = rankWithAverageTiebreak(vec);
    EXPECT_EQ(expected, ranks);
}

TEST(RankWithAverageTiebreakTest, MixedTies) {
    std::vector<float> vec{4.0, 4.0, 1.0, 2.0, 2.0, 3.0, 3.0};
    std::vector<float> expected{6.5, 6.5, 1.0, 2.5, 2.5, 4.5, 4.5};
    auto ranks = rankWithAverageTiebreak(vec);
    EXPECT_EQ(expected, ranks);
}

TEST(RankWithAverageTiebreakTest, LargeDataSet) {
    const size_t large_size = 100'000;
    std::vector<float> vec(large_size, 0.0);
    std::vector<float> expected(large_size, 50000.5);
    auto ranks = rankWithAverageTiebreak(vec);
    EXPECT_EQ(expected, ranks);
}

// ---- pearsonsR ----

TEST(PearsonsRTest, PerfectPositive) {
    std::vector<float> x_vec = {1, 2, 3, 4, 5};
    std::vector<float> y_vec = {1, 2, 3, 4, 5};
    EXPECT_FLOAT_EQ(pearsonsR(x_vec, y_vec), 1.0);
}

TEST(PearsonsRTest, PerfectNegative) {
    std::vector<float> x_vec = {1, 2, 3, 4, 5};
    std::vector<float> y_vec = {5, 4, 3, 2, 1};
    EXPECT_FLOAT_EQ(pearsonsR(x_vec, y_vec), -1.0);
}

TEST(PearsonsRTest, NoCorrelation) {
    std::vector<float> x_vec = {1, 2, 3, 4, 5};
    std::vector<float> y_vec = {2, 2, 2.1, 2, 2};  // not defined for constant
    EXPECT_NEAR(pearsonsR(x_vec, y_vec), 0.0, 0.001);
}

TEST(PearsonsRTest, ModeratePositive) {
    std::vector<float> x_vec = {1, 2, 3, 4, 5};
    std::vector<float> y_vec = {2, 3, 4, 4, 6};
    EXPECT_TRUE(pearsonsR(x_vec, y_vec) > 0.0 && pearsonsR(x_vec, y_vec) < 1.0);
}

TEST(PearsonsRTest, ModerateNegative) {
    std::vector<float> x_vec = {1, 2, 3, 4, 5};
    std::vector<float> y_vec = {8, 7, 5, 4, 3};
    EXPECT_TRUE(pearsonsR(x_vec, y_vec) < 0.0 && pearsonsR(x_vec, y_vec) > -1.0);
}

// ---- spearmansRho ----

TEST(SpearmansRhoTest, PerfectPositive) {
    std::mt19937 rnd(seed);

    std::vector<float> x_ranked{1, 2, 3, 4, 5};
    std::vector<float> y_ranked{1, 2, 3, 4, 5};
    float expected = 1.0f; // Perfect correlation
    EXPECT_NEAR(expected, spearmansRho(x_ranked, y_ranked), 0.0001f);
}

TEST(SpearmansRhoTest, PerfectNegative) {
    std::mt19937 rnd(seed);

    std::vector<float> x_ranked{1, 2, 3, 4, 5};
    std::vector<float> y_ranked{5, 4, 3, 2, 1};
    float expected = -1.0f; // Perfect negative correlation
    EXPECT_NEAR(expected, spearmansRho(x_ranked, y_ranked), 0.0001f);
}

TEST(SpearmansRhoTest, NoCorrelation) {
    std::mt19937 rnd(seed);

    // Generate large vectors
    const size_t large_size = 100'000;
    std::vector<float> x_ranked(large_size);
    std::vector<float> y_ranked(large_size);

    // Fill them with a sequence of numbers
    std::iota(x_ranked.begin(), x_ranked.end(), 0u);
    std::iota(y_ranked.begin(), y_ranked.end(), 0u);

    // Shuffle one of the vectors to simulate no correlation
    std::shuffle(x_ranked.begin(), x_ranked.end(), rnd);

    // Since this is random, we allow a bit more tolerance here
    EXPECT_NEAR(0.f, spearmansRho(x_ranked, y_ranked), 0.01f);
}

// ---- lchoose ----

TEST(LChooseTest, Basic) {
  EXPECT_NEAR(0.0, lChoose(5, 0), 1e-9); // n choose 0 should always be 1, log(1) is 0
  EXPECT_NEAR(0.0, lChoose(5, 5), 1e-9); // n choose n should always be 1, log(1) is 0
  EXPECT_NEAR(std::log(10.0), lChoose(5, 2), 1e-9); // 5 choose 2 is 10, log(10) should be the result
  EXPECT_NEAR(std::log(9.77449461715677e103), lChoose(350, 175), 1e-9);  // google'd result
}

// ---- rightTailBinomialP ----

TEST(RightTailBinomialPTest, KnownValues) {
    EXPECT_NEAR(0.5, rightTailBinomialP(9, 5, 0.5), 1e-5);
    EXPECT_NEAR(0.0107421875, rightTailBinomialP(10, 9, 0.5), 1e-5);
    EXPECT_NEAR(0.00, rightTailBinomialP(350, 349, 0.01), 1e-5);
}

TEST(RightTailBinomialPTest, ExtremeValues) {
    EXPECT_NEAR(1.0, rightTailBinomialP(10, 0, 0.5), 1e-5); 
    EXPECT_NEAR(0.0009765625, rightTailBinomialP(10, 10, 0.5), 1e-5); 
}

// ---- lRightTailBinomialP ----

TEST(LRightTailBinomialPTest, KnownValues) {
    EXPECT_NEAR(std::log(0.5), lRightTailBinomialP(9, 5, 0.5), 1e-5);
    EXPECT_NEAR(std::log(0.0107421875), lRightTailBinomialP(10, 9, 0.5), 1e-5);
    EXPECT_NEAR(-200, lRightTailBinomialP(200, 200, 1./std::exp(1)), 1e-5);
    EXPECT_NEAR(-65534, lRightTailBinomialP(65534U, 65534U, 1./std::exp(1)), 1e-5);
}

// ---- OLS ----
TEST(OLSTest, KnownPositive) {
    std::vector<float> x = {1, 2, 3, 4, 5};
    std::vector<float> y = {2, 4, 6, 8, 10}; // Perfect linear relationship y = 2x
    auto [m, b] = OLS(x, y);
    EXPECT_NEAR(2.0, m, 1e-5); // Expect slope of 2
    EXPECT_NEAR(0.0, b, 1e-5); // Expect intercept of 0
}

TEST(OLSTest, KnownNegative) {
    std::vector<float> x = {1, 2, 3, 4, 5};
    std::vector<float> y = {-1, -2, -3, -4, -5};  // Perfect negative cor
    auto [m, b] = OLS(x, y);
    EXPECT_NEAR(-1.0, m, 1e-5); // Expect slope of 1
    EXPECT_NEAR(0.0, b, 1e-5); // Expect intercept of 0
}

// ---- copulaTransform ----
TEST(CopulaTransformTest, KnownValues) {
    std::mt19937 rnd(seed); // Fixed seed for reproducibility
    std::vector<float> data = {5, 2, 3, 4, 1};
    auto transformed = copulaTransform(data, rnd);
    std::vector<float> expected {5./5, 2./5, 3./5, 4./5, 1./5};
    EXPECT_EQ(5, transformed.size()); // Expect same size
    EXPECT_EQ(expected, transformed);

    data = {2*5, 2*2, 2*3, 2*4, 2*1};
    transformed = copulaTransform(data, rnd);
    EXPECT_EQ(expected, transformed);
}

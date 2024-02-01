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
    const size_t large_size = 1'000'000;
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
    const size_t large_size = 1'000'000;
    std::vector<float> vec(large_size, 0.0);
    std::vector<float> expected(large_size, 500000.5);
    auto ranks = rankWithAverageTiebreak(vec);
    EXPECT_EQ(expected, ranks);
}

// ---- spearmansRho ----

TEST(AlgorithmsTest, SpearmansRhoPerfectPositive) {
    std::mt19937 rnd(seed);

    std::vector<float> x_ranked{1, 2, 3, 4, 5};
    std::vector<float> y_ranked{1, 2, 3, 4, 5};
    float expected = 1.0f; // Perfect correlation
    EXPECT_NEAR(expected, spearmansRho(x_ranked, y_ranked), 0.0001f);
}

TEST(AlgorithmsTest, SpearmansRhoPerfectNegative) {
    std::mt19937 rnd(seed);

    std::vector<float> x_ranked{1, 2, 3, 4, 5};
    std::vector<float> y_ranked{5, 4, 3, 2, 1};
    float expected = -1.0f; // Perfect negative correlation
    EXPECT_NEAR(expected, spearmansRho(x_ranked, y_ranked), 0.0001f);
}

TEST(AlgorithmsTest, SpearmansRhoNoCorrelation) {
    std::mt19937 rnd(seed);

    // Generate large vectors
    const size_t large_size = 5'000'000;
    std::vector<float> x_ranked(large_size);
    std::vector<float> y_ranked(large_size);

    // Fill them with a sequence of numbers
    std::iota(x_ranked.begin(), x_ranked.end(), 0u);
    std::iota(y_ranked.begin(), y_ranked.end(), 0u);

    // Shuffle one of the vectors to simulate no correlation
    std::shuffle(x_ranked.begin(), x_ranked.end(), rnd);

    // Since this is random, we allow a bit more tolerance here
    EXPECT_NEAR(0.f, spearmansRho(x_ranked, y_ranked), 0.001f);
}

// ---- lchoose ----

TEST(AlgorithmsTest, LchooseBasic) {
  EXPECT_NEAR(0.0, lchoose(5, 0), 1e-9); // n choose 0 should always be 1, log(1) is 0
  EXPECT_NEAR(0.0, lchoose(5, 5), 1e-9); // n choose n should always be 1, log(1) is 0
  EXPECT_NEAR(std::log(10.0), lchoose(5, 2), 1e-9); // 5 choose 2 is 10, log(10) should be the result
  EXPECT_NEAR(std::log(9.77449461715677e103), lchoose(350, 175), 1e-9);  // google'd result
}

// ---- rightTailBinomialP ----

TEST(AlgorithmsTest, RightTailBinomialPKnownValues) {
    EXPECT_NEAR(0.5, rightTailBinomialP(9, 5, 0.5), 1e-5);
    EXPECT_NEAR(0.0107421875, rightTailBinomialP(10, 9, 0.5), 1e-5);
    EXPECT_NEAR(0.00, rightTailBinomialP(350, 349, 0.01), 1e-5);
}

TEST(AlgorithmsTest, RightTailBinomialPExtremeValues) {
    EXPECT_NEAR(1.0, rightTailBinomialP(10, 0, 0.5), 1e-5); 
    EXPECT_NEAR(0.0009765625, rightTailBinomialP(10, 10, 0.5), 1e-5); 
}

// ---- lRightTailBinomialP ----

TEST(AlgorithmsTest, LogRightTailBinomialPKnownValues) {
    EXPECT_NEAR(std::log(0.5), lRightTailBinomialP(9, 5, 0.5), 1e-5);
    EXPECT_NEAR(std::log(0.0107421875), lRightTailBinomialP(10, 9, 0.5), 1e-5);
    EXPECT_NEAR(-200, lRightTailBinomialP(200, 200, 1./std::exp(1)), 1e-5);
    EXPECT_NEAR(-65534, lRightTailBinomialP(65534U, 65534U, 1./std::exp(1)), 1e-5);
}

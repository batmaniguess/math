#include "../../src/goldbach/goldbach.hpp"
#include <gtest/gtest.h>

namespace {

class GoldBachTest : public ::testing::Test {
protected:
  void SetUp() override {}
};

TEST_F(GoldBachTest, HandlesValidInput) {
  auto result = math::GoldbachConjecture::verify(4);
  EXPECT_TRUE(result.verified);
  EXPECT_TRUE(result.pair.has_value());
  EXPECT_EQ(result.pair->first + result.pair->second, 4);
}

TEST_F(GoldBachTest, HandlesInvalidInput) {
  auto results = math::GoldbachConjecture::verifyRange(4, 100);
  for (const auto &result : results) {
    EXPECT_TRUE(result.verified);
    EXPECT_TRUE(result.pair.has_value());
  }
}

TEST_F(GoldBachTest, HandleLargeNumbers) {
  auto result = math::GoldbachConjecture::verify(10000);
  EXPECT_TRUE(result.verified);
  EXPECT_TRUE(result.pair.has_value());
}

}; // namespace

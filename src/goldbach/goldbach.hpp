#pragma once
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace math {

class GoldbachConjecture {
public:
  struct PrimePair {
    int32_t first;
    int32_t second;
  };

  struct Result {
    bool verified;
    std::optional<PrimePair> pair;
    std::string message;
  };

  static Result verify(int32_t number);
  static std::vector<Result> verifyRange(int32_t start, int32_t end);

  static constexpr int32_t MIN_VALUE = 4;
  static constexpr int32_t MAX_VALUE = INT32_MAX;

private:
  static bool isValidInput(int32_t number);
};
} // namespace math

#pragma once
#include <bits/stdc++.h>
#include <cstdint>
#include <vector>

namespace math::utils {

class PrimeUtils {
public:
  static bool isPrime(int32_t n);
  static std::vector<int32_t> generatePrimeInRange(int32_t start, int32_t end);

private:
  static std::vector<bool> createSieve(int32_t limit);
};

} // namespace math::utils

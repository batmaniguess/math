#include "prime_utils.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace math::utils {

bool PrimeUtils::isPrime(int32_t n) {
  if (n < 2)
    return false;
  if (n == 2)
    return true;
  if (n % 2 == 0)
    return false;

  int32_t sqrtN = static_cast<int32_t>(std::sqrt(n));
  for (int32_t i = 3; i <= sqrtN; i += 2)
    if (n % i == 0)
      return false;

  return true;
}

std::vector<int32_t> PrimeUtils::generatePrimeInRange(int32_t start,
                                                      int32_t end) {
  std::vector<int32_t> primes;
  start = std::max(start, 2);

  auto sieve = createSieve(end);

  for (int32_t i = start; i <= end; ++i) {
    if (sieve[i])
      primes.push_back(i);
  }
  return primes;
}

std::vector<bool> PrimeUtils::createSieve(int32_t limit) {
  std::vector<bool> sieve(limit + 1, true);
  if (limit >= 0)
    sieve[0] = false;
  if (limit >= 1)
    sieve[1] = false;

  int32_t sqrtLimit = static_cast<int32_t>(std::sqrt(limit));
  for (int32_t i = 2; i < sqrtLimit; ++i) {
    if (sieve[i]) {
      for (int32_t j = i * i; j <= limit; j += i)
        sieve[j] = false;
    }
  }
  return sieve;
}

} // namespace math::utils

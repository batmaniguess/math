#include "goldbach.hpp"
#include "../utils/prime_utils.hpp"
#include <optional>
#include <sstream>
#include <vector>

namespace math {
bool GoldbachConjecture::isValidInput(int32_t number) {
  return number >= MIN_VALUE && number <= MAX_VALUE && number % 2 == 0;
}

GoldbachConjecture::Result GoldbachConjecture::verify(int32_t number) {
  if (!isValidInput(number)) {
    std::stringstream ss;
    ss << "Invalid input: " << number
       << ". Number must be even and >= " << MIN_VALUE
       << " and <= " << MAX_VALUE;
    return Result{false, std::nullopt, ss.str()};
  }
  auto primes = utils::PrimeUtils::generatePrimeInRange(2, number / 2);

  for (int32_t prime1 : primes) {
    int32_t prime2 = number - prime1;
    if (utils::PrimeUtils::isPrime(prime2)) {
      return Result{true, PrimePair{prime1, prime2},
                    "Successfully found prime pair"};
    }
  }
  std::stringstream ss;
  ss << "No prime pair found that sums to " << number;
  return Result{false, std::nullopt, ss.str()};
}

std::vector<GoldbachConjecture::Result>
GoldbachConjecture::verifyRange(int32_t start, int32_t end) {
  std::vector<Result> results;

  if (start % 2 != 0)
    start++;

  for (int32_t number = start; number <= end; number += 2)
    results.push_back(verify(number));

  return results;
}

} // namespace math

#include "goldbach/goldbach.hpp"
#include <iostream>

int main() {
  std::cout << "Goldbach Conjecture Verifier\n";
  std::cout << "==================================";

  while (true) {
    std::cout << "Enter an even number to verify (0 to exit): ";
    int32_t number;
    std::cin >> number;

    if (number == 0)
      break;

    auto result = math::GoldbachConjecture::verify(number);

    if (result.verified) {
      std::cout << "\nSuccess! Found prime pair for " << number << ":\n";
      std::cout << number << " = " << result.pair->first << " + "
                << result.pair->second << "\n";

    } else {
      std::cout << "\nVerification failed: " << result.message << "\n";
    }
    std::cout << "\n";
  }
  return 0;
}

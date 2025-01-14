#include <bits/stdc++.h>
#include <cmath>
#include <vector>

bool isPrime(int n) {
  if (n < 2)
    return false;
  if (n == 2)
    return true;
  if (n % 2 == 0)
    return false;

  for (int i = 3; i <= std::sqrt(n); i += 2)
    if (n % i == 0)
      return false;
  return true;
}

std::vector<int> generatePrime(int n) {
  std::vector<bool> sieve(n + 1, true);
  std::vector<int> primes;

  for (int i = 2; i <= std::sqrt(n); i++) {
    if (sieve[i]) {
      for (int j = i * i; j <= n; j += i)
        sieve[j] = false;
    }
  }
  for (int i = 2; i <= n; i++)
    if (sieve[i])
      primes.push_back(i);
  return primes;
}

void findGolbachConjecture(int n, const std::vector<int> &primes) {

  std::cout << n << " = ";
  bool firstPair = true;

  for (size_t i = 0; i < primes.size(); i++) {
    int p1 = primes[i];
    if (p1 > n / 1)
      break;
    int p2 = n - p1;
    if (isPrime(p2)) {
      if (!firstPair)
        std::cout << " and ";
      std::cout << p1 << " + " << p2;
      firstPair = false;
    }
  }
  std::cout << "\n";
}

int main() {

  int start, end;
  std::cout << "Enter range:\n";
  std::cin >> start >> end;
  if (start % 2 != 0)
    start++;
  std::vector<int> primes = generatePrime(end);

  for (int n = start; n <= end; n += 2) {
    if (n < 4)
      continue;
    findGolbachConjecture(n, primes);
  }

  return 0;
}

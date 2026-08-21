#include <cstdint>
#include <iostream>

namespace {

using Integer = std::uint64_t;

Integer power(Integer base, unsigned exponent) {
  Integer answer = 1;
  while (exponent-- != 0) {
    answer *= base;
  }
  return answer;
}

}  // namespace

int main() {
  // The GAP program finds three distinct projective regular orbits.  Each
  // orbit for the maximal scalar group has the order shown below.
  constexpr Integer number_of_regular_orbits = 3;
  constexpr Integer suz_order = 448345497600ULL;
  constexpr Integer nonzero_scalars = 15;

  const Integer maximal_group_order = nonzero_scalars * suz_order;
  const Integer regular_lower_bound =
      number_of_regular_orbits * maximal_group_order;
  const Integer density_threshold = power(16, 11) - 1;

  std::cout << "3.Suz on GF(16)^12\n"
            << "distinct regular projective orbits: "
            << number_of_regular_orbits << '\n'
            << "lower bound for regular vectors:    "
            << regular_lower_bound << '\n'
            << "density threshold:                  "
            << density_threshold << '\n';

  if (regular_lower_bound <= density_threshold) {
    std::cerr << "The density inequality was not obtained.\n";
    return 1;
  }

  std::cout
      << "Conclusion: every vector is a sum of two regular vectors.\n";
}

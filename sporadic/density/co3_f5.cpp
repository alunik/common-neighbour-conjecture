#include <array>
#include <cstdint>
#include <iostream>
#include <string_view>

namespace {

using Integer = std::uint64_t;

// For each prime-order class of Co3, the GAP calculation supplies the class
// size and the dimensions of its eigenspaces for the four nonzero elements
// of GF(5).  A nonregular vector must occur in one of these eigenspaces.
struct ClassData {
  std::string_view name;
  Integer size;
  std::array<unsigned, 4> eigenspace_dimensions;
};

constexpr std::array<ClassData, 12> classes{{
    {"2A", 170775, {15, 0, 0, 8}},
    {"2B", 2608200, {11, 0, 0, 12}},
    {"3A", 1416800, {5, 0, 0, 0}},
    {"3B", 17001600, {11, 0, 0, 0}},
    {"3C", 109296000, {7, 0, 0, 0}},
    {"5A", 330511104, {5, 0, 0, 0}},
    {"5B", 1652555520, {7, 0, 0, 0}},
    {"7A", 11803968000ULL, {5, 0, 0, 0}},
    {"11A", 22534848000ULL, {3, 0, 0, 0}},
    {"11B", 22534848000ULL, {3, 0, 0, 0}},
    {"23A", 21555072000ULL, {1, 0, 0, 0}},
    {"23B", 21555072000ULL, {1, 0, 0, 0}},
}};

Integer power(Integer base, unsigned exponent) {
  Integer answer = 1;
  while (exponent-- != 0) {
    answer *= base;
  }
  return answer;
}

}  // namespace

int main() {
  Integer nonregular_upper_bound = 0;

  for (const ClassData& conjugacy_class : classes) {
    Integer vectors_in_eigenspaces = 0;
    for (unsigned dimension : conjugacy_class.eigenspace_dimensions) {
      vectors_in_eigenspaces += power(5, dimension) - 1;
    }
    nonregular_upper_bound +=
        conjugacy_class.size * vectors_in_eigenspaces;
  }

  const Integer nonzero_vectors = power(5, 23) - 1;
  const Integer regular_lower_bound =
      nonzero_vectors - nonregular_upper_bound;
  const Integer density_threshold = power(5, 22) - 1;

  std::cout << "Co3 on GF(5)^23\n"
            << "upper bound for nonregular vectors: "
            << nonregular_upper_bound << '\n'
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

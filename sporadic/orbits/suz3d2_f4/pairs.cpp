// Enumerate stabiliser orbits on the second binary coordinate. Together with
// the five first-coordinate orbits this gives all orbits on GF(4)^24.

#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr unsigned kDim = 24;
constexpr std::uint32_t kSpaceSize = std::uint32_t{1} << kDim;
constexpr std::uint64_t kPairSpaceSize =
    std::uint64_t{kSpaceSize} * kSpaceSize;

struct Matrix {
  std::uint32_t row[kDim]{};
};

struct Case {
  std::uint32_t representative = 0;
  std::uint64_t first_orbit_length = 0;
  std::uint64_t stabilizer_order = 0;
  std::vector<Matrix> generators;
};

std::uint32_t Act(std::uint32_t vector, const Matrix& matrix) {
  std::uint32_t image = 0;
  while (vector != 0) {
    const unsigned position =
        static_cast<unsigned>(__builtin_ctz(vector));
    image ^= matrix.row[position];
    vector &= vector - 1;
  }
  return image;
}

std::vector<Case> ReadCases(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open " + path);
  }
  std::vector<Case> cases;
  std::string token;
  while (input >> token) {
    if (token != "CASE") {
      continue;
    }
    Case item;
    std::size_t generator_count = 0;
    std::uint64_t support_size = 0;
    if (!(input >> item.representative >> item.first_orbit_length >>
          item.stabilizer_order >> generator_count >> support_size)) {
      throw std::runtime_error("invalid CASE header");
    }
    if (generator_count == 0) {
      throw std::runtime_error("empty stabilizer generating set");
    }
    for (std::size_t index = 0; index < generator_count; ++index) {
      if (!(input >> token) || token != "MATRIX") {
        throw std::runtime_error("missing MATRIX record");
      }
      Matrix matrix;
      for (unsigned position = 0; position < kDim; ++position) {
        std::uint64_t row = 0;
        if (!(input >> row) || row >= kSpaceSize) {
          throw std::runtime_error("invalid packed matrix row");
        }
        matrix.row[position] = static_cast<std::uint32_t>(row);
      }
      item.generators.push_back(matrix);
    }
    if (!(input >> token) || token != "ENDCASE") {
      throw std::runtime_error("missing ENDCASE record");
    }
    cases.push_back(std::move(item));
  }
  if (cases.size() != 5) {
    throw std::runtime_error("expected exactly five first-vector cases");
  }
  return cases;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: pair_census STABILIZERS_DATA ORBITS_TSV\n";
    return 2;
  }
  try {
    const std::vector<Case> cases = ReadCases(argv[1]);
    std::ofstream details(argv[2]);
    if (!details) {
      throw std::runtime_error("cannot create detailed orbit table");
    }
    details << "orbit_index\tfirst_rep\tsecond_rep\tfirst_orbit_length"
               "\tsecond_orbit_length\tpair_orbit_length"
               "\tpair_stabilizer_order\tregular\n";

    std::vector<std::uint8_t> seen(kSpaceSize);
    std::vector<std::uint32_t> queue(kSpaceSize);
    std::uint64_t total_pair_orbits = 0;
    std::uint64_t total_regular_pair_orbits = 0;
    std::uint64_t weighted_coverage = 0;

    std::cout << "SUZ3D2_F4_ORDERED_BINARY_PAIR_CENSUS_V1\n";
    std::cout << "binary_space_size " << kSpaceSize << "\n";
    for (const Case& item : cases) {
      std::fill(seen.begin(), seen.end(), std::uint8_t{0});
      std::uint64_t covered = 0;
      std::uint64_t orbit_count = 0;
      std::uint64_t maximum = 0;
      std::map<std::uint64_t, std::uint64_t> distribution;
      std::vector<std::uint32_t> regular_second_representatives;

      for (std::uint32_t seed = 0; seed < kSpaceSize; ++seed) {
        if (seen[seed] != 0) {
          continue;
        }
        std::size_t head = 0;
        std::size_t tail = 0;
        seen[seed] = 1;
        queue[tail++] = seed;
        while (head < tail) {
          const std::uint32_t point = queue[head++];
          for (const Matrix& generator : item.generators) {
            const std::uint32_t image = Act(point, generator);
            if (seen[image] == 0) {
              seen[image] = 1;
              queue[tail++] = image;
            }
          }
        }
        const std::uint64_t second_orbit_length = tail;
        if (item.stabilizer_order % second_orbit_length != 0) {
          throw std::runtime_error(
              "second orbit length does not divide stabilizer order");
        }
        const std::uint64_t pair_orbit_length =
            item.first_orbit_length * second_orbit_length;
        const std::uint64_t pair_stabilizer_order =
            item.stabilizer_order / second_orbit_length;
        const bool regular = pair_stabilizer_order == 1;
        ++orbit_count;
        ++total_pair_orbits;
        covered += second_orbit_length;
        weighted_coverage += pair_orbit_length;
        maximum = std::max(maximum, second_orbit_length);
        ++distribution[second_orbit_length];
        if (regular) {
          regular_second_representatives.push_back(seed);
          ++total_regular_pair_orbits;
        }
        details << total_pair_orbits << '\t' << item.representative << '\t'
                << seed << '\t' << item.first_orbit_length << '\t'
                << second_orbit_length << '\t' << pair_orbit_length << '\t'
                << pair_stabilizer_order << '\t'
                << (regular ? 1 : 0) << '\n';
      }

      std::cout << "case rep " << item.representative
                << " first_orbit " << item.first_orbit_length
                << " stabilizer_order " << item.stabilizer_order
                << " generators " << item.generators.size() << "\n";
      std::cout << "orbit_count " << orbit_count
                << " covered " << covered
                << " maximum " << maximum
                << " minimum_pair_stabilizer "
                << item.stabilizer_order / maximum << "\n";
      std::cout << "distribution length multiplicity\n";
      for (const auto& [length, multiplicity] : distribution) {
        std::cout << length << " " << multiplicity << "\n";
      }
      std::cout << "regular_second_representatives";
      for (const std::uint32_t representative :
           regular_second_representatives) {
        std::cout << " " << representative;
      }
      std::cout << "\n";
      if (covered != kSpaceSize) {
        throw std::runtime_error(
            "stabilizer orbits do not cover the second-vector space");
      }
      std::cout << "case_result "
                << (regular_second_representatives.empty()
                        ? "NO_REGULAR_SECOND_VECTOR"
                        : "REGULAR_SECOND_VECTORS_EXIST")
                << "\n";
    }
    if (weighted_coverage != kPairSpaceSize) {
      throw std::runtime_error(
          "weighted pair orbits do not cover the ordered-pair space");
    }
    std::cout << "all_five_first_vector_orbits_covered true\n";
    std::cout << "weighted_pair_space_covered " << weighted_coverage << "\n";
    std::cout << "total_pair_orbits " << total_pair_orbits << "\n";
    std::cout << "total_regular_pair_orbits "
              << total_regular_pair_orbits << "\n";
    std::cout << "Pair-orbit enumeration finished.\n";
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
  return 0;
}

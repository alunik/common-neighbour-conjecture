// For each first-coordinate Co3-orbit, enumerate the orbits of its stabiliser
// on the second binary coordinate. These are the Co3-orbits on GF(4)^22 after
// writing a vector as a + omega*b.

#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr unsigned kDim = 22;
constexpr uint32_t kSpaceSize = uint32_t{1} << kDim;

struct Matrix {
  uint32_t row[kDim]{};
};

struct Case {
  uint32_t representative = 0;
  uint64_t first_orbit_length = 0;
  uint64_t stabilizer_order = 0;
  std::vector<Matrix> generators;
};

uint32_t act(uint32_t v, const Matrix& g) {
  uint32_t image = 0;
  while (v != 0) {
    const unsigned i = static_cast<unsigned>(__builtin_ctz(v));
    image ^= g.row[i];
    v &= v - 1;
  }
  return image;
}

std::vector<Case> read_cases(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::vector<Case> cases;
  std::string token;
  while (in >> token) {
    if (token != "CASE") continue;
    Case c;
    size_t generator_count = 0;
    uint64_t support_size = 0;
    if (!(in >> c.representative >> c.first_orbit_length >>
          c.stabilizer_order >> generator_count >> support_size)) {
      throw std::runtime_error("invalid CASE header");
    }
    if (generator_count == 0) {
      throw std::runtime_error("empty stabilizer generating set");
    }
    for (size_t j = 0; j < generator_count; ++j) {
      if (!(in >> token) || token != "MATRIX") {
        throw std::runtime_error("missing MATRIX record");
      }
      Matrix g;
      for (unsigned i = 0; i < kDim; ++i) {
        uint64_t x = 0;
        if (!(in >> x) || x >= kSpaceSize) {
          throw std::runtime_error("invalid packed matrix row");
        }
        g.row[i] = static_cast<uint32_t>(x);
      }
      c.generators.push_back(g);
    }
    if (!(in >> token) || token != "ENDCASE") {
      throw std::runtime_error("missing ENDCASE record");
    }
    cases.push_back(std::move(c));
  }
  if (cases.size() != 5) {
    throw std::runtime_error("expected exactly five first-vector cases");
  }
  return cases;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: pair_census STABILIZERS_DATA\n";
    return 2;
  }
  try {
    const auto cases = read_cases(argv[1]);
    std::vector<uint8_t> seen(kSpaceSize);
    std::vector<uint32_t> queue(kSpaceSize);

    std::cout << "CO3_F4_ORDERED_BINARY_PAIR_CENSUS_V1\n";
    std::cout << "space_size " << kSpaceSize << "\n";
    uint64_t total_pair_orbits = 0;
    uint64_t total_regular_pair_orbits = 0;

    for (const Case& c : cases) {
      std::fill(seen.begin(), seen.end(), uint8_t{0});
      uint64_t covered = 0;
      uint64_t orbit_count = 0;
      uint64_t maximum = 0;
      std::map<uint64_t, uint64_t> distribution;
      std::vector<uint32_t> regular_second_representatives;

      for (uint32_t seed = 0; seed < kSpaceSize; ++seed) {
        if (seen[seed]) continue;
        size_t head = 0;
        size_t tail = 0;
        seen[seed] = 1;
        queue[tail++] = seed;
        while (head < tail) {
          const uint32_t x = queue[head++];
          for (const Matrix& g : c.generators) {
            const uint32_t y = act(x, g);
            if (!seen[y]) {
              seen[y] = 1;
              queue[tail++] = y;
            }
          }
        }
        const uint64_t length = tail;
        if (c.stabilizer_order % length != 0) {
          throw std::runtime_error(
              "orbit length does not divide certified stabilizer order");
        }
        ++orbit_count;
        covered += length;
        if (length > maximum) maximum = length;
        if (length == c.stabilizer_order) {
          regular_second_representatives.push_back(seed);
        }
        ++distribution[length];
      }

      std::cout << "case rep " << c.representative
                << " first_orbit " << c.first_orbit_length
                << " stabilizer_order " << c.stabilizer_order
                << " generators " << c.generators.size() << "\n";
      std::cout << "orbit_count " << orbit_count
                << " covered " << covered
                << " maximum " << maximum
                << " minimum_pair_stabilizer "
                << c.stabilizer_order / maximum << "\n";
      std::cout << "distribution length multiplicity\n";
      for (const auto& [length, multiplicity] : distribution) {
        std::cout << length << " " << multiplicity << "\n";
      }
      std::cout << "regular_second_representatives";
      for (const uint32_t representative : regular_second_representatives) {
        std::cout << " " << representative;
      }
      std::cout << "\n";
      if (covered != kSpaceSize) {
        throw std::runtime_error("stabilizer orbits do not cover GF(2)^22");
      }
      std::cout << "case_result "
                << (regular_second_representatives.empty()
                        ? "NO_REGULAR_SECOND_VECTOR"
                        : "REGULAR_SECOND_VECTORS_EXIST")
                << "\n";
      total_pair_orbits += orbit_count;
      total_regular_pair_orbits += regular_second_representatives.size();
    }
    std::cout << "all_five_first_vector_orbits_covered true\n";
    std::cout << "total_pair_orbits " << total_pair_orbits << "\n";
    std::cout << "total_regular_pair_orbits "
              << total_regular_pair_orbits << "\n";
    std::cout << "Pair-orbit enumeration finished.\n";
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

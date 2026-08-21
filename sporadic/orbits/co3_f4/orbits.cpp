// Enumerate Co3-orbits on the binary 22-dimensional module. A vector is
// stored as a 22-bit integer, so applying a matrix uses shifts and XORs.

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr unsigned kDim = 22;
constexpr uint32_t kSpaceSize = uint32_t{1} << kDim;

struct Matrix {
  uint32_t row[kDim]{};
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

std::vector<Matrix> read_matrices(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::vector<Matrix> matrices;
  std::string token;
  while (in >> token) {
    if (token != "MATRIX") continue;
    Matrix g;
    for (unsigned i = 0; i < kDim; ++i) {
      uint64_t x = 0;
      if (!(in >> x) || x >= kSpaceSize) {
        throw std::runtime_error("invalid packed matrix row");
      }
      g.row[i] = static_cast<uint32_t>(x);
    }
    matrices.push_back(g);
  }
  if (matrices.size() != 2) {
    throw std::runtime_error("expected exactly two matrices");
  }
  return matrices;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: orbit_census INITIAL_DATA\n";
    return 2;
  }
  try {
    const auto generators = read_matrices(argv[1]);
    std::vector<uint8_t> seen(kSpaceSize, 0);
    std::vector<uint32_t> queue(kSpaceSize);
    uint64_t covered = 0;
    uint64_t orbit_number = 0;

    std::cout << "CO3_F4_BINARY_VECTOR_ORBITS_V1\n";
    std::cout << "space_size " << kSpaceSize << "\n";
    std::cout << "orbit rep length\n";
    for (uint32_t seed = 0; seed < kSpaceSize; ++seed) {
      if (seen[seed]) continue;
      size_t head = 0;
      size_t tail = 0;
      seen[seed] = 1;
      queue[tail++] = seed;
      while (head < tail) {
        const uint32_t x = queue[head++];
        for (const Matrix& g : generators) {
          const uint32_t y = act(x, g);
          if (!seen[y]) {
            seen[y] = 1;
            queue[tail++] = y;
          }
        }
      }
      ++orbit_number;
      covered += tail;
      std::cout << orbit_number << " " << seed << " " << tail << "\n";
    }
    std::cout << "orbit_count " << orbit_number << "\n";
    std::cout << "covered " << covered << "\n";
    if (covered != kSpaceSize) {
      throw std::runtime_error("orbit lengths do not cover the space");
    }
    std::cout << "Orbit enumeration finished.\n";
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

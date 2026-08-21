// Enumerate 3.Suz.2-orbits on the binary 24-dimensional module. Packed
// 24-bit vectors make the full traversal small enough for a direct BFS.

#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr unsigned kDim = 24;
constexpr std::uint32_t kSpaceSize = std::uint32_t{1} << kDim;

struct Matrix {
  std::uint32_t row[kDim]{};
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

std::vector<Matrix> ReadMatrices(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open " + path);
  }
  std::vector<Matrix> matrices;
  std::string token;
  while (input >> token) {
    if (token != "MATRIX") {
      continue;
    }
    Matrix matrix;
    for (unsigned position = 0; position < kDim; ++position) {
      std::uint64_t row = 0;
      if (!(input >> row) || row >= kSpaceSize) {
        throw std::runtime_error("invalid packed matrix row");
      }
      matrix.row[position] = static_cast<std::uint32_t>(row);
    }
    matrices.push_back(matrix);
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
    const std::vector<Matrix> generators = ReadMatrices(argv[1]);
    std::vector<std::uint8_t> seen(kSpaceSize, 0);
    std::vector<std::uint32_t> queue(kSpaceSize);
    std::uint64_t covered = 0;
    std::uint64_t orbit_number = 0;

    std::cout << "SUZ3D2_F4_BINARY_VECTOR_ORBITS_V1\n";
    std::cout << "space_size " << kSpaceSize << "\n";
    std::cout << "orbit rep length\n";
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
        for (const Matrix& generator : generators) {
          const std::uint32_t image = Act(point, generator);
          if (seen[image] == 0) {
            seen[image] = 1;
            queue[tail++] = image;
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
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
  return 0;
}

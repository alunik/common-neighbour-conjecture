// Exhaustive vector-orbit BFS for an outer involution fixed space of
// 3.Fi22.2 <= GL(54,2).
//
// The GAP input contains the exact C_H(x)-action on a 27-dimensional
// fixed-space basis, as well as each basis row in ATLAS GF(2)^54
// coordinates and in the reference GF(4)^27 coordinates for the index-two
// core.  All 2^27 vectors, including zero, are partitioned into orbits.
//
// Usage:
//   outer_orbit_bfs INPUT OUTPUT THREADS

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

constexpr unsigned kDimension = 27;
constexpr unsigned kAmbientDimension = 54;
constexpr std::uint32_t kVectorCount = std::uint32_t{1} << kDimension;
constexpr std::uint32_t kChunkMask = (std::uint32_t{1} << 9) - 1;
constexpr std::size_t kVisitedWords = kVectorCount / 64;

[[noreturn]] void Fail(const std::string& message) {
  throw std::runtime_error(message);
}

std::uint64_t ParseUnsigned(
    const std::string& text, const std::string& label) {
  if (text.empty()) Fail("empty " + label);
  std::size_t consumed = 0;
  const unsigned long long value =
      std::stoull(text, &consumed, 10);
  if (consumed != text.size()) Fail("invalid " + label);
  return static_cast<std::uint64_t>(value);
}

std::vector<std::string> Split(const std::string& line) {
  std::istringstream input(line);
  std::vector<std::string> fields;
  std::string field;
  while (input >> field) fields.push_back(field);
  return fields;
}

std::uint64_t ParseBinary(
    const std::string& digits, unsigned expected,
    const std::string& label) {
  if (digits.size() != expected) {
    Fail(label + " has the wrong length");
  }
  std::uint64_t result = 0;
  for (unsigned position = 0; position < expected; ++position) {
    if (digits[position] != '0' && digits[position] != '1') {
      Fail(label + " is not binary");
    }
    if (digits[position] == '1') {
      result |= std::uint64_t{1} << position;
    }
  }
  return result;
}

std::array<std::uint8_t, kDimension> ParseF4(
    const std::string& digits) {
  if (digits.size() != kDimension) {
    Fail("GF(4) basis row has the wrong length");
  }
  std::array<std::uint8_t, kDimension> result{};
  for (unsigned position = 0; position < kDimension; ++position) {
    if (digits[position] < '0' || digits[position] > '3') {
      Fail("GF(4) basis row has an invalid digit");
    }
    result[position] =
        static_cast<std::uint8_t>(digits[position] - '0');
  }
  return result;
}

unsigned BinaryRank(std::vector<std::uint64_t> rows, unsigned width) {
  unsigned rank = 0;
  for (unsigned column = 0; column < width; ++column) {
    unsigned pivot = rank;
    while (pivot < rows.size() &&
           ((rows[pivot] >> column) & 1U) == 0) {
      ++pivot;
    }
    if (pivot == rows.size()) continue;
    std::swap(rows[rank], rows[pivot]);
    for (unsigned row = 0; row < rows.size(); ++row) {
      if (row != rank && ((rows[row] >> column) & 1U)) {
        rows[row] ^= rows[rank];
      }
    }
    ++rank;
  }
  return rank;
}

using Matrix = std::array<std::uint32_t, kDimension>;

Matrix Inverse(const Matrix& matrix) {
  std::array<std::uint64_t, kDimension> augmented{};
  for (unsigned row = 0; row < kDimension; ++row) {
    augmented[row] = matrix[row] |
        (std::uint64_t{1} << (kDimension + row));
  }
  for (unsigned column = 0; column < kDimension; ++column) {
    unsigned pivot = column;
    while (pivot < kDimension &&
           ((augmented[pivot] >> column) & 1U) == 0) {
      ++pivot;
    }
    if (pivot == kDimension) Fail("singular restricted generator");
    std::swap(augmented[column], augmented[pivot]);
    for (unsigned row = 0; row < kDimension; ++row) {
      if (row != column &&
          ((augmented[row] >> column) & 1U)) {
        augmented[row] ^= augmented[column];
      }
    }
  }
  Matrix inverse{};
  const std::uint64_t mask =
      (std::uint64_t{1} << kDimension) - 1;
  for (unsigned row = 0; row < kDimension; ++row) {
    if ((augmented[row] & mask) !=
        (std::uint64_t{1} << row)) {
      Fail("restricted inverse reduction failed");
    }
    inverse[row] = static_cast<std::uint32_t>(
        augmented[row] >> kDimension);
  }
  return inverse;
}

struct InputData {
  std::string classId;
  unsigned classPosition = 0;
  std::uint64_t classSize = 0;
  std::uint64_t centralizerOrder = 0;
  std::array<std::uint64_t, kAmbientDimension> representative{};
  std::array<std::uint64_t, kDimension> ambientBasis{};
  std::array<std::array<std::uint8_t, kDimension>,
             kDimension> f4Basis{};
  std::vector<Matrix> generators;
};

InputData Load(const std::string& path) {
  std::ifstream input(path);
  if (!input) Fail("cannot open " + path);
  InputData data;
  unsigned declaredDimension = 0;
  unsigned declaredAmbient = 0;
  unsigned declaredGenerators = 0;
  unsigned representativeRows = 0;
  unsigned ambientBasisRows = 0;
  unsigned f4BasisRows = 0;
  enum class Section {
    kRecords, kRepresentative, kAmbientBasis, kF4Basis, kGenerator
  };
  Section section = Section::kRecords;
  Matrix currentGenerator{};
  unsigned generatorRows = 0;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    const std::vector<std::string> fields = Split(line);
    if (section == Section::kRepresentative &&
        representativeRows < kAmbientDimension) {
      data.representative[representativeRows++] =
          ParseBinary(fields.at(0), kAmbientDimension,
                      "representative row");
      if (representativeRows == kAmbientDimension) {
        section = Section::kRecords;
      }
      continue;
    }
    if (section == Section::kAmbientBasis &&
        ambientBasisRows < kDimension) {
      data.ambientBasis[ambientBasisRows++] =
          ParseBinary(fields.at(0), kAmbientDimension,
                      "ambient basis row");
      if (ambientBasisRows == kDimension) {
        section = Section::kRecords;
      }
      continue;
    }
    if (section == Section::kF4Basis &&
        f4BasisRows < kDimension) {
      data.f4Basis[f4BasisRows++] = ParseF4(fields.at(0));
      if (f4BasisRows == kDimension) {
        section = Section::kRecords;
      }
      continue;
    }
    if (section == Section::kGenerator &&
        generatorRows < kDimension) {
      currentGenerator[generatorRows++] =
          static_cast<std::uint32_t>(
              ParseBinary(fields.at(0), kDimension,
                          "restricted generator row"));
      if (generatorRows == kDimension) {
        data.generators.push_back(currentGenerator);
        currentGenerator.fill(0);
        generatorRows = 0;
        section = Section::kRecords;
      }
      continue;
    }

    if (fields[0] == "SCHEMA") {
      if (fields.size() != 2 ||
          fields[1] != "FI22_OUTER_ACTION_V1") {
        Fail("unexpected input schema");
      }
    } else if (fields[0] == "CLASS") {
      if (fields.size() != 2) Fail("malformed CLASS record");
      data.classId = fields[1];
    } else if (fields[0] == "CLASS_POSITION") {
      data.classPosition = static_cast<unsigned>(
          ParseUnsigned(fields.at(1), "class position"));
    } else if (fields[0] == "CLASS_SIZE") {
      data.classSize =
          ParseUnsigned(fields.at(1), "class size");
    } else if (fields[0] == "CENT") {
      data.centralizerOrder =
          ParseUnsigned(fields.at(1), "centralizer order");
    } else if (fields[0] == "DIM") {
      declaredDimension = static_cast<unsigned>(
          ParseUnsigned(fields.at(1), "dimension"));
    } else if (fields[0] == "AMBIENT") {
      declaredAmbient = static_cast<unsigned>(
          ParseUnsigned(fields.at(1), "ambient dimension"));
    } else if (fields[0] == "NGENS") {
      declaredGenerators = static_cast<unsigned>(
          ParseUnsigned(fields.at(1), "generator count"));
    } else if (fields[0] == "REPRESENTATIVE_POWER") {
      if (fields.size() != 2) {
        Fail("malformed REPRESENTATIVE_POWER record");
      }
    } else if (fields[0] == "REPRESENTATIVE_MATRIX") {
      section = Section::kRepresentative;
    } else if (fields[0] == "BASIS_H") {
      section = Section::kAmbientBasis;
    } else if (fields[0] == "BASIS_F4") {
      section = Section::kF4Basis;
    } else if (fields[0] == "GEN") {
      if (fields.size() != 4 || fields[2] != "ADJUST") {
        Fail("malformed GEN record");
      }
      const unsigned expectedIndex =
          static_cast<unsigned>(data.generators.size() + 1);
      if (ParseUnsigned(fields[1], "generator index") !=
          expectedIndex) {
        Fail("nonconsecutive generator index");
      }
      section = Section::kGenerator;
    } else {
      Fail("unknown input record: " + fields[0]);
    }
  }
  if (section != Section::kRecords ||
      representativeRows != kAmbientDimension ||
      ambientBasisRows != kDimension ||
      f4BasisRows != kDimension ||
      declaredDimension != kDimension ||
      declaredAmbient != kAmbientDimension ||
      declaredGenerators != data.generators.size() ||
      data.classId.empty() || data.classPosition == 0 ||
      data.classSize == 0 || data.centralizerOrder == 0 ||
      data.generators.empty()) {
    Fail("incomplete fixed-action input");
  }
  std::vector<std::uint64_t> basisRows(
      data.ambientBasis.begin(), data.ambientBasis.end());
  if (BinaryRank(basisRows, kAmbientDimension) != kDimension) {
    Fail("ambient fixed-space basis is dependent");
  }
  for (const Matrix& generator : data.generators) {
    std::vector<std::uint64_t> rows(
        generator.begin(), generator.end());
    if (BinaryRank(rows, kDimension) != kDimension) {
      Fail("restricted generator is singular");
    }
  }
  return data;
}

struct GeneratorTable {
  std::array<std::array<std::uint32_t, 512>, 3> chunk{};

  explicit GeneratorTable(const Matrix& matrix) {
    for (unsigned block = 0; block < 3; ++block) {
      chunk[block][0] = 0;
      for (unsigned mask = 1; mask < 512; ++mask) {
        const unsigned bit = std::countr_zero(mask);
        chunk[block][mask] =
            chunk[block][mask & (mask - 1)] ^
            matrix[block * 9 + bit];
      }
    }
  }

  std::uint32_t Apply(std::uint32_t vector) const {
    return chunk[0][vector & kChunkMask] ^
        chunk[1][(vector >> 9) & kChunkMask] ^
        chunk[2][(vector >> 18) & kChunkMask];
  }
};

std::string CoordinateDigits(std::uint32_t vector) {
  std::string result(kDimension, '0');
  for (unsigned position = 0; position < kDimension; ++position) {
    if ((vector >> position) & 1U) result[position] = '1';
  }
  return result;
}

std::string AmbientDigits(
    std::uint32_t vector, const InputData& data) {
  std::uint64_t ambient = 0;
  while (vector != 0) {
    const unsigned bit = std::countr_zero(vector);
    ambient ^= data.ambientBasis[bit];
    vector &= vector - 1;
  }
  std::string result(kAmbientDimension, '0');
  for (unsigned position = 0;
       position < kAmbientDimension; ++position) {
    if ((ambient >> position) & 1U) result[position] = '1';
  }
  return result;
}

std::string F4Digits(
    std::uint32_t vector, const InputData& data) {
  std::array<std::uint8_t, kDimension> f4{};
  while (vector != 0) {
    const unsigned bit = std::countr_zero(vector);
    for (unsigned position = 0;
         position < kDimension; ++position) {
      f4[position] ^= data.f4Basis[bit][position];
    }
    vector &= vector - 1;
  }
  std::string result(kDimension, '0');
  for (unsigned position = 0; position < kDimension; ++position) {
    result[position] =
        static_cast<char>('0' + f4[position]);
  }
  return result;
}

class OrbitEnumerator {
 public:
  OrbitEnumerator(
      const InputData& data, unsigned threads)
      : data_(data),
        threadCount_(std::max(1U, threads)),
        visited_(std::make_unique<
                 std::atomic<std::uint64_t>[]>(kVisitedWords)) {
    for (std::size_t word = 0; word < kVisitedWords; ++word) {
      visited_[word].store(0, std::memory_order_relaxed);
    }
    std::vector<Matrix> matrices = data.generators;
    const std::size_t originalCount = matrices.size();
    for (std::size_t index = 0; index < originalCount; ++index) {
      const Matrix inverse = Inverse(matrices[index]);
      if (std::find(matrices.begin(), matrices.end(), inverse) ==
          matrices.end()) {
        matrices.push_back(inverse);
      }
    }
    tables_.reserve(matrices.size());
    for (const Matrix& matrix : matrices) {
      tables_.emplace_back(matrix);
    }
    std::cerr << "generators original=" << originalCount
              << " with_inverses=" << tables_.size() << "\n";
  }

  void Run(const std::string& outputPath) {
    std::ofstream output(outputPath);
    if (!output) Fail("cannot create " + outputPath);
    output << "# FI22_OUTER_ORBITS_V1"
           << " class=" << data_.classId
           << " class_position=" << data_.classPosition
           << " class_size=" << data_.classSize
           << " centralizer=" << data_.centralizerOrder
           << " fixed_dimension=" << kDimension << "\n";
    output << "# ORBIT representative_index orbit_size"
              " fixed_coordinates atlas_h54 core_f4\n";

    std::uint64_t orbitCount = 0;
    std::uint64_t visitedCount = 0;
    const auto start = std::chrono::steady_clock::now();
    for (std::uint32_t seed = 0; seed < kVectorCount; ++seed) {
      if (TestAndSet(seed)) continue;
      const std::uint64_t orbitSize = Explore(seed);
      if (data_.centralizerOrder % orbitSize != 0) {
        Fail("orbit size does not divide the exact centralizer order");
      }
      ++orbitCount;
      visitedCount += orbitSize;
      output << "ORBIT " << seed << ' ' << orbitSize << ' '
             << CoordinateDigits(seed) << ' '
             << AmbientDigits(seed, data_) << ' '
             << F4Digits(seed, data_) << '\n';
      if (orbitCount % 100 == 0) {
        const auto seconds =
            std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::steady_clock::now() - start).count();
        std::cerr << "progress class=" << data_.classId
                  << " orbits=" << orbitCount
                  << " visited=" << visitedCount
                  << " elapsed=" << seconds << "s\n";
      }
    }
    if (visitedCount != kVectorCount) {
      Fail("visited-vector total mismatch");
    }
    output << "TOTAL " << orbitCount << ' '
           << visitedCount << '\n';
    output.close();
    const auto seconds =
        std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - start).count();
    std::cerr << "DONE class=" << data_.classId
              << " orbits=" << orbitCount
              << " vectors=" << visitedCount
              << " elapsed=" << seconds << "s\n";
  }

 private:
  bool TestAndSet(std::uint32_t vector) {
    const std::uint64_t bit =
        std::uint64_t{1} << (vector & 63);
    const std::uint64_t old =
        visited_[vector >> 6].fetch_or(
            bit, std::memory_order_relaxed);
    return (old & bit) != 0;
  }

  std::uint64_t Explore(std::uint32_t seed) {
    std::vector<std::uint32_t> current{seed};
    std::vector<std::uint32_t> next;
    std::uint64_t size = 1;
    while (!current.empty()) {
      next.clear();
      if (threadCount_ == 1 || current.size() < 16384) {
        for (std::uint32_t vector : current) {
          for (const GeneratorTable& table : tables_) {
            const std::uint32_t image = table.Apply(vector);
            if (!TestAndSet(image)) next.push_back(image);
          }
        }
      } else {
        const unsigned workers = std::min<unsigned>(
            threadCount_, static_cast<unsigned>(current.size()));
        std::vector<std::vector<std::uint32_t>> parts(workers);
        std::vector<std::thread> threads;
        threads.reserve(workers);
        const std::size_t chunk =
            (current.size() + workers - 1) / workers;
        for (unsigned worker = 0; worker < workers; ++worker) {
          threads.emplace_back([&, worker]() {
            const std::size_t begin = worker * chunk;
            const std::size_t end =
                std::min(current.size(), begin + chunk);
            std::vector<std::uint32_t>& mine = parts[worker];
            for (std::size_t index = begin; index < end; ++index) {
              const std::uint32_t vector = current[index];
              for (const GeneratorTable& table : tables_) {
                const std::uint32_t image = table.Apply(vector);
                if (!TestAndSet(image)) mine.push_back(image);
              }
            }
          });
        }
        for (std::thread& thread : threads) thread.join();
        std::size_t total = 0;
        for (const auto& part : parts) total += part.size();
        next.reserve(total);
        for (auto& part : parts) {
          next.insert(next.end(), part.begin(), part.end());
        }
      }
      size += next.size();
      current.swap(next);
    }
    return size;
  }

  const InputData& data_;
  unsigned threadCount_;
  std::unique_ptr<std::atomic<std::uint64_t>[]> visited_;
  std::vector<GeneratorTable> tables_;
};

unsigned ParseThreads(const char* text) {
  const std::uint64_t value =
      ParseUnsigned(text, "thread count");
  if (value == 0 || value > 256) Fail("invalid thread count");
  return static_cast<unsigned>(value);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 4) {
    std::fprintf(stderr,
        "usage: %s INPUT OUTPUT THREADS\n", argv[0]);
    return 2;
  }
  try {
    const InputData data = Load(argv[1]);
    const unsigned threads = ParseThreads(argv[3]);
    OrbitEnumerator enumerator(data, threads);
    enumerator.Run(argv[2]);
  } catch (const std::exception& error) {
    std::fprintf(stderr, "ERROR: %s\n", error.what());
    return 1;
  }
  return 0;
}

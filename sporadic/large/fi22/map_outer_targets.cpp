// Re-express the 743 exact H54 outer-involution targets in the GF(4)^27
// coordinates used by the H-native WL result.
//
// The two original GAP exporters found module intertwiners independently,
// so their GF(4) coordinates can differ by a unit in End_K(V) = GF(4).
// The 54-bit ATLAS coordinates identify the vectors.  This program inverts the
// exported WL bridge over GF(2), checks every round trip, and writes a target
// table with exactly aligned GF(4) coordinates.  The program determines the
// scalar relating the two intertwiners and checks that it is uniform.
//
// Usage:
//   map_outer_targets ORIGINAL F4_TO_H54 INNER_SCREEN OUTPUT

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

constexpr int kCoordinates = 27;
constexpr int kBits = 54;
constexpr int kTargets = 743;
constexpr int kInnerTargets = 5363;

[[noreturn]] void Fail(const std::string& message) {
  std::cerr << "FATAL: " << message << '\n';
  std::exit(1);
}

std::vector<std::string> Split(const std::string& line) {
  std::istringstream input(line);
  std::vector<std::string> fields;
  std::string field;
  while (input >> field) fields.push_back(field);
  return fields;
}

std::uint64_t Bits(const std::string& digits) {
  if (digits.size() != kBits) Fail("a binary vector has the wrong length");
  std::uint64_t result = 0;
  for (int position = 0; position < kBits; ++position) {
    if (digits[position] == '1') {
      result |= std::uint64_t{1} << position;
    } else if (digits[position] != '0') {
      Fail("a binary vector contains a non-bit");
    }
  }
  return result;
}

std::string BitDigits(std::uint64_t value) {
  std::string result(kBits, '0');
  for (int position = 0; position < kBits; ++position) {
    if ((value >> position) & 1U) result[position] = '1';
  }
  return result;
}

bool F4Digits(const std::string& digits) {
  return digits.size() == kCoordinates &&
      std::all_of(digits.begin(), digits.end(),
          [](char digit) { return digit >= '0' && digit <= '3'; });
}

char TimesUnit(char digit, int power) {
  // Reference polynomial encoding: 1 -> 1, 2 -> z, 3 -> z^2.
  static constexpr char image[3][4] = {
      {'0', '1', '2', '3'},
      {'0', '2', '3', '1'},
      {'0', '3', '1', '2'}};
  if (digit < '0' || digit > '3') Fail("invalid GF(4) digit");
  if (power < 0 || power > 2) Fail("invalid GF(4) unit");
  return image[power][digit - '0'];
}

std::string TimesUnit(std::string vector, int power) {
  std::transform(vector.begin(), vector.end(), vector.begin(),
                 [power](char digit) { return TimesUnit(digit, power); });
  return vector;
}

struct Bridge {
  std::uint64_t basis[kBits]{};
  std::uint64_t inverseRows[kBits]{};

  explicit Bridge(const std::string& path) {
    std::ifstream input(path);
    if (!input) Fail("cannot open WL bridge " + path);
    std::string line;
    if (!std::getline(input, line) ||
        line != "coordinate\tdigit\th54_bits") {
      Fail("WL bridge header mismatch");
    }
    std::uint64_t contribution[kCoordinates][4]{};
    bool seen[kCoordinates][4]{};
    int rows = 0;
    while (std::getline(input, line)) {
      if (line.empty()) continue;
      const auto fields = Split(line);
      if (fields.size() != 3) Fail("malformed WL bridge row");
      const int coordinate = std::stoi(fields[0]);
      const int digit = std::stoi(fields[1]);
      if (coordinate < 1 || coordinate > kCoordinates ||
          digit < 1 || digit > 3 ||
          seen[coordinate - 1][digit]) {
        Fail("duplicate or out-of-range WL bridge row");
      }
      contribution[coordinate - 1][digit] = Bits(fields[2]);
      seen[coordinate - 1][digit] = true;
      ++rows;
    }
    if (rows != 3 * kCoordinates) Fail("WL bridge row count mismatch");
    for (int coordinate = 0; coordinate < kCoordinates; ++coordinate) {
      for (int digit = 1; digit <= 3; ++digit) {
        if (!seen[coordinate][digit]) Fail("incomplete WL bridge");
      }
      if ((contribution[coordinate][1] ^
           contribution[coordinate][2]) !=
          contribution[coordinate][3]) {
        Fail("WL bridge is not GF(2)-additive");
      }
      basis[2 * coordinate] = contribution[coordinate][1];
      basis[2 * coordinate + 1] = contribution[coordinate][2];
    }

    // Gauss-Jordan reduction.  inverseRows[column] is the combination of
    // original bridge-basis rows that produces the corresponding unit row.
    std::uint64_t reduced[kBits];
    for (int row = 0; row < kBits; ++row) {
      reduced[row] = basis[row];
      inverseRows[row] = std::uint64_t{1} << row;
    }
    for (int column = 0; column < kBits; ++column) {
      int pivot = column;
      while (pivot < kBits &&
             ((reduced[pivot] >> column) & 1U) == 0) {
        ++pivot;
      }
      if (pivot == kBits) Fail("WL bridge has rank below 54");
      std::swap(reduced[column], reduced[pivot]);
      std::swap(inverseRows[column], inverseRows[pivot]);
      for (int row = 0; row < kBits; ++row) {
        if (row != column && ((reduced[row] >> column) & 1U)) {
          reduced[row] ^= reduced[column];
          inverseRows[row] ^= inverseRows[column];
        }
      }
    }
    for (int row = 0; row < kBits; ++row) {
      if (reduced[row] != (std::uint64_t{1} << row)) {
        Fail("WL bridge reduction did not reach the identity");
      }
    }
  }

  std::uint64_t Coordinates(std::uint64_t target) const {
    std::uint64_t coordinates = 0;
    for (int bit = 0; bit < kBits; ++bit) {
      if ((target >> bit) & 1U) coordinates ^= inverseRows[bit];
    }
    return coordinates;
  }

  std::uint64_t Encode(std::uint64_t coordinates) const {
    std::uint64_t target = 0;
    for (int bit = 0; bit < kBits; ++bit) {
      if ((coordinates >> bit) & 1U) target ^= basis[bit];
    }
    return target;
  }

  std::string Decode(std::uint64_t target) const {
    const std::uint64_t coordinates = Coordinates(target);
    if (Encode(coordinates) != target) Fail("WL bridge round trip failed");
    std::string digits(kCoordinates, '0');
    for (int coordinate = 0; coordinate < kCoordinates; ++coordinate) {
      const int low = (coordinates >> (2 * coordinate)) & 1U;
      const int high = (coordinates >> (2 * coordinate + 1)) & 1U;
      digits[coordinate] = static_cast<char>('0' + low + 2 * high);
    }
    return digits;
  }
};

std::unordered_set<std::string> ReadInnerTargets(const std::string& path) {
  std::ifstream input(path);
  if (!input) Fail("cannot open inner witness screen " + path);
  std::unordered_set<std::string> targets;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') continue;
    const auto fields = Split(line);
    if (fields.size() != 14 || !F4Digits(fields[4])) {
      Fail("malformed inner witness row");
    }
    targets.insert(fields[4]);
  }
  if (targets.size() != kInnerTargets) {
    Fail("inner target count is not 5363");
  }
  return targets;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cerr << "usage: " << argv[0]
              << " ORIGINAL F4_TO_H54 INNER_SCREEN OUTPUT\n";
    return 1;
  }
  const Bridge bridge(argv[2]);
  const auto innerTargets = ReadInnerTargets(argv[3]);
  std::ifstream input(argv[1]);
  if (!input) Fail("cannot open original outer target table");
  std::ofstream output(argv[4]);
  if (!output) Fail("cannot create remapped outer target table");
  output << "# FI22_OUTER_TARGETS_IN_CORE_COORDINATES_V1\n";
  output << "# tid classes representative_indices orbit_sizes"
            " atlas_h54 core_f4 overlaps_inner\n";

  std::unordered_set<std::string> remappedTargets;
  std::string line;
  int targetCount = 0;
  int scalarChecks = 0;
  int innerOverlaps = 0;
  int scalarPower = -1;
  while (std::getline(input, line)) {
    if (line.rfind("TARGET ", 0) != 0) continue;
    const auto fields = Split(line);
    if (fields.size() != 8 || fields[0] != "TARGET" ||
        std::stoi(fields[1]) != targetCount ||
        fields[5].size() != kBits || !F4Digits(fields[6]) ||
        fields[7] != "0") {
      Fail("malformed original outer target row");
    }
    const std::uint64_t h54 = Bits(fields[5]);
    const std::string remapped = bridge.Decode(h54);
    if (remapped == std::string(kCoordinates, '0')) {
      Fail("zero outer target");
    }
    if (scalarPower < 0) {
      for (int power = 0; power < 3; ++power) {
        if (remapped == TimesUnit(fields[6], power)) {
          scalarPower = power;
          break;
        }
      }
      if (scalarPower < 0) {
        Fail("the two GF(4) coordinate maps do not differ by a scalar");
      }
    }
    if (remapped != TimesUnit(fields[6], scalarPower)) {
      Fail("the scalar relating the two GF(4) coordinate maps is not uniform");
    }
    ++scalarChecks;
    if (!remappedTargets.insert(remapped).second) {
      Fail("duplicate target after WL bridge remap");
    }
    const bool overlapsInner = innerTargets.count(remapped) != 0;
    if (overlapsInner) ++innerOverlaps;
    output << "TARGET " << targetCount << ' ' << fields[2] << ' '
           << fields[3] << ' ' << fields[4] << ' '
           << BitDigits(h54) << ' ' << remapped << ' '
           << (overlapsInner ? 1 : 0) << '\n';
    ++targetCount;
  }
  if (targetCount != kTargets || scalarChecks != kTargets ||
      remappedTargets.size() != kTargets || innerOverlaps != 0) {
    Fail("remapped outer target census mismatch");
  }
  output << "# done targets=743 scalar_power=" << scalarPower
            << " scalar_checks=743 unique=743"
            " inner_overlaps=0\n";
  std::cout << "Mapped 743 outer targets; the coordinate maps differ by z^"
            << scalarPower << ".\n";
  return 0;
}

// Exact merge of the three outer-involution fixed-space orbit inventories
// for H = 3.Fi22.2, with an overlap check against the existing inner-class
// target inventory.
//
// Each outer orbit representative is recorded in both the ATLAS H54
// coordinates and the reference F4^27 core coordinates.  The two coordinate
// maps are checked to be one-to-one on every record seen.  Exact vectors,
// rather than orbit labels, are then deduplicated across the three outer
// involution families.
//
// Usage:
//   outer_target_merge OUT1 OUT2 OUT3 INNER_SCREEN OUTPUT

#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr std::uint64_t kFixedSpaceSize = std::uint64_t{1} << 27;
constexpr std::size_t kF4Length = 27;
constexpr std::size_t kH54Length = 54;
constexpr std::size_t kExpectedInnerTargets = 5363;

[[noreturn]] void Fail(const std::string& message) {
  throw std::runtime_error(message);
}

std::vector<std::string> Split(const std::string& line) {
  std::istringstream input(line);
  std::vector<std::string> fields;
  std::string field;
  while (input >> field) fields.push_back(field);
  return fields;
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

void CheckDigits(
    const std::string& digits, std::size_t length,
    char maximum, const std::string& label) {
  if (digits.size() != length) {
    Fail(label + " has the wrong length");
  }
  for (const char digit : digits) {
    if (digit < '0' || digit > maximum) {
      Fail(label + " has an invalid digit");
    }
  }
}

struct Occurrence {
  std::string classId;
  std::uint64_t representativeIndex = 0;
  std::uint64_t orbitSize = 0;
};

struct Target {
  std::string h54;
  std::string f4;
  std::vector<Occurrence> occurrences;
};

struct InventorySummary {
  std::string classId;
  std::uint64_t rows = 0;
  std::uint64_t nonzeroRows = 0;
  std::uint64_t coveredVectors = 0;
};

InventorySummary ReadOuter(
    const std::string& path,
    std::map<std::string, Target>& targets,
    std::unordered_map<std::string, std::string>& h54ToF4,
    std::unordered_map<std::string, std::string>& f4ToH54) {
  std::ifstream input(path);
  if (!input) Fail("cannot open " + path);
  InventorySummary summary;
  bool sawHeader = false;
  bool sawTotal = false;
  bool sawZero = false;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    if (line.rfind("# FI22_OUTER_ORBITS_V1 ", 0) == 0) {
      if (sawHeader) Fail("duplicate inventory header in " + path);
      sawHeader = true;
      const std::vector<std::string> fields = Split(line);
      for (const std::string& field : fields) {
        if (field.rfind("class=", 0) == 0) {
          summary.classId = field.substr(6);
        }
      }
      continue;
    }
    if (line[0] == '#') continue;
    const std::vector<std::string> fields = Split(line);
    if (fields.empty()) continue;
    if (fields[0] == "ORBIT") {
      if (fields.size() != 6) {
        Fail("malformed ORBIT row in " + path);
      }
      const std::uint64_t representativeIndex =
          ParseUnsigned(fields[1], "representative index");
      const std::uint64_t orbitSize =
          ParseUnsigned(fields[2], "orbit size");
      CheckDigits(fields[3], kF4Length, '1',
                  "fixed coordinates");
      CheckDigits(fields[4], kH54Length, '1',
                  "ATLAS H54 vector");
      CheckDigits(fields[5], kF4Length, '3',
                  "reference F4 vector");
      ++summary.rows;
      summary.coveredVectors += orbitSize;
      if (representativeIndex == 0) {
        if (sawZero || orbitSize != 1 ||
            fields[3] != std::string(kF4Length, '0') ||
            fields[4] != std::string(kH54Length, '0') ||
            fields[5] != std::string(kF4Length, '0')) {
          Fail("invalid zero orbit in " + path);
        }
        sawZero = true;
        continue;
      }
      ++summary.nonzeroRows;

      const auto hInsert = h54ToF4.emplace(fields[4], fields[5]);
      if (!hInsert.second && hInsert.first->second != fields[5]) {
        Fail("one H54 vector has two reference F4 encodings");
      }
      const auto fInsert = f4ToH54.emplace(fields[5], fields[4]);
      if (!fInsert.second && fInsert.first->second != fields[4]) {
        Fail("one reference F4 vector has two H54 encodings");
      }

      Target& target = targets[fields[5]];
      if (target.f4.empty()) {
        target.h54 = fields[4];
        target.f4 = fields[5];
      } else if (target.h54 != fields[4]) {
        Fail("exact-vector deduplication bridge mismatch");
      }
      target.occurrences.push_back(
          Occurrence{summary.classId, representativeIndex, orbitSize});
      continue;
    }
    if (fields[0] == "TOTAL") {
      if (fields.size() != 3 || sawTotal) {
        Fail("malformed TOTAL row in " + path);
      }
      if (ParseUnsigned(fields[1], "total orbit count") !=
              summary.rows ||
          ParseUnsigned(fields[2], "total vector count") !=
              summary.coveredVectors) {
        Fail("TOTAL row does not match inventory rows in " + path);
      }
      sawTotal = true;
      continue;
    }
    Fail("unknown inventory row in " + path);
  }
  if (!sawHeader || !sawTotal || !sawZero ||
      summary.classId.empty() ||
      summary.coveredVectors != kFixedSpaceSize) {
    Fail("incomplete outer inventory " + path);
  }
  return summary;
}

std::set<std::string> ReadInnerTargets(const std::string& path) {
  std::ifstream input(path);
  if (!input) Fail("cannot open " + path);
  std::set<std::string> result;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') continue;
    const std::vector<std::string> fields = Split(line);
    if (fields.size() < 5) {
      Fail("malformed inner witness-screen row");
    }
    CheckDigits(fields[4], kF4Length, '3',
                "inner reference F4 target");
    result.insert(fields[4]);
  }
  if (result.size() != kExpectedInnerTargets) {
    Fail("inner target inventory has " +
         std::to_string(result.size()) +
         " exact vectors, expected " +
         std::to_string(kExpectedInnerTargets));
  }
  return result;
}

std::string JoinOccurrenceField(
    const std::vector<Occurrence>& occurrences,
    int field) {
  std::string result;
  for (std::size_t index = 0;
       index < occurrences.size(); ++index) {
    if (index != 0) result += ',';
    if (field == 0) {
      result += occurrences[index].classId;
    } else if (field == 1) {
      result +=
          std::to_string(occurrences[index].representativeIndex);
    } else {
      result += std::to_string(occurrences[index].orbitSize);
    }
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 6) {
    std::cerr
        << "usage: " << argv[0]
        << " OUT1 OUT2 OUT3 INNER_SCREEN OUTPUT\n";
    return 2;
  }
  try {
    std::map<std::string, Target> targets;
    std::unordered_map<std::string, std::string> h54ToF4;
    std::unordered_map<std::string, std::string> f4ToH54;
    std::array<InventorySummary, 3> inventories;
    for (std::size_t index = 0; index < inventories.size(); ++index) {
      inventories[index] = ReadOuter(
          argv[index + 1], targets, h54ToF4, f4ToH54);
      const std::string expectedClass =
          "OUT" + std::to_string(index + 1);
      if (inventories[index].classId != expectedClass) {
        Fail("outer inventories are not in OUT1,OUT2,OUT3 order");
      }
    }
    const std::set<std::string> innerTargets =
        ReadInnerTargets(argv[4]);

    std::uint64_t rawNonzero = 0;
    for (const InventorySummary& inventory : inventories) {
      rawNonzero += inventory.nonzeroRows;
    }
    std::uint64_t overlapInner = 0;
    std::uint64_t crossFamilyDuplicateRows = 0;
    for (const auto& [f4, target] : targets) {
      if (innerTargets.contains(f4)) ++overlapInner;
      crossFamilyDuplicateRows += target.occurrences.size() - 1;
    }
    if (rawNonzero != targets.size() + crossFamilyDuplicateRows) {
      Fail("outer exact-vector accounting mismatch");
    }

    std::ofstream output(argv[5]);
    if (!output) Fail("cannot create output");
    output << "# FI22_OUTER_TARGETS_V1\n";
    output << "# tid classes representative_indices orbit_sizes"
              " atlas_h54 core_f4 overlaps_inner\n";
    std::uint64_t targetId = 0;
    for (const auto& [f4, target] : targets) {
      output << "TARGET " << targetId++ << ' '
             << JoinOccurrenceField(target.occurrences, 0) << ' '
             << JoinOccurrenceField(target.occurrences, 1) << ' '
             << JoinOccurrenceField(target.occurrences, 2) << ' '
             << target.h54 << ' ' << target.f4 << ' '
             << (innerTargets.contains(f4) ? 1 : 0) << '\n';
    }
    const std::uint64_t unionNonregular =
        innerTargets.size() + targets.size() - overlapInner;
    output << "TOTAL"
           << " raw_outer_nonzero " << rawNonzero
           << " unique_outer " << targets.size()
           << " cross_family_duplicate_rows "
           << crossFamilyDuplicateRows
           << " inner_unique " << innerTargets.size()
           << " outer_inner_overlap " << overlapInner
           << " union_nonregular " << unionNonregular
           << '\n';
    output.close();

    std::cout << "OUTER_MERGE";
    for (const InventorySummary& inventory : inventories) {
      std::cout << ' ' << inventory.classId
                << "_orbits=" << inventory.rows
                << ' ' << inventory.classId
                << "_nonzero=" << inventory.nonzeroRows;
    }
    std::cout << " raw_outer_nonzero=" << rawNonzero
              << " unique_outer=" << targets.size()
              << " cross_family_duplicate_rows="
              << crossFamilyDuplicateRows
              << " inner_unique=" << innerTargets.size()
              << " outer_inner_overlap=" << overlapInner
              << " union_nonregular=" << unionNonregular
              << " coordinate_bridge=BIJECTIVE_ON_RECORDS\n";
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
  return 0;
}

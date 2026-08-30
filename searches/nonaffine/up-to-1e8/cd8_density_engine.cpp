// Exact regular-suborbit density engine for the 39 compound-diagonal CD8
// actions.  The state space is (A5^2)^2, represented by four A5 elements.

#define main diagonal_starstar_small_embedded_main
#include "../../diagonal/shared/diagonal_starstar_small.cpp"
#undef main

#include <fstream>
#include <set>
#include <sstream>
#include <tuple>

namespace {

using CDCode = std::uint32_t;
using CDDigits = std::array<Element, 4>;

struct SignedPermutation6 {
  std::array<std::uint8_t, 6> target{};
  std::array<std::uint8_t, 6> outer{};

  bool operator<(const SignedPermutation6& other) const {
    return std::tie(target, outer) < std::tie(other.target, other.outer);
  }
};

struct CDCase {
  int id = 0;
  int quotient_order = 0;
  int stabiliser_order = 0;
  std::vector<SignedPermutation6> generators;
};

std::vector<std::string> split(const std::string& text, char delimiter) {
  std::vector<std::string> result;
  std::stringstream stream(text);
  std::string item;
  while (std::getline(stream, item, delimiter)) result.push_back(item);
  return result;
}

int assignment_int(const std::string& field, const std::string& name) {
  const std::string prefix = name + "=";
  if (field.rfind(prefix, 0) != 0) {
    throw std::runtime_error("expected " + prefix);
  }
  return std::stoi(field.substr(prefix.size()));
}

std::array<std::uint8_t, 6> six_values(
    const std::string& field, const std::string& name, int minimum,
    int maximum) {
  const std::string prefix = name + "=";
  if (field.rfind(prefix, 0) != 0) {
    throw std::runtime_error("expected " + prefix);
  }
  const auto pieces = split(field.substr(prefix.size()), ',');
  if (pieces.size() != 6) throw std::runtime_error("expected six values");
  std::array<std::uint8_t, 6> answer{};
  std::array<std::uint8_t, 7> seen{};
  for (int i = 0; i < 6; ++i) {
    const int value = std::stoi(pieces[i]);
    if (value < minimum || value > maximum) {
      throw std::runtime_error("descriptor value outside range");
    }
    answer[i] = static_cast<std::uint8_t>(value - minimum);
    if (name == "perm") {
      if (seen[value]) throw std::runtime_error("descriptor is not a permutation");
      seen[value] = 1;
    }
  }
  return answer;
}

CDCase load_case(const std::string& path, int selected) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open descriptor file");
  std::string line;
  if (!std::getline(input, line) || line != "CD8_DESCRIPTORS_V1") {
    throw std::runtime_error("bad descriptor header");
  }
  CDCase answer;
  int expected_generators = -1;
  bool complete = false;
  bool terminal = false;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    const auto fields = split(line, '|');
    if (fields[0] == "CD_CASE") {
      if (fields.size() != 6) throw std::runtime_error("bad CD_CASE row");
      const int id = assignment_int(fields[1], "case");
      if (id != selected) continue;
      answer.id = id;
      answer.quotient_order = assignment_int(fields[2], "B");
      answer.stabiliser_order = assignment_int(fields[3], "H");
      if (assignment_int(fields[4], "degree") != 12960000) {
        throw std::runtime_error("bad descriptor degree");
      }
      expected_generators = assignment_int(fields[5], "generators");
    } else if (fields[0] == "CD_GEN") {
      if (fields.size() != 5) throw std::runtime_error("bad CD_GEN row");
      const int id = assignment_int(fields[1], "case");
      if (id != selected) continue;
      if (answer.id != selected) throw std::runtime_error("generator before case");
      const int sequence = assignment_int(fields[2], "sequence");
      if (sequence != static_cast<int>(answer.generators.size()) + 1) {
        throw std::runtime_error("generator sequence mismatch");
      }
      SignedPermutation6 generator;
      generator.target = six_values(fields[3], "perm", 1, 6);
      generator.outer = six_values(fields[4], "outer", 0, 1);
      answer.generators.push_back(generator);
    } else if (fields[0] == "CD_CASE_COMPLETE") {
      if (fields.size() != 4) throw std::runtime_error("bad case completion");
      const int id = assignment_int(fields[1], "case");
      if (id != selected) continue;
      if (assignment_int(fields[2], "B") != answer.quotient_order ||
          assignment_int(fields[3], "signed_closure") != answer.quotient_order) {
        throw std::runtime_error("descriptor closure mismatch");
      }
      complete = true;
    } else if (fields[0] == "CD8_DESCRIPTORS_COMPLETE") {
      if (fields.size() != 3 || assignment_int(fields[1], "cases") != 39 ||
          assignment_int(fields[2], "degree") != 12960000) {
        throw std::runtime_error("bad descriptor terminal row");
      }
      terminal = true;
    } else {
      throw std::runtime_error("unknown descriptor row");
    }
  }
  if (answer.id != selected || expected_generators < 1 ||
      expected_generators != static_cast<int>(answer.generators.size()) ||
      !complete || !terminal || answer.stabiliser_order != 3600 * answer.quotient_order) {
    throw std::runtime_error("selected descriptor is incomplete");
  }
  return answer;
}

SignedPermutation6 signed_identity() {
  SignedPermutation6 answer;
  for (int i = 0; i < 6; ++i) answer.target[i] = i;
  return answer;
}

SignedPermutation6 signed_multiply(
    const SignedPermutation6& left, const SignedPermutation6& right) {
  SignedPermutation6 answer;
  for (int source = 0; source < 6; ++source) {
    answer.target[source] = right.target[left.target[source]];
    answer.outer[source] = static_cast<std::uint8_t>(
        left.outer[source] ^ right.outer[left.target[source]]);
  }
  return answer;
}

int signed_closure_order(const std::vector<SignedPermutation6>& generators) {
  std::set<SignedPermutation6> closure{signed_identity()};
  std::vector<SignedPermutation6> queue{signed_identity()};
  for (std::size_t head = 0; head < queue.size(); ++head) {
    for (const auto& generator : generators) {
      const auto product = signed_multiply(queue[head], generator);
      if (closure.insert(product).second) queue.push_back(product);
    }
  }
  return static_cast<int>(closure.size());
}

constexpr CDCode kCDDegree = 12960000;

CDDigits decode_cd(CDCode code) {
  CDDigits answer{};
  for (int i = 0; i < 4; ++i) {
    answer[i] = static_cast<Element>(code % 60);
    code /= 60;
  }
  return answer;
}

CDCode encode_cd(const CDDigits& digits) {
  CDCode answer = 0;
  CDCode place = 1;
  for (Element value : digits) {
    answer += place * value;
    place *= 60;
  }
  return answer;
}

CDCode apply_inner(
    const FiniteGroup& group, CDCode code, int component,
    Element conjugator) {
  CDDigits digits = decode_cd(code);
  const Element inverse = group.inverse[conjugator];
  const int offset = 2 * component;
  for (int i = 0; i < 2; ++i) {
    digits[offset + i] = group.multiply(
        group.multiply(inverse, digits[offset + i]), conjugator);
  }
  return encode_cd(digits);
}

CDCode apply_signed(
    const FiniteGroup& group, CDCode code,
    const SignedPermutation6& action) {
  const CDDigits digits = decode_cd(code);
  std::array<Element, 6> input{group.identity, digits[0], digits[1],
                               group.identity, digits[2], digits[3]};
  std::array<Element, 6> output{};
  for (int source = 0; source < 6; ++source) {
    Element value = input[source];
    if (action.outer[source]) value = group.outer_map[value];
    output[action.target[source]] = value;
  }
  CDDigits normalised{};
  for (int component = 0; component < 2; ++component) {
    const int base = 3 * component;
    const Element inverse = group.inverse[output[base]];
    normalised[2 * component] = group.multiply(inverse, output[base + 1]);
    normalised[2 * component + 1] = group.multiply(inverse, output[base + 2]);
  }
  return encode_cd(normalised);
}

struct ArgumentsCD {
  std::string descriptor;
  int selected = 0;
  bool graph = false;
};

ArgumentsCD parse_cd_arguments(int argc, char** argv) {
  ArgumentsCD answer;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--descriptor" && index + 1 < argc) {
      answer.descriptor = argv[++index];
    } else if (argument == "--case" && index + 1 < argc) {
      answer.selected = std::stoi(argv[++index]);
    } else if (argument == "--graph") {
      answer.graph = true;
    } else {
      throw std::runtime_error(
          "usage: cd8_density_engine --descriptor FILE --case 1..39 [--graph]");
    }
  }
  if (answer.descriptor.empty() || answer.selected < 1 || answer.selected > 39) {
    throw std::runtime_error("bad CD8 arguments");
  }
  return answer;
}

CDCode relative_cd(const FiniteGroup& group, CDCode point, CDCode origin) {
  const CDDigits point_digits = decode_cd(point);
  const CDDigits origin_digits = decode_cd(origin);
  CDDigits answer{};
  for (int index = 0; index < 4; ++index) {
    answer[index] = group.multiply(
        point_digits[index], group.inverse[origin_digits[index]]);
  }
  return encode_cd(answer);
}

void checksum_mix(std::uint64_t& checksum, CDCode value) {
  for (int byte = 0; byte < 4; ++byte) {
    checksum ^= static_cast<std::uint8_t>(value >> (8 * byte));
    checksum *= 1099511628211ULL;
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_cd_arguments(argc, argv);
    const CDCase spec = load_case(arguments.descriptor, arguments.selected);
    if (signed_closure_order(spec.generators) != spec.quotient_order) {
      throw std::runtime_error("C++ signed quotient closure mismatch");
    }
    FiniteGroup group = alternating_group(5);
    verify_group(group);
    if (group.order != 60 || group.outer_order != 2) {
      throw std::runtime_error("bad A5 model");
    }

    std::vector<std::int32_t> orbit(kCDDegree, -1);
    std::vector<std::uint8_t> regular(kCDDegree, 0);
    std::vector<CDCode> queue;
    std::vector<CDCode> representatives;
    std::uint64_t regular_points = 0;
    std::uint64_t regular_orbits = 0;
    std::uint64_t orbit_count = 0;
    std::uint64_t transition_checks = 0;
    for (CDCode start = 0; start < kCDDegree; ++start) {
      if (orbit[start] >= 0) continue;
      representatives.push_back(start);
      const std::int32_t orbit_id = static_cast<std::int32_t>(orbit_count++);
      queue.clear();
      queue.push_back(start);
      orbit[start] = orbit_id;
      for (std::size_t head = 0; head < queue.size(); ++head) {
        const CDCode point = queue[head];
        for (int component = 0; component < 2; ++component) {
          for (Element generator : group.generators) {
            const CDCode target = apply_inner(group, point, component, generator);
            ++transition_checks;
            if (orbit[target] < 0) {
              orbit[target] = orbit_id;
              queue.push_back(target);
            }
          }
        }
        for (const auto& generator : spec.generators) {
          const CDCode target = apply_signed(group, point, generator);
          ++transition_checks;
          if (orbit[target] < 0) {
            orbit[target] = orbit_id;
            queue.push_back(target);
          }
        }
      }
      if (queue.size() > static_cast<std::size_t>(spec.stabiliser_order) ||
          spec.stabiliser_order % static_cast<int>(queue.size()) != 0) {
        throw std::runtime_error("orbit size does not divide H order");
      }
      if (queue.size() == static_cast<std::size_t>(spec.stabiliser_order)) {
        ++regular_orbits;
        regular_points += queue.size();
        for (CDCode point : queue) regular[point] = 1;
      }
    }
    if (regular_orbits == 0) {
      throw std::runtime_error("no regular orbit to certify H order");
    }
    const bool density = 2 * regular_points > kCDDegree;
    if (arguments.graph && !density) {
      bool common_all = true;
      CDCode counterexample_rep = 0;
      std::uint64_t intersection_checks = 0;
      std::uint64_t witness_checksum = 1469598103934665603ULL;
      std::uint64_t representatives_checked = 0;
      for (CDCode representative : representatives) {
        bool found = false;
        for (CDCode witness = 0; witness < kCDDegree; ++witness) {
          if (!regular[witness]) continue;
          ++intersection_checks;
          const CDCode relative = relative_cd(group, witness, representative);
          if (regular[relative]) {
            checksum_mix(witness_checksum, representative);
            checksum_mix(witness_checksum, witness);
            checksum_mix(witness_checksum, relative);
            found = true;
            break;
          }
        }
        ++representatives_checked;
        if (!found) {
          common_all = false;
          counterexample_rep = representative;
          break;
        }
      }
      std::cout << "CD8_GRAPH_V1\n"
                << "CD8_GRAPH_RESULT|case=" << spec.id
                << "|degree=" << kCDDegree
                << "|B=" << spec.quotient_order
                << "|H=" << spec.stabiliser_order
                << "|quotient_closure=" << signed_closure_order(spec.generators)
                << "|H_orbits=" << orbit_count
                << "|regular_orbits=" << regular_orbits
                << "|regular_points=" << regular_points
                << "|orbit_reps_checked=" << representatives_checked
                << "|common_all=" << (common_all ? "true" : "false")
                << "|counterexample_rep=" << counterexample_rep
                << "|witness_checksum=" << witness_checksum
                << "|intersection_checks=" << intersection_checks
                << "|transition_checks=" << transition_checks << '\n'
                << "CD8_GRAPH_COMPLETE|case=" << spec.id
                << "|common_neighbour=" << (common_all ? "true" : "false")
                << "|method=exact_H_orbit_intersections\n";
      return 0;
    }
    std::cout << "CD8_DENSITY_V1\n"
              << "CD8_RESULT|case=" << spec.id
              << "|degree=" << kCDDegree
              << "|B=" << spec.quotient_order
              << "|H=" << spec.stabiliser_order
              << "|quotient_closure=" << signed_closure_order(spec.generators)
              << "|H_orbits=" << orbit_count
              << "|regular_orbits=" << regular_orbits
              << "|regular_points=" << regular_points
              << "|density_num=" << regular_points
              << "|density_den=" << kCDDegree
              << "|gt_half=" << (density ? "true" : "false")
              << "|transition_checks=" << transition_checks << '\n'
              << "CD8_COMPLETE|case=" << spec.id
              << "|common_neighbour=" << (density ? "true" : "unresolved")
              << "|method=exact_regular_density\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "CD8_ERROR " << error.what() << '\n';
    return 1;
  }
}

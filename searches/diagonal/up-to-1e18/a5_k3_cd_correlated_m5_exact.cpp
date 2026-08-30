// Exact common-neighbour witnesses for the residual correlated A5,k=3,m=5
// quotient classes.  Magma supplies one representative of every quotient
// orbit by a double-coset calculation; this engine independently rebuilds
// the quotient group and the 55-colour/77-target incidence relation.

#define A5_CORRELATED_M3_EMBEDDED
#include "a5_k3_cd_correlated_m3_exact.cpp"
#undef A5_CORRELATED_M3_EMBEDDED

#include <bit>
#include <charconv>
#include <limits>

namespace {

using Tuple5 = std::array<std::uint8_t, 5>;

struct W5 {
  std::array<std::uint8_t, 5> components{};
  Tuple5 top{};
};

struct Case5 {
  int number = 0;
  int m = 0;
  int order = 0;
  std::uint64_t regular_points = 0;
  std::vector<W5> generators;
  std::vector<Tuple5> x_reps;
  std::vector<Tuple5> y_reps;
};

std::string field_value(const std::string& line, const std::string& key) {
  const std::string needle = "|" + key + "=";
  const auto begin = line.find(needle);
  if (begin == std::string::npos) return {};
  const auto value_begin = begin + needle.size();
  const auto end = line.find('|', value_begin);
  return line.substr(value_begin, end == std::string::npos
                                      ? std::string::npos : end - value_begin);
}

std::uint64_t unsigned_field(const std::string& line, const std::string& key) {
  const std::string value = field_value(line, key);
  if (value.empty()) throw std::runtime_error("missing field " + key);
  std::uint64_t answer = 0;
  const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), answer);
  if (error != std::errc{} || end != value.data() + value.size()) {
    throw std::runtime_error("bad integer field " + key);
  }
  return answer;
}

Tuple5 tuple5_from_numbers(const std::string& text, int width, int offset = 0) {
  const auto values = parse_numbers(text);
  if (values.size() != static_cast<std::size_t>(width + offset)) {
    throw std::runtime_error("bad tuple descriptor");
  }
  Tuple5 answer{0,1,2,3,4};
  for (int i = 0; i < width; ++i) answer[i] = values[i + offset] - 1;
  return answer;
}

std::vector<Case5> parse_representatives(
    const std::string& path, const CoreData& core) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open orbit-representative transcript");
  std::map<int, Case5> cases;
  std::string line;
  int width = 0;
  while (std::getline(input, line)) {
    if (line.starts_with("ORBIT_REP_SHAPE|")) {
      width = unsigned_field(line, "m");
      if (width < 3 || width > 5) throw std::runtime_error("bad compound width");
    } else if (line.starts_with("ORBIT_REP_CASE|")) {
      if (width == 0) throw std::runtime_error("case before shape");
      Case5 item;
      item.number = unsigned_field(line, "case");
      item.m = width;
      item.order = unsigned_field(line, "M");
      item.regular_points = unsigned_field(line, "regular");
      cases[item.number] = item;
    } else if (line.starts_with("ORBIT_REP_GENERATOR|")) {
      const int number = unsigned_field(line, "case");
      if (!cases.contains(number)) throw std::runtime_error("generator before case");
      const auto top_begin = line.find("|top=");
      const auto component_begin = line.find("|components=");
      if (top_begin == std::string::npos || component_begin == std::string::npos) {
        throw std::runtime_error("bad generator row");
      }
      W5 generator;
      generator.top = tuple5_from_numbers(
          line.substr(top_begin, component_begin - top_begin), width);
      generator.components.fill(core.s3_index.at(SmallPerm{0,1,2}));
      std::string components = line.substr(component_begin + 12);
      if (components.front() != '[' || components.back() != ']') {
        throw std::runtime_error("bad components bracket");
      }
      components = components.substr(1, components.size() - 2);
      std::stringstream stream(components);
      std::string component;
      int coordinate = 0;
      while (std::getline(stream, component, ';')) {
        const auto values = parse_numbers(component);
        if (values.size() != 4 || coordinate >= width) {
          throw std::runtime_error("bad component descriptor");
        }
        SmallPerm permutation{};
        for (int i = 0; i < 3; ++i) permutation[i] = values[i + 1] - 1;
        generator.components[coordinate++] =
            values[0] * 6 + core.s3_index.at(permutation);
      }
      if (coordinate != width) throw std::runtime_error("bad component count");
      cases[number].generators.push_back(generator);
    } else if (line.starts_with("X_REGULAR_REP|")) {
      const int number = unsigned_field(line, "case");
      cases.at(number).x_reps.push_back(
          tuple5_from_numbers(field_value(line, "tuple"), width));
    } else if (line.starts_with("Y_TARGET_REP|")) {
      const int number = unsigned_field(line, "case");
      cases.at(number).y_reps.push_back(
          tuple5_from_numbers(field_value(line, "tuple"), width));
    } else if (line.starts_with("ORBIT_REP_CASE_COMPLETE|")) {
      const int number = unsigned_field(line, "case");
      const auto& item = cases.at(number);
      if (item.x_reps.size() != unsigned_field(line, "x") ||
          item.y_reps.size() != unsigned_field(line, "y")) {
        throw std::runtime_error("representative census mismatch");
      }
    }
  }
  std::vector<Case5> answer;
  for (auto& [number, item] : cases) {
    if (item.generators.empty() || item.x_reps.empty() || item.y_reps.empty()) {
      throw std::runtime_error("incomplete residual case");
    }
    answer.push_back(std::move(item));
  }
  return answer;
}

Tuple5 compose5(const Tuple5& left, const Tuple5& right, int width) {
  Tuple5 answer{};
  for (int i = 0; i < 5; ++i) answer[i] = i;
  for (int i = 0; i < width; ++i) answer[i] = right[left[i]];
  return answer;
}

W5 multiply5(
    const W5& left, const W5& right, const CoreData& core, int width) {
  W5 answer;
  answer.top = compose5(left.top, right.top, width);
  answer.components.fill(core.s3_index.at(SmallPerm{0,1,2}));
  for (int i = 0; i < width; ++i) {
    answer.components[i] = r_product(
        left.components[i], right.components[left.top[i]], core);
  }
  return answer;
}

std::uint64_t key5(const W5& element, int width) {
  std::uint64_t key = width;
  for (int i = 0; i < width; ++i) key = key * 5 + element.top[i];
  for (int i = 0; i < width; ++i) key = key * 12 + element.components[i];
  return key;
}

std::vector<W5> closure5(const Case5& item, const CoreData& core) {
  W5 identity;
  identity.top = {0,1,2,3,4};
  identity.components.fill(core.s3_index.at(SmallPerm{0,1,2}));
  std::vector<W5> elements = {identity};
  std::unordered_map<std::uint64_t, int> seen;
  seen.reserve(static_cast<std::size_t>(item.order) * 5 / 4);
  seen.emplace(key5(identity, item.m), 0);
  for (std::size_t head = 0; head < elements.size(); ++head) {
    for (const auto& generator : item.generators) {
      const W5 product = multiply5(elements[head], generator, core, item.m);
      const auto [it, inserted] = seen.emplace(key5(product, item.m), elements.size());
      if (inserted) elements.push_back(product);
    }
  }
  if (elements.size() != static_cast<std::size_t>(item.order)) {
    throw std::runtime_error("quotient closure order mismatch");
  }
  return elements;
}

Tuple5 act5(
    const Tuple5& input, const W5& element, const CoreData& core, int width) {
  Tuple5 output{};
  for (int i = 0; i < width; ++i) {
    output[element.top[i]] = core.x_action[element.components[i]][input[i]];
  }
  return output;
}

constexpr std::array<std::size_t, 5> powers = {
    1, 55, 55*55, 55*55*55, 55ULL*55*55*55};

std::vector<std::uint64_t> regular_prefix_masks(
    const Case5& item, const CoreData& core, const std::vector<W5>& group) {
  std::vector<std::uint64_t> masks(powers[item.m - 1], 0);
  for (const auto& representative : item.x_reps) {
    for (const auto& element : group) {
      const Tuple5 image = act5(representative, element, core, item.m);
      std::size_t prefix = 0;
      for (int i = 0; i < item.m - 1; ++i) prefix += image[i] * powers[i];
      masks[prefix] |= std::uint64_t{1} << image[item.m - 1];
    }
  }
  std::uint64_t count = 0;
  for (auto mask : masks) count += std::popcount(mask);
  if (count != item.regular_points ||
      count != item.x_reps.size() * static_cast<std::uint64_t>(item.order)) {
    throw std::runtime_error("regular-set census mismatch: got " +
        std::to_string(count) + " expected " +
        std::to_string(item.regular_points));
  }
  return masks;
}

bool box_search_recursive(
    int depth, int width, const std::array<int,4>& order,
    const std::array<std::uint64_t,5>& allowed,
    std::size_t prefix, const std::vector<std::uint64_t>& regular,
    Tuple5* witness, std::uint64_t* probes) {
  if (depth == width - 1) {
    ++*probes;
    const std::uint64_t possible = regular[prefix] & allowed[width - 1];
    if (!possible) return false;
    witness->at(width - 1) = std::countr_zero(possible);
    return true;
  }
  const int coordinate = order[depth];
  std::uint64_t values = allowed[coordinate];
  while (values) {
    const int value = std::countr_zero(values);
    values &= values - 1;
    witness->at(coordinate) = value;
    if (box_search_recursive(depth + 1, width, order, allowed,
                             prefix + value * powers[coordinate], regular,
                             witness, probes)) return true;
  }
  return false;
}

bool target_witness(
    const Tuple5& target, const Case5& item, const CoreData& core,
    const std::vector<W5>& group, const std::vector<std::uint64_t>& regular,
    Tuple5* left, Tuple5* right, std::size_t* left_index,
    std::uint64_t* probes, std::uint64_t* expanded_left_tests) {
  auto try_left = [&](const Tuple5& u, std::size_t index) {
    std::array<std::uint64_t,5> allowed{};
    for (int i = 0; i < item.m; ++i) allowed[i] = core.incidence[target[i]][u[i]];
    std::array<int,4> order = {0,1,2,3};
    std::sort(order.begin(), order.begin() + item.m - 1, [&](int a, int b) {
      return std::popcount(allowed[a]) < std::popcount(allowed[b]);
    });
    Tuple5 v{};
    if (box_search_recursive(0, item.m, order, allowed, 0, regular, &v, probes)) {
      *left = u;
      *right = v;
      *left_index = index;
      return true;
    }
    return false;
  };
  for (std::size_t index = 0; index < item.x_reps.size(); ++index) {
    if (try_left(item.x_reps[index], index)) return true;
  }
  // A target representative need not be aligned with the selected regular
  // source-orbit representative.  The fallback is exhaustive: act every
  // source representative by every quotient element.  It is normally entered
  // only for a tiny number of target orbits.
  for (std::size_t index = 0; index < item.x_reps.size(); ++index) {
    for (std::size_t element = 1; element < group.size(); ++element) {
      ++*expanded_left_tests;
      if (try_left(act5(item.x_reps[index], group[element], core, item.m), index)) {
        return true;
      }
    }
  }
  return false;
}

std::uint64_t mix_witness(
    std::uint64_t hash, const Tuple5& target, const Tuple5& left,
    const Tuple5& right, std::size_t left_index, int width) {
  auto mix = [&](std::uint64_t value) {
    hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
  };
  mix(left_index);
  for (int i = 0; i < width; ++i) { mix(target[i]); mix(left[i]); mix(right[i]); }
  return hash;
}

void audit5(const std::string& mode, const Case5& item, const CoreData& core) {
  const auto group = closure5(item, core);
  const auto regular = regular_prefix_masks(item, core, group);
  std::uint64_t probes = 0;
  std::uint64_t expanded_left_tests = 0;
  std::uint64_t digest = 0xcbf29ce484222325ULL;
  std::size_t maximum_left_index = 0;
  for (std::size_t target_number = 0;
       target_number < item.y_reps.size(); ++target_number) {
    Tuple5 left{}, right{};
    std::size_t left_index = 0;
    if (!target_witness(item.y_reps[target_number], item, core, group, regular,
                        &left, &right, &left_index, &probes,
                        &expanded_left_tests)) {
      std::cout << "CORRELATED_M5_DEFECT|mode=" << mode
                << "|case=" << item.number
                << "|target_orbit=" << target_number + 1 << "\n";
      throw std::runtime_error("target orbit has no certified witness");
    }
    maximum_left_index = std::max(maximum_left_index, left_index);
    digest = mix_witness(
        digest, item.y_reps[target_number], left, right, left_index, item.m);
  }
  std::cout << "CORRELATED_M5_RESULT|mode=" << mode
            << "|case=" << item.number
            << "|m=" << item.m
            << "|quotient_order=" << item.order
            << "|regular_points=" << item.regular_points
            << "|regular_orbits=" << item.x_reps.size()
            << "|target_orbits=" << item.y_reps.size()
            << "|maximum_left_rep=" << maximum_left_index + 1
            << "|box_probes=" << probes
            << "|expanded_left_tests=" << expanded_left_tests
            << "|witness_digest=" << digest
            << "|status=PASS\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 3) {
      throw std::runtime_error("usage: engine FULL_TRANSCRIPT S3_TRANSCRIPT");
    }
    std::cout << "ENGINE|name=a5_k3_cd_correlated_m5_exact|schema=1\n";
    const CoreData core = build_core_data();
    const auto full = std::string(argv[1]) == "-"
        ? std::vector<Case5>{} : parse_representatives(argv[1], core);
    const auto s3 = std::string(argv[2]) == "-"
        ? std::vector<Case5>{} : parse_representatives(argv[2], core);
    if (full.empty() && s3.empty()) {
      throw std::runtime_error("unexpected residual-case census");
    }
    for (const auto& item : full) audit5("full", item, core);
    for (const auto& item : s3) audit5("s3", item, core);
    std::cout << "CORRELATED_M5_COMPLETE|cases=" << full.size() + s3.size()
              << "|failures=0|status=ALL_PASS\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ERROR|" << error.what() << '\n';
    return 1;
  }
}

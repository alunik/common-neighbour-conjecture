// Exact common-neighbour check for the maximal correlated quotient classes
// below (Out(A5) x S3) wr S3 that are genuinely base-two.
//
// Input is the immutable transcript produced by cd_maximal_exact_density.m.
// The 3600-point component is reduced exactly to 55 regular inner-diagonal
// colours and 77 target orbit types.  Quotient orbits are then enumerated on
// 55^3 and 77^3; no large primitive permutation action is constructed.

#define main diagonal_starstar_small_embedded_main
#include "../shared/diagonal_starstar_small.cpp"
#undef main
#include "../shared/coupled_engine.hpp"

#include <array>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <unordered_map>

namespace {

using SmallPerm = std::array<std::uint8_t, 3>;

struct RElement {
  std::uint8_t outer = 0;
  std::uint8_t s3 = 0;
};

struct WElement {
  std::array<std::uint8_t, 3> components{};
  std::uint8_t top = 0;
};

struct QuotientCase {
  int number = 0;
  int order = 0;
  std::uint64_t regular_points = 0;
  std::vector<WElement> generators;
};

struct CoreData {
  int x_size = 0;
  int y_size = 0;
  std::vector<std::vector<std::uint64_t>> incidence;
  std::array<std::vector<std::uint8_t>, 12> x_action;
  std::array<std::vector<std::uint8_t>, 12> y_action;
  std::array<SmallPerm, 6> s3_permutations{};
  std::map<SmallPerm, int> s3_index;
};

SmallPerm compose_small(const SmallPerm& left, const SmallPerm& right) {
  SmallPerm answer{};
  for (int i = 0; i < 3; ++i) answer[i] = right[left[i]];
  return answer;
}

std::vector<int> parse_numbers(const std::string& text) {
  std::vector<int> answer;
  static const std::regex number("[0-9]+");
  for (auto it = std::sregex_iterator(text.begin(), text.end(), number);
       it != std::sregex_iterator(); ++it) {
    answer.push_back(std::stoi(it->str()));
  }
  return answer;
}

CoreData build_core_data() {
  const FiniteGroup group = alternating_group(5);
  verify_group(group);
  const StateSpace space = make_state_space(group.order, 3);
  const ActionMaps maps = build_action_maps(group, space, 3);
  std::vector<const std::vector<Code>*> inner_actions;
  for (std::size_t i = 0; i < group.generators.size(); ++i) {
    inner_actions.push_back(&maps.base[i]);
  }
  const CoupledOrbitData orbit =
      coupled_stabiliser_orbits(space, inner_actions, group.order);
  std::vector<int> regular_oids;
  std::vector<int> regular_id(space.count, -1);
  for (int oid = 0; oid < static_cast<int>(orbit.points.size()); ++oid) {
    if (orbit.points[oid].size() != static_cast<std::size_t>(group.order)) {
      continue;
    }
    const int colour = static_cast<int>(regular_oids.size());
    regular_oids.push_back(oid);
    for (Code point : orbit.points[oid]) regular_id[point] = colour;
  }
  if (regular_oids.size() != 55 || orbit.points.size() != 77) {
    throw std::runtime_error("A5 core census mismatch");
  }

  CoreData data;
  data.x_size = static_cast<int>(regular_oids.size());
  data.y_size = static_cast<int>(orbit.points.size());
  data.incidence.assign(data.y_size,
                        std::vector<std::uint64_t>(data.x_size, 0));
  for (int y = 0; y < data.y_size; ++y) {
    const Code target = orbit.points[y][0];
    Digits inverse_digits{};
    for (int coordinate = 0; coordinate < space.coordinates; ++coordinate) {
      inverse_digits[coordinate] =
          group.inverse[space.digits[target][coordinate]];
    }
    for (int left = 0; left < data.x_size; ++left) {
      for (Code lambda : orbit.points[regular_oids[left]]) {
        Digits quotient_digits{};
        for (int coordinate = 0; coordinate < space.coordinates; ++coordinate) {
          quotient_digits[coordinate] = group.multiply(
              space.digits[lambda][coordinate], inverse_digits[coordinate]);
        }
        const int right = regular_id[encode_state(space, quotient_digits)];
        if (right >= 0) data.incidence[y][left] |= std::uint64_t{1} << right;
      }
    }
  }

  auto induced = [&](const std::vector<Code>& point_map,
                     bool regular) {
    const int count = regular ? data.x_size : data.y_size;
    std::vector<std::uint8_t> answer(count);
    for (int i = 0; i < count; ++i) {
      const int oid = regular ? regular_oids[i] : i;
      const int image_oid = orbit.id[point_map[orbit.points[oid][0]]];
      if (regular) {
        const auto found =
            std::find(regular_oids.begin(), regular_oids.end(), image_oid);
        if (found == regular_oids.end()) {
          throw std::runtime_error("R did not preserve regular core orbits");
        }
        answer[i] = found - regular_oids.begin();
      } else {
        answer[i] = image_oid;
      }
    }
    return answer;
  };
  const auto x_c3 = induced(maps.base[group.generators.size()], true);
  const auto x_t2 = induced(maps.symmetric, true);
  const auto x_outer = induced(maps.outer, true);
  const auto y_c3 = induced(maps.base[group.generators.size()], false);
  const auto y_t2 = induced(maps.symmetric, false);
  const auto y_outer = induced(maps.outer, false);

  struct S3Action {
    SmallPerm permutation;
    std::vector<std::uint8_t> x;
    std::vector<std::uint8_t> y;
  };
  std::vector<std::uint8_t> x_identity(data.x_size);
  std::vector<std::uint8_t> y_identity(data.y_size);
  std::iota(x_identity.begin(), x_identity.end(), 0);
  std::iota(y_identity.begin(), y_identity.end(), 0);
  std::vector<S3Action> s3 = {{{0,1,2}, x_identity, y_identity}};
  const std::vector<S3Action> generators = {
      {{1,2,0}, x_c3, y_c3}, {{1,0,2}, x_t2, y_t2}};
  std::set<SmallPerm> seen = {{0,1,2}};
  auto compose_map = [](const std::vector<std::uint8_t>& left,
                        const std::vector<std::uint8_t>& right) {
    std::vector<std::uint8_t> answer(left.size());
    for (std::size_t i = 0; i < left.size(); ++i) {
      answer[i] = right[left[i]];
    }
    return answer;
  };
  for (std::size_t head = 0; head < s3.size(); ++head) {
    for (const auto& generator : generators) {
      const SmallPerm permutation =
          compose_small(s3[head].permutation, generator.permutation);
      if (!seen.insert(permutation).second) continue;
      s3.push_back({permutation,
                    compose_map(s3[head].x, generator.x),
                    compose_map(s3[head].y, generator.y)});
    }
  }
  if (s3.size() != 6) throw std::runtime_error("S3 closure mismatch");
  for (int i = 0; i < 6; ++i) {
    data.s3_permutations[i] = s3[i].permutation;
    data.s3_index.emplace(s3[i].permutation, i);
    for (int outer = 0; outer < 2; ++outer) {
      const int rid = outer * 6 + i;
      data.x_action[rid] = outer
          ? compose_map(s3[i].x, x_outer) : s3[i].x;
      data.y_action[rid] = outer
          ? compose_map(s3[i].y, y_outer) : s3[i].y;
    }
  }

  // Exact covariance of every incidence edge under the three R generators.
  for (int rid : {1, data.s3_index.at({1,0,2}), 6}) {
    for (int y = 0; y < data.y_size; ++y) {
      for (int left = 0; left < data.x_size; ++left) {
        for (int right = 0; right < data.x_size; ++right) {
          const bool edge = (data.incidence[y][left] >> right) & 1U;
          const bool image_edge =
              (data.incidence[data.y_action[rid][y]]
                             [data.x_action[rid][left]] >>
               data.x_action[rid][right]) & 1U;
          if (edge != image_edge) {
            throw std::runtime_error("core incidence covariance failed");
          }
        }
      }
    }
  }
  return data;
}

std::vector<QuotientCase> parse_cases(
    const std::string& path, const CoreData& core) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open density transcript");
  std::map<int, QuotientCase> cases;
  std::string line;
  const std::regex case_re(
      R"re(^EXACT_DENSITY_CASE\|case=([0-9]+)\|M=([0-9]+).*\|regular=([0-9]+).*\|base_two=true\|)re");
  const std::regex generator_re(
      R"re(^EXACT_DENSITY_GENERATOR\|case=([0-9]+)\|generator=([0-9]+)\|top=\[([^\]]+)\]\|components=\[(.*)\]$)re");
  const std::regex frontier_case_re(
      R"re(^FRONTIER_BASE2\|node=([0-9]+).*\|order=([0-9]+)\|.*$)re");
  const std::regex frontier_generator_re(
      R"re(^FRONTIER_GENERATOR\|node=([0-9]+)\|generator=([0-9]+)\|top=\[([^\]]+)\]\|components=\[(.*)\]$)re");
  while (std::getline(input, line)) {
    std::smatch match;
    if (std::regex_search(line, match, case_re)) {
      QuotientCase item;
      item.number = std::stoi(match[1]);
      item.order = std::stoi(match[2]);
      item.regular_points = std::stoull(match[3]);
      cases[item.number] = item;
      continue;
    }
    if (std::regex_match(line, match, frontier_case_re)) {
      QuotientCase item;
      item.number = std::stoi(match[1]);
      item.order = std::stoi(match[2]);
      item.regular_points = 0;
      cases[item.number] = item;
      continue;
    }
    if (!std::regex_match(line, match, generator_re) &&
        !std::regex_match(line, match, frontier_generator_re)) continue;
    const int number = std::stoi(match[1]);
    if (!cases.contains(number)) continue;
    const auto top_values = parse_numbers(match[3]);
    if (top_values.size() != 3) throw std::runtime_error("bad top descriptor");
    SmallPerm top{};
    for (int i = 0; i < 3; ++i) top[i] = top_values[i] - 1;
    WElement generator;
    generator.top = core.s3_index.at(top);

    std::stringstream component_stream(match[4]);
    std::string component;
    int coordinate = 0;
    while (std::getline(component_stream, component, ';')) {
      const auto values = parse_numbers(component);
      if (values.size() != 4 || coordinate >= 3) {
        throw std::runtime_error("bad component descriptor");
      }
      SmallPerm s3{};
      for (int i = 0; i < 3; ++i) s3[i] = values[i + 1] - 1;
      generator.components[coordinate++] =
          values[0] * 6 + core.s3_index.at(s3);
    }
    if (coordinate != 3) throw std::runtime_error("wrong component count");
    cases[number].generators.push_back(generator);
  }
  std::vector<QuotientCase> answer;
  for (auto& [number, item] : cases) {
    if (item.generators.empty()) {
      throw std::runtime_error("base-two case has no generators");
    }
    answer.push_back(std::move(item));
  }
  if (answer.empty()) throw std::runtime_error("no base-two quotient cases");
  return answer;
}

int s3_product(int left, int right, const CoreData& core) {
  const SmallPerm product = compose_small(
      core.s3_permutations[left], core.s3_permutations[right]);
  return core.s3_index.at(product);
}

std::uint8_t r_product(
    std::uint8_t left, std::uint8_t right, const CoreData& core) {
  return ((left / 6) ^ (right / 6)) * 6 +
         s3_product(left % 6, right % 6, core);
}

WElement wreath_product(
    const WElement& left, const WElement& right, const CoreData& core) {
  WElement answer;
  answer.top = s3_product(left.top, right.top, core);
  const auto& left_top = core.s3_permutations[left.top];
  for (int i = 0; i < 3; ++i) {
    answer.components[i] = r_product(
        left.components[i], right.components[left_top[i]], core);
  }
  return answer;
}

std::uint32_t wreath_key(const WElement& element) {
  std::uint32_t key = element.top;
  for (int i = 0; i < 3; ++i) key = key * 12 + element.components[i];
  return key;
}

std::vector<WElement> group_closure(
    const QuotientCase& item, const CoreData& core) {
  WElement identity;
  identity.top = core.s3_index.at({0,1,2});
  identity.components.fill(0 * 6 + identity.top);
  std::vector<WElement> elements = {identity};
  std::unordered_map<std::uint32_t, int> seen = {{wreath_key(identity), 0}};
  for (std::size_t head = 0; head < elements.size(); ++head) {
    for (const auto& generator : item.generators) {
      const WElement product =
          wreath_product(elements[head], generator, core);
      const auto [it, inserted] =
          seen.emplace(wreath_key(product), elements.size());
      if (inserted) elements.push_back(product);
    }
  }
  if (static_cast<int>(elements.size()) != item.order) {
    throw std::runtime_error("quotient closure order mismatch");
  }
  return elements;
}

std::uint32_t act_tuple(
    std::uint32_t code, int alphabet, const WElement& element,
    const CoreData& core, bool regular_alphabet) {
  std::array<int, 3> input{};
  for (int i = 0; i < 3; ++i) {
    input[i] = code % alphabet;
    code /= alphabet;
  }
  std::array<int, 3> output{};
  const auto& top = core.s3_permutations[element.top];
  for (int i = 0; i < 3; ++i) {
    const auto& action = regular_alphabet
        ? core.x_action[element.components[i]]
        : core.y_action[element.components[i]];
    output[top[i]] = action[input[i]];
  }
  return output[0] + alphabet * (output[1] + alphabet * output[2]);
}

struct OrbitCensus {
  std::vector<std::uint8_t> regular;
  std::vector<std::uint32_t> representatives;
  int regular_orbits = 0;
  std::uint64_t regular_points = 0;
};

OrbitCensus orbit_census(
    int alphabet, const QuotientCase& item, const CoreData& core,
    bool regular_alphabet) {
  const std::uint32_t degree = alphabet * alphabet * alphabet;
  std::vector<std::int32_t> orbit_id(degree, -1);
  OrbitCensus census;
  if (regular_alphabet) census.regular.assign(degree, 0);
  std::vector<std::uint32_t> queue;
  for (std::uint32_t start = 0; start < degree; ++start) {
    if (orbit_id[start] >= 0) continue;
    const int oid = static_cast<int>(census.representatives.size());
    census.representatives.push_back(start);
    queue.clear();
    queue.push_back(start);
    orbit_id[start] = oid;
    for (std::size_t head = 0; head < queue.size(); ++head) {
      for (const auto& generator : item.generators) {
        const auto image = act_tuple(
            queue[head], alphabet, generator, core, regular_alphabet);
        if (orbit_id[image] >= 0) continue;
        orbit_id[image] = oid;
        queue.push_back(image);
      }
    }
    if (regular_alphabet && static_cast<int>(queue.size()) == item.order) {
      ++census.regular_orbits;
      census.regular_points += queue.size();
      for (auto point : queue) census.regular[point] = 1;
    }
  }
  return census;
}

std::array<int, 3> decode(std::uint32_t code, int alphabet) {
  std::array<int, 3> answer{};
  for (int i = 0; i < 3; ++i) {
    answer[i] = code % alphabet;
    code /= alphabet;
  }
  return answer;
}

bool target_has_common_neighbour(
    const std::array<int, 3>& target, const CoreData& core,
    const std::vector<std::uint8_t>& regular,
    const std::vector<std::uint64_t>& regular_last_masks) {
  const int x = core.x_size;
  for (std::uint32_t u_code = 0; u_code < regular.size(); ++u_code) {
    if (!regular[u_code]) continue;
    const auto u = decode(u_code, x);
    const std::uint64_t a0 = core.incidence[target[0]][u[0]];
    const std::uint64_t a1 = core.incidence[target[1]][u[1]];
    const std::uint64_t a2 = core.incidence[target[2]][u[2]];
    for (int v0 = 0; v0 < x; ++v0) {
      if (((a0 >> v0) & 1U) == 0) continue;
      for (int v1 = 0; v1 < x; ++v1) {
        if (((a1 >> v1) & 1U) == 0) continue;
        if (regular_last_masks[v0 + x * v1] & a2) return true;
      }
    }
  }
  return false;
}

bool audit_case(const QuotientCase& item, const CoreData& core) {
  const auto elements = group_closure(item, core);
  const OrbitCensus x_census =
      orbit_census(core.x_size, item, core, true);
  if (item.regular_points != 0 &&
      x_census.regular_points != item.regular_points) {
    throw std::runtime_error("regular-point count disagrees with Magma");
  }
  const OrbitCensus y_census =
      orbit_census(core.y_size, item, core, false);
  std::vector<std::uint64_t> regular_last_masks(
      core.x_size * core.x_size, 0);
  for (std::uint32_t code = 0; code < x_census.regular.size(); ++code) {
    if (!x_census.regular[code]) continue;
    const auto tuple = decode(code, core.x_size);
    regular_last_masks[tuple[0] + core.x_size * tuple[1]] |=
        std::uint64_t{1} << tuple[2];
  }
  std::uint64_t checked = 0;
  for (auto representative : y_census.representatives) {
    ++checked;
    const auto target = decode(representative, core.y_size);
    if (!target_has_common_neighbour(
            target, core, x_census.regular, regular_last_masks)) {
      std::cout << "CORRELATED_DEFECT|case=" << item.number
                << "|target_code=" << representative
                << "|target=[" << target[0] + 1 << ','
                << target[1] + 1 << ',' << target[2] + 1 << "]\n";
      return false;
    }
  }
  std::cout << "CORRELATED_RESULT|case=" << item.number
            << "|quotient_order=" << item.order
            << "|group_elements=" << elements.size()
            << "|regular_points=" << x_census.regular_points
            << "|regular_orbits=" << x_census.regular_orbits
            << "|target_orbits=" << y_census.representatives.size()
            << "|targets_checked=" << checked
            << "|status=PASS\n";
  return true;
}

}  // namespace

#ifndef A5_CORRELATED_M3_EMBEDDED
int main(int argc, char** argv) {
  try {
    if (argc != 2) {
      std::cerr << "usage: a5_k3_cd_correlated_m3_exact TRANSCRIPT\n";
      return 1;
    }
    std::cout << "ENGINE|name=a5_k3_cd_correlated_m3_exact|schema=1\n";
    const CoreData core = build_core_data();
    const auto cases = parse_cases(argv[1], core);
    int failures = 0;
    for (const auto& item : cases) failures += !audit_case(item, core);
    std::cout << "CORRELATED_RUN_RESULT|cases=" << cases.size()
              << "|failures=" << failures
              << "|status=" << (failures == 0 ? "ALL_PASS" : "FAIL")
              << '\n';
    return failures == 0 ? 0 : 2;
  } catch (const std::exception& error) {
    std::cerr << "ERROR|" << error.what() << '\n';
    return 1;
  }
}
#endif

// Exact common-neighbour lift for full CD wreath products built from each
// A5,k=3 diagonal component.  For a component point stabiliser H, every
// target H-orbit gives an r by r incidence matrix on the regular H-orbits.
// A common neighbour in L wr P is exactly a choice of one incidence edge in
// each coordinate whose two endpoint colourings both distinguish P.

#define main diagonal_starstar_small_embedded_main
#include "../shared/diagonal_starstar_small.cpp"
#undef main
#include "../shared/coupled_engine.hpp"

#include <map>
#include <set>
#include <unordered_map>

namespace {

using Mask = std::uint32_t;

struct ComponentSpec {
  std::string label;
  int quotient_order;
  bool symmetric;
  bool outer;
  bool sign_graph;
};

struct ComponentData {
  std::string label;
  int r = 0;
  int target_orbits = 0;
  std::vector<std::vector<Mask>> matrices;
  bool all_perfect = false;
  std::uint64_t matrix_hash = 1469598103934665603ULL;
};

struct TopGroup {
  std::string label;
  int degree;
  int distinguishing_number;
  std::vector<std::vector<int>> elements;
};

std::vector<int> compose_permutations(
    const std::vector<int>& left, const std::vector<int>& right) {
  std::vector<int> answer(left.size());
  for (std::size_t point = 0; point < left.size(); ++point) {
    answer[point] = right[left[point]];
  }
  return answer;
}

TopGroup make_top(
    const std::string& label, int degree, int distinguishing_number,
    const std::vector<std::vector<int>>& generators,
    int expected_order) {
  std::set<std::vector<int>> seen;
  std::vector<std::vector<int>> queue;
  std::vector<int> identity(degree);
  std::iota(identity.begin(), identity.end(), 0);
  seen.insert(identity);
  queue.push_back(identity);
  for (std::size_t head = 0; head < queue.size(); ++head) {
    for (const auto& generator : generators) {
      const auto product = compose_permutations(queue[head], generator);
      if (seen.insert(product).second) queue.push_back(product);
    }
  }
  if (static_cast<int>(queue.size()) != expected_order) {
    throw std::runtime_error(label + " top-group order mismatch");
  }
  return {label, degree, distinguishing_number, queue};
}

std::vector<TopGroup> top_groups() {
  auto cycle = [](int degree) {
    std::vector<int> p(degree);
    for (int i = 0; i < degree; ++i) p[i] = (i + 1) % degree;
    return p;
  };
  auto transposition = [](int degree) {
    std::vector<int> p(degree);
    std::iota(p.begin(), p.end(), 0);
    std::swap(p[0], p[1]);
    return p;
  };
  return {
      make_top("C3", 3, 2, {cycle(3)}, 3),
      make_top("S3", 3, 3, {cycle(3), transposition(3)}, 6),
      make_top("C4", 4, 2, {cycle(4)}, 4),
      make_top("V4", 4, 2, {{1,0,3,2}, {2,3,0,1}}, 4),
      make_top("D8", 4, 3, {cycle(4), {0,3,2,1}}, 8),
      make_top("A4", 4, 3, {{1,2,0,3}, {1,0,3,2}}, 12),
      make_top("S4", 4, 4, {cycle(4), transposition(4)}, 24),
      make_top("C5", 5, 2, {cycle(5)}, 5),
      make_top("D10", 5, 3, {cycle(5), {0,4,3,2,1}}, 10),
      make_top("F20", 5, 3, {cycle(5), {0,2,4,1,3}}, 20),
      make_top("A5", 5, 4, {cycle(5), {1,2,0,3,4}}, 60),
      make_top("S5", 5, 5, {cycle(5), transposition(5)}, 120),
  };
}

std::vector<const std::vector<Code>*> component_actions(
    const ActionMaps& maps, const ComponentSpec& spec,
    std::vector<Code>& combined) {
  std::vector<const std::vector<Code>*> actions;
  for (const auto& action : maps.base) actions.push_back(&action);
  if (spec.sign_graph) {
    if (!coupled_maps_commute(maps.symmetric, maps.outer)) {
      throw std::runtime_error("outer and odd-top actions do not commute");
    }
    combined = coupled_compose_maps(maps.symmetric, maps.outer);
    actions.push_back(&combined);
  } else {
    if (spec.symmetric) actions.push_back(&maps.symmetric);
    if (spec.outer) actions.push_back(&maps.outer);
  }
  return actions;
}

ComponentData build_component_data(
    const FiniteGroup& group, const StateSpace& space,
    const ActionMaps& maps, const ComponentSpec& spec) {
  std::vector<Code> combined;
  const auto actions = component_actions(maps, spec, combined);
  const int stabiliser_order = group.order * spec.quotient_order;
  const CoupledOrbitData orbit =
      coupled_stabiliser_orbits(space, actions, stabiliser_order);
  std::vector<int> regular_oids;
  std::vector<int> lambda_of_point(space.count, -1);
  for (int oid = 0; oid < static_cast<int>(orbit.points.size()); ++oid) {
    if (orbit.points[oid].size() !=
        static_cast<std::size_t>(stabiliser_order)) continue;
    const int lambda = static_cast<int>(regular_oids.size());
    regular_oids.push_back(oid);
    for (Code point : orbit.points[oid]) lambda_of_point[point] = lambda;
  }

  ComponentData data;
  data.label = spec.label;
  data.r = static_cast<int>(regular_oids.size());
  data.target_orbits = static_cast<int>(orbit.points.size());
  std::map<std::vector<Mask>, int> unique;
  for (const auto& target_orbit : orbit.points) {
    const Code target = target_orbit[0];
    Digits inverse_digits{};
    for (int coordinate = 0; coordinate < space.coordinates; ++coordinate) {
      inverse_digits[coordinate] =
          group.inverse[space.digits[target][coordinate]];
    }
    std::vector<Mask> rows(data.r, 0);
    for (int left = 0; left < data.r; ++left) {
      for (Code lambda : orbit.points[regular_oids[left]]) {
        Digits quotient_digits{};
        for (int coordinate = 0; coordinate < space.coordinates; ++coordinate) {
          quotient_digits[coordinate] = group.multiply(
              space.digits[lambda][coordinate], inverse_digits[coordinate]);
        }
        const int right =
            lambda_of_point[encode_state(space, quotient_digits)];
        if (right >= 0) rows[left] |= Mask{1} << right;
      }
    }
    if (!unique.contains(rows)) {
      unique.emplace(rows, static_cast<int>(data.matrices.size()));
      data.matrices.push_back(rows);
    }
  }
  data.all_perfect = true;
  for (const auto& rows : data.matrices) {
    std::vector<std::vector<std::uint8_t>> adjacency(
        data.r, std::vector<std::uint8_t>(data.r, 0));
    for (int i = 0; i < data.r; ++i) {
      for (int j = 0; j < data.r; ++j) {
        adjacency[i][j] = (rows[i] >> j) & 1U;
        data.matrix_hash ^= adjacency[i][j] + 31U * i + 131U * j;
        data.matrix_hash *= 1099511628211ULL;
      }
    }
    if (coupled_maximum_matching(adjacency).size != data.r) {
      data.all_perfect = false;
    }
  }
  std::cout << "COMPONENT_MATRIX_CENSUS|label=" << data.label
            << "|r=" << data.r
            << "|target_orbits=" << data.target_orbits
            << "|unique_matrices=" << data.matrices.size()
            << "|all_perfect=" << data.all_perfect
            << "|matrix_hash=" << data.matrix_hash << '\n';
  return data;
}

bool distinguishes(
    const std::vector<std::uint8_t>& colours, const TopGroup& top) {
  for (std::size_t index = 1; index < top.elements.size(); ++index) {
    const auto& permutation = top.elements[index];
    bool fixed = true;
    for (int point = 0; point < top.degree; ++point) {
      if (colours[point] != colours[permutation[point]]) {
        fixed = false;
        break;
      }
    }
    if (fixed) return false;
  }
  return true;
}

void generate_distinguishing_colourings(
    const TopGroup& top, int colours,
    std::vector<std::vector<std::uint8_t>>& output) {
  std::uint64_t count = 1;
  for (int i = 0; i < top.degree; ++i) count *= colours;
  for (std::uint64_t code = 0; code < count; ++code) {
    std::uint64_t value = code;
    std::vector<std::uint8_t> colouring(top.degree);
    for (int i = 0; i < top.degree; ++i) {
      colouring[i] = value % colours;
      value /= colours;
    }
    if (distinguishes(colouring, top)) output.push_back(colouring);
  }
}

bool canonical_profile(
    const std::vector<int>& profile, const TopGroup& top) {
  for (const auto& permutation : top.elements) {
    std::vector<int> image(top.degree);
    for (int point = 0; point < top.degree; ++point) {
      image[point] = profile[permutation[point]];
    }
    if (image < profile) return false;
  }
  return true;
}

bool right_colouring_exists(
    std::uint64_t key, const std::vector<Mask>& masks,
    const std::vector<std::vector<std::uint8_t>>& distinguishing,
    std::unordered_map<std::uint64_t, bool>& cache) {
  const auto found = cache.find(key);
  if (found != cache.end()) return found->second;
  for (const auto& right : distinguishing) {
    bool allowed = true;
    for (std::size_t i = 0; i < masks.size(); ++i) {
      if (((masks[i] >> right[i]) & 1U) == 0) {
        allowed = false;
        break;
      }
    }
    if (allowed) {
      cache.emplace(key, true);
      return true;
    }
  }
  cache.emplace(key, false);
  return false;
}

bool profile_has_common_neighbour(
    const ComponentData& component,
    const std::vector<int>& profile,
    const std::vector<std::vector<std::uint8_t>>& distinguishing,
    std::unordered_map<std::uint64_t, bool>& right_cache) {
  for (const auto& left : distinguishing) {
    std::vector<Mask> masks(profile.size());
    std::uint64_t key = 0;
    for (std::size_t i = 0; i < profile.size(); ++i) {
      masks[i] = component.matrices[profile[i]][left[i]];
      key |= static_cast<std::uint64_t>(masks[i]) << (component.r * i);
    }
    if (right_colouring_exists(key, masks, distinguishing, right_cache)) {
      return true;
    }
  }
  return false;
}

bool audit_wreath(const ComponentData& component, const TopGroup& top) {
  if (component.r < top.distinguishing_number) {
    std::cout << "WREATH_RESULT|component=" << component.label
              << "|top=" << top.label
              << "|m=" << top.degree
              << "|r=" << component.r
              << "|D=" << top.distinguishing_number
              << "|status=BASE_SIZE_GT_2\n";
    return true;
  }
  if (component.all_perfect && component.r >= 2 * top.degree - 1) {
    std::cout << "WREATH_RESULT|component=" << component.label
              << "|top=" << top.label
              << "|m=" << top.degree
              << "|r=" << component.r
              << "|D=" << top.distinguishing_number
              << "|status=PASS"
              << "|method=greedy_rainbow_perfect_matchings\n";
    return true;
  }

  std::vector<std::vector<std::uint8_t>> distinguishing;
  generate_distinguishing_colourings(top, component.r, distinguishing);
  if (distinguishing.empty()) {
    throw std::runtime_error("distinguishing-colouring census is empty");
  }
  const int types = static_cast<int>(component.matrices.size());
  std::uint64_t profiles = 1;
  for (int i = 0; i < top.degree; ++i) profiles *= types;
  std::uint64_t canonical = 0;
  std::unordered_map<std::uint64_t, bool> right_cache;
  for (std::uint64_t code = 0; code < profiles; ++code) {
    std::uint64_t value = code;
    std::vector<int> profile(top.degree);
    for (int i = 0; i < top.degree; ++i) {
      profile[i] = value % types;
      value /= types;
    }
    if (!canonical_profile(profile, top)) continue;
    ++canonical;
    if (!profile_has_common_neighbour(
            component, profile, distinguishing, right_cache)) {
      std::cout << "WREATH_DEFECT|component=" << component.label
                << "|top=" << top.label << "|profile=[";
      for (int i = 0; i < top.degree; ++i) {
        if (i) std::cout << ',';
        std::cout << profile[i] + 1;
      }
      std::cout << "]|certificate=no_distinguishing_incidence_selection\n";
      return false;
    }
  }
  std::cout << "WREATH_RESULT|component=" << component.label
            << "|top=" << top.label
            << "|m=" << top.degree
            << "|r=" << component.r
            << "|D=" << top.distinguishing_number
            << "|matrix_types=" << types
            << "|canonical_profiles=" << canonical
            << "|distinguishing_colourings=" << distinguishing.size()
            << "|right_cache=" << right_cache.size()
            << "|status=PASS|method=exact_incidence_csp\n";
  return true;
}

}  // namespace

int sealed_m5_main() {
  try {
    std::cout << "ENGINE|name=a5_k3_cd_wreath_exact|schema=1\n";
    coupled_matching_self_test();
    const FiniteGroup group = alternating_group(5);
    verify_group(group);
    const StateSpace space = make_state_space(group.order, 3);
    const ActionMaps maps = build_action_maps(group, space, 3);
    const std::vector<ComponentSpec> specs = {
        {"A5_k3_C3", 3, false, false, false},
        {"A5_k3_sign_graph", 6, false, false, true},
        {"A5_k3_C2xC3", 6, false, true, false},
        {"A5_k3_S3", 6, true, false, false},
        {"A5_k3_C2xS3", 12, true, true, false},
    };
    std::vector<ComponentData> components;
    for (const auto& spec : specs) {
      components.push_back(build_component_data(group, space, maps, spec));
    }
    const auto tops = top_groups();
    int cases = 0;
    int failures = 0;
    for (const auto& component : components) {
      for (const auto& top : tops) {
        ++cases;
        failures += !audit_wreath(component, top);
      }
    }
    std::cout << "WREATH_RUN_RESULT|cases=" << cases
              << "|failures=" << failures
              << "|status=" << (failures == 0 ? "ALL_PASS" : "FAIL")
              << '\n';
    return failures == 0 ? 0 : 2;
  } catch (const std::exception& error) {
    std::cerr << "ERROR|" << error.what() << '\n';
    return 1;
  }
}


// Narrow degree-six endpoint for the five A5,k=3 component actions.
int main() {
  try {
    std::cout << "ENGINE|name=a5_k3_cd_s6_exact|schema=1\n";
    coupled_matching_self_test();
    const FiniteGroup group = alternating_group(5);
    verify_group(group);
    const StateSpace space = make_state_space(group.order, 3);
    const ActionMaps maps = build_action_maps(group, space, 3);
    const std::vector<ComponentSpec> specs = {
        {"A5_k3_C3", 3, false, false, false},
        {"A5_k3_sign_graph", 6, false, false, true},
        {"A5_k3_C2xC3", 6, false, true, false},
        {"A5_k3_S3", 6, true, false, false},
        {"A5_k3_C2xS3", 12, true, true, false},
    };
    auto z = [](std::initializer_list<int> one_based) {
      std::vector<int> result;
      for (int value : one_based) result.push_back(value - 1);
      return result;
    };
    struct TopSpec {
      const char* label;
      int order;
      std::vector<std::vector<int>> generators;
    };
    // GAP's complete TransitiveGroup(6,i), i=1,...,16, census.  Auditing
    // every conjugacy type means that a base-two subgroup below the full
    // S6 envelope is covered without a recursive subgroup search.
    const std::vector<TopSpec> top_specs = {
      {"6T1",6,{z({2,3,4,5,6,1})}},
      {"6T2",6,{z({3,4,5,6,1,2}),z({4,3,2,1,6,5})}},
      {"6T3",12,{z({2,3,4,5,6,1}),z({4,3,2,1,6,5})}},
      {"6T4",12,{z({3,4,5,6,1,2}),z({4,5,3,1,2,6})}},
      {"6T5",18,{z({1,4,3,6,5,2}),z({4,5,6,1,2,3})}},
      {"6T6",24,{z({3,4,5,6,1,2}),z({1,2,6,4,5,3})}},
      {"6T7",24,{z({3,4,5,6,1,2}),z({4,5,3,1,2,6}),z({5,4,3,2,1,6})}},
      {"6T8",24,{z({3,4,5,6,1,2}),z({4,5,3,1,2,6}),z({5,4,6,2,1,3})}},
      {"6T9",36,{z({1,4,3,6,5,2}),z({4,5,6,1,2,3}),z({5,4,3,2,1,6})}},
      {"6T10",36,{z({4,1,6,5,2,3}),z({1,4,3,6,5,2})}},
      {"6T11",48,{z({3,4,5,6,1,2}),z({1,2,6,4,5,3}),z({5,4,3,2,1,6})}},
      {"6T12",60,{z({2,3,4,6,5,1}),z({4,2,3,1,6,5})}},
      {"6T13",72,{z({1,4,3,6,5,2}),z({1,4,3,2,5,6}),z({4,5,6,1,2,3})}},
      {"6T14",120,{z({2,3,4,6,5,1}),z({2,1,4,3,6,5})}},
      {"6T15",360,{z({2,3,4,5,1,6}),z({1,2,3,5,6,4})}},
      {"6T16",720,{z({2,3,4,5,6,1}),z({2,1,3,4,5,6})}},
    };
    std::vector<TopGroup> tops;
    for (const auto& spec : top_specs) {
      int distinguishing_number = 7;
      for (int colours = 1; colours <= 6; ++colours) {
        TopGroup candidate = make_top(spec.label, 6, colours,
                                      spec.generators, spec.order);
        std::vector<std::vector<std::uint8_t>> colourings;
        generate_distinguishing_colourings(candidate, colours, colourings);
        if (!colourings.empty()) {
          distinguishing_number = colours;
          break;
        }
      }
      tops.push_back(make_top(spec.label, 6, distinguishing_number,
                              spec.generators, spec.order));
      std::cout << "TOP_CENSUS|label=" << spec.label
                << "|order=" << spec.order
                << "|D=" << distinguishing_number << '\n';
    }
    int cases = 0;
    int failures = 0;
    for (const auto& spec : specs) {
      const ComponentData component =
          build_component_data(group, space, maps, spec);
      for (const auto& top : tops) {
        ++cases;
        failures += !audit_wreath(component, top);
      }
    }
    std::cout << "S6_RUN_RESULT|cases=" << cases
              << "|failures=" << failures
              << "|status=" << (failures == 0 ? "ALL_PASS" : "FAIL")
              << '\n';
    return failures == 0 ? 0 : 2;
  } catch (const std::exception& error) {
    std::cerr << "ERROR|" << error.what() << '\n';
    return 1;
  }
}

// Exact common-neighbour and stronger star-star checker for the five
// primitive simple-diagonal groups with socle PSL(2,7)^4.
//
// The action is on Diag(T)\T^4, represented by T^3.  The point stabiliser is
// H = T.Q, where Q is one of
//
//   A4, C2 x A4, S4, C2 x S4, graph(sign:S4 -> C2).
//
// We reuse the audited PSL(2,7) multiplication table and normalized diagonal
// action maps from research/diagonal/diagonal_starstar_small.cpp.  Only the
// case selection, exact quotient validation, concise orbit census and
// parallel target scan are new here.

#define main diagonal_starstar_small_embedded_main
#include "../../diagonal/shared/diagonal_starstar_small.cpp"
#undef main

#include <atomic>
#include <limits>
#include <set>
#include <sstream>
#include <tuple>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Permutation4 = std::array<std::uint8_t, 4>;

struct QuotientElement {
  bool outer = false;
  Permutation4 permutation{};

  bool operator<(const QuotientElement& other) const {
    return std::tie(outer, permutation) <
           std::tie(other.outer, other.permutation);
  }
};

Permutation4 identity_permutation4() {
  return Permutation4{0, 1, 2, 3};
}

Permutation4 to_permutation4(const std::vector<int>& input) {
  if (input.size() != 4) {
    throw std::runtime_error("top permutation does not have degree four");
  }
  Permutation4 answer{};
  std::array<std::uint8_t, 4> seen{};
  for (int point = 0; point < 4; ++point) {
    if (input[point] < 0 || input[point] >= 4 || seen[input[point]]) {
      throw std::runtime_error("invalid top permutation");
    }
    answer[point] = static_cast<std::uint8_t>(input[point]);
    seen[input[point]] = 1;
  }
  return answer;
}

Permutation4 compose_permutation4(
    const Permutation4& left, const Permutation4& right) {
  Permutation4 answer{};
  for (int point = 0; point < 4; ++point) {
    answer[point] = right[left[point]];
  }
  return answer;
}

QuotientElement quotient_multiply(
    const QuotientElement& left, const QuotientElement& right) {
  // Out(T) is central in Out(T) x S4.
  return QuotientElement{
      static_cast<bool>(left.outer != right.outer),
      compose_permutation4(left.permutation, right.permutation)};
}

std::set<QuotientElement> quotient_closure(
    const std::vector<QuotientElement>& generators) {
  const QuotientElement identity{false, identity_permutation4()};
  std::set<QuotientElement> closure{identity};
  std::vector<QuotientElement> queue{identity};
  for (std::size_t head = 0; head < queue.size(); ++head) {
    for (const QuotientElement& generator : generators) {
      const QuotientElement product =
          quotient_multiply(queue[head], generator);
      if (closure.insert(product).second) queue.push_back(product);
    }
  }
  return closure;
}

std::set<Permutation4> top_projection(
    const std::set<QuotientElement>& quotient) {
  std::set<Permutation4> answer;
  for (const auto& element : quotient) {
    answer.insert(element.permutation);
  }
  return answer;
}

std::uint8_t image_mask(
    std::uint8_t mask, const Permutation4& permutation) {
  std::uint8_t answer = 0;
  for (int point = 0; point < 4; ++point) {
    if (mask & (1U << point)) {
      answer |= static_cast<std::uint8_t>(1U << permutation[point]);
    }
  }
  return answer;
}

bool primitive_degree_four(const std::set<Permutation4>& group) {
  if (group.empty()) return false;
  std::uint8_t orbit = 0;
  for (const auto& permutation : group) {
    orbit |= static_cast<std::uint8_t>(1U << permutation[0]);
  }
  if (orbit != 0x0f) return false;

  // A nontrivial block may be translated to contain point zero.
  for (std::uint8_t block = 1; block < 0x0f; ++block) {
    if (!(block & 1U) || block == 1U) continue;
    bool is_block = true;
    for (const auto& permutation : group) {
      const std::uint8_t image = image_mask(block, permutation);
      if (image != block && (image & block) != 0) {
        is_block = false;
        break;
      }
    }
    if (is_block) return false;
  }
  return true;
}

enum class QuotientKind {
  A4,
  C2A4,
  S4,
  C2S4,
  SignGraphS4,
};

struct CaseSpec {
  int id = -1;
  QuotientKind kind = QuotientKind::A4;
  const char* label = nullptr;
  const char* quotient = nullptr;
  const char* top = nullptr;
  int quotient_order = 0;
  int stabiliser_order = 0;
  bool symmetric_top = false;
  bool independent_outer = false;
  bool sign_graph = false;
};

constexpr int kSimpleOrder = 168;
constexpr int kDiagonalLength = 4;
constexpr Code kExpectedDegree = 4741632;

const std::array<CaseSpec, 5> kCases = {{
    {0, QuotientKind::A4, "PSL27_k4_A4", "A4", "A4", 12, 2016,
     false, false, false},
    {1, QuotientKind::C2A4, "PSL27_k4_C2xA4", "C2xA4", "A4", 24,
     4032, false, true, false},
    {2, QuotientKind::S4, "PSL27_k4_S4", "S4", "S4", 24, 4032,
     true, false, false},
    {3, QuotientKind::C2S4, "PSL27_k4_C2xS4", "C2xS4", "S4", 48,
     8064, true, true, false},
    {4, QuotientKind::SignGraphS4, "PSL27_k4_sign_graph", "graph_sign_S4",
     "S4", 24, 4032, true, false, true},
}};

std::vector<QuotientElement> quotient_generators(const CaseSpec& spec) {
  std::vector<QuotientElement> generators;
  for (const auto& generator : alternating_top_generators(4)) {
    generators.push_back(
        QuotientElement{false, to_permutation4(generator)});
  }
  const Permutation4 odd =
      to_permutation4(symmetric_top_generator(4));
  if (spec.symmetric_top && !spec.sign_graph) {
    generators.push_back(QuotientElement{false, odd});
  }
  if (spec.independent_outer) {
    generators.push_back(
        QuotientElement{true, identity_permutation4()});
  }
  if (spec.sign_graph) {
    generators.push_back(QuotientElement{true, odd});
  }
  return generators;
}

void validate_case_metadata(const CaseSpec& spec) {
  const auto quotient = quotient_closure(quotient_generators(spec));
  const auto top = top_projection(quotient);
  const int expected_top_order = spec.symmetric_top ? 24 : 12;
  if (static_cast<int>(quotient.size()) != spec.quotient_order ||
      static_cast<int>(top.size()) != expected_top_order ||
      !primitive_degree_four(top) ||
      spec.stabiliser_order != kSimpleOrder * spec.quotient_order) {
    throw std::runtime_error(
        std::string("quotient metadata validation failed for ") +
        spec.label);
  }

  int outer_kernel = 0;
  for (const auto& element : quotient) {
    if (element.permutation == identity_permutation4()) ++outer_kernel;
  }
  const int expected_outer_kernel =
      spec.independent_outer ? 2 : 1;
  if (outer_kernel != expected_outer_kernel) {
    throw std::runtime_error(
        std::string("outer kernel validation failed for ") + spec.label);
  }
  if (spec.sign_graph) {
    for (const auto& element : quotient) {
      int inversions = 0;
      for (int first = 0; first < 4; ++first) {
        for (int second = first + 1; second < 4; ++second) {
          inversions +=
              element.permutation[first] > element.permutation[second];
        }
      }
      if (element.outer != static_cast<bool>(inversions % 2)) {
        throw std::runtime_error("sign graph parity validation failed");
      }
    }
  }

  std::cout << "CLASS_VALIDATION|case=" << spec.id
            << "|label=" << spec.label
            << "|T=PSL(2,7)|T_order=" << kSimpleOrder
            << "|k=" << kDiagonalLength
            << "|degree=" << kExpectedDegree
            << "|quotient=" << spec.quotient
            << "|quotient_order=" << quotient.size()
            << "|top=" << spec.top
            << "|top_order=" << top.size()
            << "|top_primitive=YES"
            << "|outer_kernel_order=" << outer_kernel
            << "|H_order=" << spec.stabiliser_order
            << "|status=PASS\n";
}

struct SdArguments {
  int selected_case = -1;
  int threads = 1;
  bool metadata_only = false;
  bool census_only = false;
  bool progress = false;
  bool verbose_targets = false;
};

SdArguments parse_sd_arguments(int argc, char** argv) {
  SdArguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--case") {
      if (++index >= argc) {
        throw std::runtime_error("--case requires an integer 0..4");
      }
      arguments.selected_case = std::stoi(argv[index]);
    } else if (argument == "--threads") {
      if (++index >= argc) {
        throw std::runtime_error("--threads requires an integer 1..32");
      }
      arguments.threads = std::stoi(argv[index]);
    } else if (argument == "--metadata") {
      arguments.metadata_only = true;
    } else if (argument == "--census") {
      arguments.census_only = true;
    } else if (argument == "--progress") {
      arguments.progress = true;
    } else if (argument == "--verbose-targets") {
      arguments.verbose_targets = true;
    } else if (argument == "--help") {
      std::cout
          << "usage: psl27_k4_starstar [--case 0..4] [--threads 1..32] "
             "[--metadata|--census] [--progress] [--verbose-targets]\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown argument: " + argument);
    }
  }
  if (arguments.selected_case < -1 || arguments.selected_case > 4) {
    throw std::runtime_error("--case must lie in 0..4");
  }
  if (arguments.threads < 1 || arguments.threads > 32) {
    throw std::runtime_error("--threads must lie in 1..32");
  }
  if (arguments.metadata_only && arguments.census_only) {
    throw std::runtime_error("--metadata and --census are mutually exclusive");
  }
  return arguments;
}

struct ExplicitOrbitData {
  std::vector<std::int32_t> id;
  std::vector<std::vector<Code>> points;
};

ExplicitOrbitData stabiliser_orbits_explicit(
    const StateSpace& space,
    const std::vector<const std::vector<Code>*>& actions,
    int stabiliser_order) {
  ExplicitOrbitData orbit;
  orbit.id.assign(space.count, -1);
  std::vector<Code> queue;
  for (Code start = 0; start < space.count; ++start) {
    if (orbit.id[start] >= 0) continue;
    const std::int32_t oid =
        static_cast<std::int32_t>(orbit.points.size());
    queue.clear();
    queue.push_back(start);
    orbit.id[start] = oid;
    for (std::size_t head = 0; head < queue.size(); ++head) {
      for (const auto* action : actions) {
        const Code target = (*action)[queue[head]];
        if (orbit.id[target] >= 0) continue;
        orbit.id[target] = oid;
        queue.push_back(target);
      }
    }
    if (queue.size() > static_cast<std::size_t>(stabiliser_order) ||
        stabiliser_order % static_cast<int>(queue.size()) != 0) {
      throw std::runtime_error(
          "H-orbit size does not divide asserted stabiliser order");
    }
    orbit.points.push_back(queue);
  }
  return orbit;
}

std::vector<Code> compose_maps(
    const std::vector<Code>& left, const std::vector<Code>& right) {
  if (left.size() != right.size()) {
    throw std::runtime_error("action-map degree mismatch");
  }
  std::vector<Code> answer(left.size());
  for (Code point = 0; point < left.size(); ++point) {
    answer[point] = right[left[point]];
  }
  return answer;
}

bool maps_commute(
    const std::vector<Code>& left, const std::vector<Code>& right) {
  if (left.size() != right.size()) return false;
  for (Code point = 0; point < left.size(); ++point) {
    if (left[right[point]] != right[left[point]]) return false;
  }
  return true;
}

struct TargetResult {
  std::size_t met = 0;
  std::uint64_t tests = 0;
  std::uint64_t reconstruction_checks = 0;
  Code first_witness = std::numeric_limits<Code>::max();
};

struct ScanResult {
  bool applicable = false;
  bool common_neighbour = false;
  bool starstar = false;
  int regular_orbits = 0;
  int h_orbits = 0;
  std::uint64_t regular_points = 0;
  std::size_t minimum_met = 0;
  int zero_common_targets = 0;
  std::uint64_t quotient_tests = 0;
  std::uint64_t reconstruction_checks = 0;
};

ScanResult run_case(
    const FiniteGroup& group,
    const StateSpace& space,
    const ActionMaps& maps,
    const CaseSpec& spec,
    const SdArguments& arguments) {
  std::vector<const std::vector<Code>*> actions;
  for (const auto& action : maps.base) actions.push_back(&action);
  std::vector<Code> combined;
  if (spec.sign_graph) {
    if (!maps_commute(maps.symmetric, maps.outer)) {
      throw std::runtime_error("outer and odd-top maps do not commute");
    }
    combined = compose_maps(maps.symmetric, maps.outer);
    Digits origin_digits{};
    origin_digits.fill(group.identity);
    verify_action_map(
        combined, encode_state(space, origin_digits),
        "PSL27_k4_combined_sign_outer");
    actions.push_back(&combined);
  } else {
    if (spec.symmetric_top) actions.push_back(&maps.symmetric);
    if (spec.independent_outer) actions.push_back(&maps.outer);
  }

  const ExplicitOrbitData orbit = stabiliser_orbits_explicit(
      space, actions, spec.stabiliser_order);
  std::vector<int> regular_oids;
  std::vector<std::uint8_t> regular(space.count, 0);
  std::uint64_t regular_points = 0;
  for (int oid = 0; oid < static_cast<int>(orbit.points.size()); ++oid) {
    if (orbit.points[oid].size() !=
        static_cast<std::size_t>(spec.stabiliser_order)) {
      continue;
    }
    regular_oids.push_back(oid);
    regular_points += orbit.points[oid].size();
    for (Code point : orbit.points[oid]) regular[point] = 1;
  }

  ScanResult result;
  result.regular_orbits = static_cast<int>(regular_oids.size());
  result.h_orbits = static_cast<int>(orbit.points.size());
  result.regular_points = regular_points;
  result.minimum_met = regular_oids.size();
  result.applicable = !regular_oids.empty();
  std::cout << "SD_CENSUS|case=" << spec.id
            << "|label=" << spec.label
            << "|degree=" << space.count
            << "|H_order=" << spec.stabiliser_order
            << "|H_orbits=" << orbit.points.size()
            << "|regular_orbits=" << regular_oids.size()
            << "|regular_points=" << regular_points
            << "|density=" << std::fixed << std::setprecision(9)
            << static_cast<double>(regular_points) / space.count
            << "|generated_H_order_certified_by_regular_orbit="
            << (regular_oids.empty() ? "NO" : "YES") << '\n';
  if (regular_oids.empty()) {
    std::cout << "SD_RESULT|case=" << spec.id
              << "|label=" << spec.label
              << "|status=NOT_APPLICABLE|reason=no_regular_orbit\n";
    return result;
  }

  std::uint64_t inverse_checks = 0;
  for (int oid : regular_oids) {
    for (Code point : orbit.points[oid]) {
      Digits inverse_digits{};
      for (int coordinate = 0; coordinate < space.coordinates; ++coordinate) {
        inverse_digits[coordinate] =
            group.inverse[space.digits[point][coordinate]];
      }
      if (!regular[encode_state(space, inverse_digits)]) {
        throw std::runtime_error("regular connection set is not inverse-stable");
      }
      ++inverse_checks;
    }
  }
  std::cout << "CONVENTION_VALIDATION|case=" << spec.id
            << "|coset_representative=(1,x1,x2,x3)"
            << "|neighbourhood=R*x"
            << "|membership=lambda*x^-1_in_R"
            << "|inverse_checks=" << inverse_checks
            << "|status=PASS\n";

  if (arguments.census_only) {
    std::cout << "SD_RESULT|case=" << spec.id
              << "|label=" << spec.label
              << "|status=CENSUS_PASS\n";
    return result;
  }

#ifdef _OPENMP
  omp_set_num_threads(arguments.threads);
#endif
  std::vector<TargetResult> target_results(orbit.points.size());
  std::atomic<int> completed{0};
  std::atomic<bool> reconstruction_failure{false};

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
  for (std::int64_t target_oid = 0;
       target_oid < static_cast<std::int64_t>(orbit.points.size());
       ++target_oid) {
    TargetResult target_result;
    const Code target = orbit.points[target_oid][0];
    Digits inverse_digits{};
    for (int coordinate = 0; coordinate < space.coordinates; ++coordinate) {
      inverse_digits[coordinate] =
          group.inverse[space.digits[target][coordinate]];
    }
    for (int oid : regular_oids) {
      bool found = false;
      for (Code lambda : orbit.points[oid]) {
        ++target_result.tests;
        Digits quotient_digits{};
        for (int coordinate = 0; coordinate < space.coordinates; ++coordinate) {
          quotient_digits[coordinate] = group.multiply(
              space.digits[lambda][coordinate], inverse_digits[coordinate]);
        }
        const Code quotient = encode_state(space, quotient_digits);
        if (!regular[quotient]) continue;
        for (int coordinate = 0; coordinate < space.coordinates; ++coordinate) {
          if (group.multiply(
                  quotient_digits[coordinate],
                  space.digits[target][coordinate]) !=
              space.digits[lambda][coordinate]) {
            reconstruction_failure.store(true);
          }
        }
        ++target_result.reconstruction_checks;
        if (target_result.first_witness ==
            std::numeric_limits<Code>::max()) {
          target_result.first_witness = lambda;
        }
        found = true;
        break;
      }
      target_result.met += found;
    }
    target_results[target_oid] = target_result;
    const int now = ++completed;
    if (arguments.progress && now % 25 == 0) {
#ifdef _OPENMP
#pragma omp critical(sd_progress_output)
#endif
      std::cerr << "SD_PROGRESS|case=" << spec.id
                << "|label=" << spec.label
                << "|targets=" << now << '/' << orbit.points.size()
                << "|threads=" << arguments.threads << '\n';
    }
  }
  if (reconstruction_failure.load()) {
    throw std::runtime_error("right-quotient reconstruction failed");
  }

  int minimum_target_oid = -1;
  result.minimum_met = regular_oids.size() + 1;
  for (int target_oid = 0;
       target_oid < static_cast<int>(target_results.size()); ++target_oid) {
    const TargetResult& target = target_results[target_oid];
    result.quotient_tests += target.tests;
    result.reconstruction_checks += target.reconstruction_checks;
    if (target.met < result.minimum_met) {
      result.minimum_met = target.met;
      minimum_target_oid = target_oid;
    }
    result.zero_common_targets += target.met == 0;
    if (arguments.verbose_targets || target.met == 0) {
      std::cout << "SD_TARGET|case=" << spec.id
                << "|H_orbit=" << target_oid + 1
                << "|target_code=" << orbit.points[target_oid][0]
                << "|met_regular_orbits=" << target.met
                << "|regular_orbits=" << regular_oids.size()
                << "|first_witness=";
      if (target.first_witness == std::numeric_limits<Code>::max()) {
        std::cout << "NONE";
      } else {
        std::cout << target.first_witness;
      }
      std::cout << "|quotient_tests=" << target.tests << '\n';
    }
  }
  result.common_neighbour = result.zero_common_targets == 0;
  result.starstar = result.minimum_met == regular_oids.size();
  if (minimum_target_oid < 0) {
    throw std::runtime_error("minimum target bookkeeping failed");
  }
  std::cout << "SD_MINIMUM|case=" << spec.id
            << "|label=" << spec.label
            << "|minimum_met=" << result.minimum_met
            << "|regular_orbits=" << regular_oids.size()
            << "|target_H_orbit=" << minimum_target_oid + 1
            << "|target_code=" << orbit.points[minimum_target_oid][0]
            << '\n';
  std::cout << "SD_RESULT|case=" << spec.id
            << "|label=" << spec.label
            << "|common_neighbour="
            << (result.common_neighbour ? "PASS" : "FAIL")
            << "|starstar=" << (result.starstar ? "PASS" : "FAIL")
            << "|zero_common_targets=" << result.zero_common_targets
            << "|quotient_tests=" << result.quotient_tests
            << "|reconstruction_checks=" << result.reconstruction_checks
            << "|threads=" << arguments.threads
            << "|status="
            << (result.common_neighbour
                    ? (result.starstar
                           ? "COMMON_NEIGHBOUR_AND_STARSTAR_PASS"
                           : "COMMON_NEIGHBOUR_PASS_STARSTAR_FAIL")
                    : "COUNTEREXAMPLE_CERTIFIED")
            << '\n';
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const SdArguments arguments = parse_sd_arguments(argc, argv);
    std::cout << "ENGINE|name=psl27_k4_starstar|schema=1"
              << "|mode="
              << (arguments.metadata_only
                      ? "metadata"
                      : (arguments.census_only ? "census" : "exact"))
              << "|threads=" << arguments.threads
              << "|openmp="
#ifdef _OPENMP
              << "YES"
#else
              << "NO"
#endif
              << "|exactness=integer_group_tables_H_orbits_and_quotients\n";

    std::vector<const CaseSpec*> selected;
    for (const CaseSpec& spec : kCases) {
      if (arguments.selected_case < 0 || arguments.selected_case == spec.id) {
        validate_case_metadata(spec);
        selected.push_back(&spec);
      }
    }
    if (arguments.metadata_only) {
      std::cout << "SD_RUN_RESULT|cases=" << selected.size()
                << "|metadata_pass=" << selected.size()
                << "|status=METADATA_PASS\n";
      return 0;
    }

    FiniteGroup group = psl2_7();
    verify_group(group);
    if (group.order != kSimpleOrder || group.outer_order != 2) {
      throw std::runtime_error("PSL(2,7) model metadata mismatch");
    }
    const StateSpace space = make_state_space(group.order, kDiagonalLength);
    if (space.count != kExpectedDegree) {
      throw std::runtime_error("PSL(2,7)^4 degree mismatch");
    }
    const ActionMaps maps = build_action_maps(group, space, kDiagonalLength);
    std::cout << "MODEL_VALIDATION|T=" << group.name
              << "|T_order=" << group.order
              << "|outer_order=" << group.outer_order
              << "|degree=" << space.count
              << "|action_maps=" << maps.checked_maps
              << "|status=PASS\n";

    int applicable = 0;
    int common_pass = 0;
    int starstar_pass = 0;
    int counterexamples = 0;
    for (const CaseSpec* spec : selected) {
      const ScanResult result =
          run_case(group, space, maps, *spec, arguments);
      applicable += result.applicable;
      if (!arguments.census_only && result.applicable) {
        common_pass += result.common_neighbour;
        starstar_pass += result.starstar;
        counterexamples += !result.common_neighbour;
      }
    }
    std::cout << "SD_RUN_RESULT|cases=" << selected.size()
              << "|applicable=" << applicable
              << "|common_pass=" << common_pass
              << "|starstar_pass=" << starstar_pass
              << "|counterexamples=" << counterexamples
              << "|status="
              << (arguments.census_only
                      ? "CENSUS_PASS"
                      : (counterexamples == 0
                             ? "COMMON_NEIGHBOUR_PASS"
                             : "COUNTEREXAMPLE_FOUND"))
              << '\n';
    return counterexamples == 0 ? 0 : 2;
  } catch (const std::exception& error) {
    std::cerr << "ERROR|" << error.what() << '\n';
    return 1;
  }
}

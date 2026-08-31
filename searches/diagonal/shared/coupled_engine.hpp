#ifndef BG_DIAGONAL_COUPLED_ENGINE_HPP
#define BG_DIAGONAL_COUPLED_ENGINE_HPP

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Included after diagonal_starstar_small.cpp.  This engine reuses its exact
// finite-group tables and diagonal action maps, but accepts an explicit list
// of point-stabiliser generators so outer/top graph subgroups can be modelled
// without accidentally adjoining either factor independently.

std::vector<Code> coupled_compose_maps(
    const std::vector<Code>& left,
    const std::vector<Code>& right) {
  if (left.size() != right.size()) {
    throw std::runtime_error("coupled map degree mismatch");
  }
  std::vector<Code> answer(left.size());
  for (Code point = 0; point < left.size(); ++point) {
    answer[point] = right[left[point]];
  }
  return answer;
}

bool coupled_maps_commute(
    const std::vector<Code>& first,
    const std::vector<Code>& second) {
  if (first.size() != second.size()) return false;
  for (Code point = 0; point < first.size(); ++point) {
    if (first[second[point]] != second[first[point]]) return false;
  }
  return true;
}

int coupled_map_order(
    const std::vector<Code>& map, int maximum_order) {
  std::vector<Code> power(map.size());
  std::iota(power.begin(), power.end(), Code{0});
  for (int exponent = 1; exponent <= maximum_order; ++exponent) {
    power = coupled_compose_maps(power, map);
    bool identity = true;
    for (Code point = 0; point < power.size(); ++point) {
      if (power[point] != point) {
        identity = false;
        break;
      }
    }
    if (identity) return exponent;
  }
  return -1;
}

struct CoupledOrbitData {
  std::vector<std::int32_t> id;
  std::vector<std::vector<Code>> points;
};

CoupledOrbitData coupled_stabiliser_orbits(
    const StateSpace& space,
    const std::vector<const std::vector<Code>*>& actions,
    int stabiliser_order) {
  CoupledOrbitData orbit;
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
          "coupled H-orbit size does not divide asserted H order");
    }
    orbit.points.push_back(queue);
  }
  return orbit;
}

void coupled_print_indices(const std::vector<int>& values) {
  std::cout << '[';
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << values[index] + 1;
  }
  std::cout << ']';
}

struct CoupledCaseResult {
  bool applicable = false;
  bool starstar = false;
  bool common_neighbour = false;
  int r = 0;
  int h_orbits = 0;
  std::uint64_t regular_points = 0;
  int minimum_met = 0;
  std::uint64_t minimum_common = 0;
  int starstar_defects = 0;
  int common_neighbour_defects = 0;
  std::uint64_t quotient_tests = 0;
};

CoupledCaseResult run_coupled_case(
    const FiniteGroup& group,
    const StateSpace& space,
    int k,
    const std::string& label,
    int quotient_order,
    const std::vector<const std::vector<Code>*>& actions,
    bool census_only,
    bool progress,
    bool verbose_targets) {
  const int stabiliser_order = group.order * quotient_order;
  const CoupledOrbitData orbit = coupled_stabiliser_orbits(
      space, actions, stabiliser_order);
  std::vector<int> regular_oids;
  std::vector<int> lambda_of_point(space.count, -1);
  std::uint64_t regular_points = 0;
  for (int oid = 0; oid < static_cast<int>(orbit.points.size()); ++oid) {
    if (orbit.points[oid].size() !=
        static_cast<std::size_t>(stabiliser_order)) {
      continue;
    }
    const int lambda = static_cast<int>(regular_oids.size());
    regular_oids.push_back(oid);
    regular_points += orbit.points[oid].size();
    for (Code point : orbit.points[oid]) {
      if (lambda_of_point[point] >= 0) {
        throw std::runtime_error(
            "coupled regular H-orbits overlap");
      }
      lambda_of_point[point] = lambda;
    }
  }

  CoupledCaseResult result;
  result.r = static_cast<int>(regular_oids.size());
  result.h_orbits = static_cast<int>(orbit.points.size());
  result.regular_points = regular_points;
  result.minimum_met = result.r;
  result.minimum_common =
      std::numeric_limits<std::uint64_t>::max();
  std::cout << "COUPLED_CASE|label=" << label
            << "|T=" << group.name
            << "|k=" << k
            << "|quotient_order=" << quotient_order
            << "|degree=" << space.count
            << "|H_order=" << stabiliser_order
            << "|H_orbits=" << orbit.points.size()
            << "|regular_orbits=" << result.r
            << "|regular_points=" << regular_points
            << "|density=" << std::fixed << std::setprecision(9)
            << static_cast<double>(regular_points) / space.count
            << '\n';
  if (result.r == 0) {
    result.minimum_common = 0;
    std::cout << "COUPLED_RESULT|label=" << label
              << "|status=NOT_APPLICABLE|reason=base_size_gt_2\n";
    return result;
  }
  result.applicable = true;
  if (census_only) {
    std::cout << "COUPLED_RESULT|label=" << label
              << "|status=CENSUS_ONLY|r=" << result.r << '\n';
    return result;
  }

  for (int target_oid = 0;
       target_oid < static_cast<int>(orbit.points.size());
       ++target_oid) {
    const Code target = orbit.points[target_oid][0];
    Digits inverse_digits{};
    for (int coordinate = 0;
         coordinate < space.coordinates; ++coordinate) {
      inverse_digits[coordinate] =
          group.inverse[space.digits[target][coordinate]];
    }
    std::vector<std::vector<std::uint8_t>> adjacency(
        result.r, std::vector<std::uint8_t>(result.r, 0));
    std::uint64_t common = 0;

    for (int left = 0; left < result.r; ++left) {
      for (Code lambda : orbit.points[regular_oids[left]]) {
        Digits quotient_digits{};
        for (int coordinate = 0;
             coordinate < space.coordinates; ++coordinate) {
          quotient_digits[coordinate] = group.multiply(
              space.digits[lambda][coordinate],
              inverse_digits[coordinate]);
        }
        const Code quotient = encode_state(space, quotient_digits);
        ++result.quotient_tests;
        const int right = lambda_of_point[quotient];
        if (right < 0) continue;
        ++common;
        adjacency[left][right] = 1;
        for (int coordinate = 0;
             coordinate < space.coordinates; ++coordinate) {
          if (group.multiply(
                  quotient_digits[coordinate],
                  space.digits[target][coordinate]) !=
              space.digits[lambda][coordinate]) {
            throw std::runtime_error(
                "coupled quotient reconstruction failed");
          }
        }
      }
    }

    int met = result.r;
    std::uint64_t edges = 0;
    std::vector<int> zero_rows;
    for (int left = 0; left < result.r; ++left) {
      const auto& row = adjacency[left];
      const int degree = static_cast<int>(
          std::count(row.begin(), row.end(), std::uint8_t{1}));
      met = std::min(met, degree);
      edges += degree;
      if (degree == 0) zero_rows.push_back(left);
    }
    result.minimum_met = std::min(result.minimum_met, met);
    result.minimum_common =
        std::min(result.minimum_common, common);
    result.starstar_defects += met == 0;
    if (!zero_rows.empty()) {
      std::cout << "COUPLED_STARSTAR_DEFECT|label=" << label
                << "|target_H_orbit=" << target_oid + 1
                << "|target_code=" << target
                << "|target=";
      print_digits(space, target);
      std::cout << "|zero_lambda_indices=";
      coupled_print_indices(zero_rows);
      std::cout << "|zero_lambda_representatives=[";
      for (std::size_t index = 0;
           index < zero_rows.size(); ++index) {
        if (index) std::cout << ',';
        std::cout
            << orbit.points[regular_oids[zero_rows[index]]][0];
      }
      std::cout << "]|certificate=no_point_of_Lambda_i*x^-1_in_R\n";
    }

    if (common == 0) {
      ++result.common_neighbour_defects;
      std::cout << "COUPLED_COMMON_NEIGHBOUR_DEFECT|label=" << label
                << "|target_H_orbit=" << target_oid + 1
                << "|target_code=" << target
                << "|target=";
      print_digits(space, target);
      std::cout << "|common_neighbours=0\n";
    }

    if (verbose_targets) {
      std::cout << "COUPLED_TARGET|label=" << label
                << "|H_orbit=" << target_oid + 1
                << "|target_code=" << target
                << "|common=" << common
                << "|incidence_edges=" << edges
                << "|minimum_row_degree=" << met
                << '\n';
    }
    if (progress && (target_oid + 1) % 25 == 0) {
      std::cerr << "COUPLED_PROGRESS|label=" << label
                << "|targets=" << target_oid + 1
                << '/' << orbit.points.size()
                << "|quotient_tests=" << result.quotient_tests
                << "|starstar_defects=" << result.starstar_defects
                << '\n';
    }
  }

  result.starstar = result.minimum_met > 0;
  result.common_neighbour = result.common_neighbour_defects == 0;
  std::cout << "COUPLED_STARSTAR_RESULT|label=" << label
            << "|status=" << (result.starstar ? "PASS" : "FAIL")
            << "|minimum_met=" << result.minimum_met
            << "|r=" << result.r
            << "|minimum_common=" << result.minimum_common
            << "|zero_row_targets=" << result.starstar_defects
            << "|quotient_tests=" << result.quotient_tests
            << '\n';
  std::cout << "COUPLED_COMMON_NEIGHBOUR_RESULT|label=" << label
            << "|status="
            << (result.common_neighbour ? "PASS" : "FAIL")
            << "|minimum_common=" << result.minimum_common
            << "|zero_common_targets="
            << result.common_neighbour_defects
            << '\n';
  return result;
}

#endif  // BG_DIAGONAL_COUPLED_ENGINE_HPP

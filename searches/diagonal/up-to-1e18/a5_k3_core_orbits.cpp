// Exact inner-diagonal core orbit census for A5,k=3.  This is the reduced
// colour set on which every compound quotient Q <= (Out(A5)xS3) wr S_m
// acts: a pair is a base for T^(3m).Q precisely when every coordinate lies
// in a regular inner-diagonal orbit and the resulting colour tuple is
// Q-regular.

#define main diagonal_starstar_small_embedded_main
#include "../shared/diagonal_starstar_small.cpp"
#undef main
#include "../shared/coupled_engine.hpp"

#include <map>
#include <set>

namespace {

std::vector<int> compose_small_permutations(
    const std::vector<int>& left, const std::vector<int>& right) {
  std::vector<int> answer(left.size());
  for (std::size_t i = 0; i < left.size(); ++i) {
    answer[i] = right[left[i]];
  }
  return answer;
}

}  // namespace

int main() {
  try {
    std::cout << "ENGINE|name=a5_k3_core_orbits|schema=1\n";
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
    std::cout << "CORE_CENSUS|T=A5|k=3|degree=" << space.count
              << "|inner_order=" << group.order
              << "|orbits=" << orbit.points.size()
              << "|regular_orbits=" << regular_oids.size()
              << "|regular_points=" << regular_oids.size() * group.order
              << "\n";

    auto induced = [&](const std::vector<Code>& point_map) {
      std::vector<int> answer(regular_oids.size());
      for (int colour = 0; colour < static_cast<int>(regular_oids.size());
           ++colour) {
        const Code representative = orbit.points[regular_oids[colour]][0];
        const int image_oid = orbit.id[point_map[representative]];
        const auto found = std::find(
            regular_oids.begin(), regular_oids.end(), image_oid);
        if (found == regular_oids.end()) {
          throw std::runtime_error("normaliser did not preserve core-regular orbits");
        }
        answer[colour] = static_cast<int>(found - regular_oids.begin());
      }
      return answer;
    };
    const auto top3 = induced(maps.base[group.generators.size()]);
    const auto top2 = induced(maps.symmetric);
    const auto outer = induced(maps.outer);
    auto print_generator = [](const std::string& name,
                              const std::vector<int>& image) {
      std::cout << "CORE_R_GENERATOR|name=" << name << "|images=[";
      for (std::size_t i = 0; i < image.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << image[i] + 1;
      }
      std::cout << "]\n";
    };
    print_generator("top3", top3);
    print_generator("top2", top2);
    print_generator("outer", outer);
    auto induced_target = [&](const std::vector<Code>& point_map) {
      std::vector<int> answer(orbit.points.size());
      for (int oid = 0; oid < static_cast<int>(orbit.points.size()); ++oid) {
        answer[oid] = orbit.id[point_map[orbit.points[oid][0]]];
      }
      return answer;
    };
    auto print_target_generator = [](const std::string& name,
                                     const std::vector<int>& image) {
      std::cout << "CORE_Y_GENERATOR|name=" << name << "|images=[";
      for (std::size_t i = 0; i < image.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << image[i] + 1;
      }
      std::cout << "]\n";
    };
    print_target_generator("top3", induced_target(maps.base[group.generators.size()]));
    print_target_generator("top2", induced_target(maps.symmetric));
    print_target_generator("outer", induced_target(maps.outer));
    struct S3Element {
      std::vector<int> coordinates;
      std::vector<int> colours;
    };
    const std::vector<int> coordinate_identity = {0, 1, 2};
    std::vector<int> colour_identity(regular_oids.size());
    std::iota(colour_identity.begin(), colour_identity.end(), 0);
    const std::vector<int> coordinate_top3 = {1, 2, 0};
    const std::vector<int> coordinate_top2 = {1, 0, 2};
    std::vector<S3Element> s3 = {
        {coordinate_identity, colour_identity}};
    std::set<std::vector<int>> seen_s3 = {coordinate_identity};
    for (std::size_t head = 0; head < s3.size(); ++head) {
      for (const auto& generator :
           std::vector<S3Element>{{coordinate_top3, top3},
                                  {coordinate_top2, top2}}) {
        const auto coordinates = compose_small_permutations(
            s3[head].coordinates, generator.coordinates);
        const auto colours = compose_small_permutations(
            s3[head].colours, generator.colours);
        if (seen_s3.insert(coordinates).second) {
          s3.push_back({coordinates, colours});
        }
      }
    }
    if (s3.size() != 6) throw std::runtime_error("S3 action closure mismatch");
    for (int outer_bit : {0, 1}) {
      for (const auto& element : s3) {
        const auto action = outer_bit
            ? compose_small_permutations(element.colours, outer)
            : element.colours;
        int fixed = 0;
        std::map<int, int> cycle_histogram;
        std::vector<std::uint8_t> seen(action.size(), 0);
        for (int start = 0; start < static_cast<int>(action.size()); ++start) {
          if (seen[start]) continue;
          int length = 0;
          for (int point = start; !seen[point]; point = action[point]) {
            seen[point] = 1;
            ++length;
          }
          ++cycle_histogram[length];
          fixed += length == 1;
        }
        std::cout << "CORE_R_ACTION|outer=" << outer_bit << "|s3=[";
        for (int i = 0; i < 3; ++i) {
          if (i) std::cout << ',';
          std::cout << element.coordinates[i] + 1;
        }
        std::cout << "]|fixed=" << fixed << "|cycles=";
        bool first = true;
        for (const auto& [length, count] : cycle_histogram) {
          if (!first) std::cout << ',';
          first = false;
          std::cout << length << '^' << count;
        }
        std::cout << '\n';
      }
    }

    std::set<std::vector<std::uint64_t>> unique;
    std::vector<int> target_multiplicity(regular_oids.size() * regular_oids.size(), 0);
    int minimum_row_degree = static_cast<int>(regular_oids.size());
    int maximum_row_degree = 0;
    int all_perfect = 1;
    int minimum_matching = static_cast<int>(regular_oids.size());
    for (const auto& target_orbit : orbit.points) {
      const Code target = target_orbit[0];
      Digits inverse_digits{};
      for (int coordinate = 0; coordinate < space.coordinates; ++coordinate) {
        inverse_digits[coordinate] =
            group.inverse[space.digits[target][coordinate]];
      }
      const int r = static_cast<int>(regular_oids.size());
      std::vector<std::uint64_t> rows(r, 0);
      for (int left = 0; left < r; ++left) {
        for (Code lambda : orbit.points[regular_oids[left]]) {
          Digits quotient_digits{};
          for (int coordinate = 0; coordinate < space.coordinates;
               ++coordinate) {
            quotient_digits[coordinate] = group.multiply(
                space.digits[lambda][coordinate],
                inverse_digits[coordinate]);
          }
          const int right = regular_id[encode_state(space, quotient_digits)];
          if (right >= 0) rows[left] |= std::uint64_t{1} << right;
        }
      }
      for (int left = 0; left < r; ++left) {
        const int row_degree = std::popcount(rows[left]);
        minimum_row_degree = std::min(minimum_row_degree, row_degree);
        maximum_row_degree = std::max(maximum_row_degree, row_degree);
        for (int right = 0; right < r; ++right) {
          target_multiplicity[left * r + right] += (rows[left] >> right) & 1U;
        }
      }
      unique.insert(rows);
    }
    for (const auto& rows : unique) {
      const int r = static_cast<int>(regular_oids.size());
      std::vector<std::vector<std::uint8_t>> adjacency(
          r, std::vector<std::uint8_t>(r, 0));
      for (int i = 0; i < r; ++i) {
        for (int j = 0; j < r; ++j) adjacency[i][j] = (rows[i] >> j) & 1U;
      }
      const int matching = coupled_maximum_matching(adjacency).size;
      minimum_matching = std::min(minimum_matching, matching);
      all_perfect &= matching == r;
    }
    std::cout << "CORE_MATRIX_CENSUS|target_orbits=" << orbit.points.size()
              << "|unique_matrices=" << unique.size()
              << "|all_perfect=" << all_perfect
              << "|minimum_matching=" << minimum_matching << '\n';
    const auto [min_target, max_target] = std::minmax_element(
        target_multiplicity.begin(), target_multiplicity.end());
    const std::uint64_t target_total = std::accumulate(
        target_multiplicity.begin(), target_multiplicity.end(), std::uint64_t{0});
    std::cout << "CORE_PAIR_TARGET_MULTIPLICITY|min=" << *min_target
              << "|max=" << *max_target << "|total=" << target_total
              << "|row_degree_min=" << minimum_row_degree
              << "|row_degree_max=" << maximum_row_degree << '\n';
    return all_perfect ? 0 : 2;
  } catch (const std::exception& error) {
    std::cerr << "ERROR|" << error.what() << '\n';
    return 1;
  }
}

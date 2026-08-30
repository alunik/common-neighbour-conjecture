// Exact (star-star) checker for the four simple-diagonal actions with
// T=A7, k=3, top A3/S3, without/with Out(A7)=C2.
//
// The A7 model is generated independently as even permutations of degree 7;
// the coverage engine is shared with diagonal_starstar_small.cpp.

#define main diagonal_starstar_small_embedded_main
#include "diagonal_starstar_small.cpp"
#undef main

#include <unordered_map>

namespace {

using DynamicPermutation = std::vector<std::uint8_t>;

std::string dynamic_key(const DynamicPermutation& permutation) {
  return std::string(
      reinterpret_cast<const char*>(permutation.data()),
      permutation.size());
}

DynamicPermutation dynamic_compose(
    const DynamicPermutation& left,
    const DynamicPermutation& right) {
  DynamicPermutation answer(left.size());
  for (std::size_t point = 0; point < left.size(); ++point) {
    answer[point] = right[left[point]];
  }
  return answer;
}

DynamicPermutation three_cycle_dynamic(
    int degree, int first, int second, int third) {
  DynamicPermutation answer(degree);
  std::iota(answer.begin(), answer.end(), std::uint8_t{0});
  answer[first] = static_cast<std::uint8_t>(second);
  answer[second] = static_cast<std::uint8_t>(third);
  answer[third] = static_cast<std::uint8_t>(first);
  return answer;
}

FiniteGroup alternating_group_dynamic(int degree) {
  if (degree < 5 || degree > 9) {
    throw std::runtime_error("dynamic alternating model supports 5<=n<=9");
  }
  DynamicPermutation identity(degree);
  std::iota(identity.begin(), identity.end(), std::uint8_t{0});
  std::vector<DynamicPermutation> generators;
  for (int point = 2; point < degree; ++point) {
    generators.push_back(three_cycle_dynamic(
        degree, 0, 1, point));
  }

  std::vector<DynamicPermutation> elements{identity};
  std::unordered_map<std::string, int> index;
  index.emplace(dynamic_key(identity), 0);
  for (std::size_t head = 0; head < elements.size(); ++head) {
    const DynamicPermutation current = elements[head];
    for (const auto& generator : generators) {
      DynamicPermutation image =
          dynamic_compose(current, generator);
      auto [position, inserted] = index.emplace(
          dynamic_key(image), static_cast<int>(elements.size()));
      if (inserted) elements.push_back(std::move(image));
    }
  }

  int expected = 1;
  for (int value = 2; value <= degree; ++value) expected *= value;
  expected /= 2;
  if (static_cast<int>(elements.size()) != expected) {
    throw std::runtime_error("dynamic alternating generators failed");
  }

  FiniteGroup group;
  group.name = "A" + std::to_string(degree);
  group.model =
      "even_permutations_degree_" + std::to_string(degree);
  group.outer_name = "full_C2_odd_S" + std::to_string(degree);
  group.order = expected;
  group.identity = 0;
  group.product.resize(
      static_cast<std::size_t>(expected) * expected);
  for (int left = 0; left < expected; ++left) {
    for (int right = 0; right < expected; ++right) {
      const auto position = index.find(dynamic_key(
          dynamic_compose(elements[left], elements[right])));
      if (position == index.end()) {
        throw std::runtime_error("dynamic alternating product escaped");
      }
      group.product[
          static_cast<std::size_t>(left) * expected + right] =
          static_cast<Element>(position->second);
    }
  }
  fill_inverses(group);
  for (const auto& generator : generators) {
    group.generators.push_back(static_cast<Element>(
        index.at(dynamic_key(generator))));
  }

  DynamicPermutation transposition = identity;
  std::swap(transposition[0], transposition[1]);
  group.outer_map.resize(expected);
  for (int value = 0; value < expected; ++value) {
    const DynamicPermutation image = dynamic_compose(
        dynamic_compose(transposition, elements[value]),
        transposition);
    group.outer_map[value] = static_cast<Element>(
        index.at(dynamic_key(image)));
  }
  group.outer_order = 2;
  return group;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    Arguments arguments;
    for (int index = 1; index < argc; ++index) {
      const std::string argument = argv[index];
      if (argument == "--progress") {
        arguments.progress = true;
      } else if (argument == "--help") {
        std::cout << "usage: diagonal_starstar_a7 [--progress]\n";
        return 0;
      } else {
        throw std::runtime_error("unknown argument: " + argument);
      }
    }
    std::cout << "ENGINE|name=diagonal_starstar_a7"
              << "|schema=1"
              << "|exactness=integer_permutation_tables_and_H_orbits\n";
    RunTotals totals;
    run_group(alternating_group_dynamic(7), 3, arguments, totals);
    std::cout << "RUN_RESULT|cases=" << totals.cases
              << "|pass=" << totals.pass
              << "|fail=" << totals.fail
              << "|pairs=" << totals.pairs
              << "|membership_tests=" << totals.membership_tests
              << "|status="
              << (totals.fail == 0 ? "ALL_PASS" : "COUNTEREXAMPLE_FOUND")
              << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ERROR|" << error.what() << '\n';
    return 1;
  }
}

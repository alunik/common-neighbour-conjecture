// Exact (star-star) checker for the currently modelled small simple-diagonal
// actions.
//
// Let beta be the point represented by (1,...,1), let H = G_beta, let
// Lambda_1,...,Lambda_r be the regular H-orbits, and put
//
//     R = Lambda_1 union ... union Lambda_r.
//
// In the regular-coordinate model Omega = T^(k-1), the Saxl neighbourhood of
// x is R*x (coordinatewise group multiplication).  Condition (star-star) is
//
//     (R*x) intersection Lambda_i is nonempty
//
// for every x and every i.  The condition is constant on H-orbits of x:
// N(x^h) = N(x)^h and every Lambda_i is H-invariant.  We therefore test one
// representative of every H-orbit.  For a representative x and lambda in
// Lambda_i, the hot membership test is
//
//     lambda in R*x  <=>  lambda*x^{-1} in R.
//
// This is an exact integer computation.  It does not construct the Saxl graph.
//
// Covered cases:
//   PSL(2,7), k=3, top A3/S3, without/with graph C2;
//   A6,       k=3, top A3/S3, without/with the standard odd-S6 C2;
//   A5,       k=4, top A4/S4, without/with the full outer C2;
//   PSL(2,8), k=3, top A3/S3, without/with field C3.
//
// Usage:
//   diagonal_starstar_small --all       # all sixteen actions (default)
//   diagonal_starstar_small --smoke     # PSL(2,7), S3, graph C2 only
//   diagonal_starstar_small --all --progress

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Element = std::uint16_t;
using Code = std::uint32_t;
using Digits = std::array<Element, 3>;

struct FiniteGroup {
  std::string name;
  std::string model;
  std::string outer_name;
  int order = 0;
  Element identity = 0;
  std::vector<Element> product;
  std::vector<Element> inverse;
  std::vector<Element> generators;
  std::vector<Element> outer_map;
  int outer_order = 1;

  Element multiply(Element left, Element right) const {
    return product[
        static_cast<std::size_t>(left) * static_cast<std::size_t>(order) +
        right];
  }
};

int factorial(int value) {
  int answer = 1;
  for (int i = 2; i <= value; ++i) answer *= i;
  return answer;
}

void fill_inverses(FiniteGroup& group) {
  const Element missing = static_cast<Element>(group.order);
  group.inverse.assign(group.order, missing);
  for (int value = 0; value < group.order; ++value) {
    for (int candidate = 0; candidate < group.order; ++candidate) {
      if (group.multiply(value, candidate) == group.identity &&
          group.multiply(candidate, value) == group.identity) {
        group.inverse[value] = static_cast<Element>(candidate);
        break;
      }
    }
    if (group.inverse[value] == missing) {
      throw std::runtime_error(group.name + " inverse-table failure");
    }
  }
}

void verify_group(const FiniteGroup& group) {
  if (group.order <= 0 ||
      group.product.size() !=
          static_cast<std::size_t>(group.order) * group.order ||
      group.inverse.size() != static_cast<std::size_t>(group.order)) {
    throw std::runtime_error(group.name + " malformed group table");
  }
  for (int value = 0; value < group.order; ++value) {
    const Element x = static_cast<Element>(value);
    if (group.multiply(group.identity, x) != x ||
        group.multiply(x, group.identity) != x ||
        group.multiply(x, group.inverse[x]) != group.identity ||
        group.multiply(group.inverse[x], x) != group.identity) {
      throw std::runtime_error(group.name + " identity/inverse check failed");
    }
  }

  std::vector<std::uint8_t> seen(group.order, 0);
  std::queue<Element> queue;
  seen[group.identity] = 1;
  queue.push(group.identity);
  int generated = 1;
  while (!queue.empty()) {
    const Element current = queue.front();
    queue.pop();
    for (Element generator : group.generators) {
      const Element target = group.multiply(current, generator);
      if (seen[target]) continue;
      seen[target] = 1;
      queue.push(target);
      ++generated;
    }
  }
  if (generated != group.order) {
    throw std::runtime_error(group.name + " generator check failed");
  }

  if (group.outer_order <= 1 ||
      group.outer_map.size() != static_cast<std::size_t>(group.order)) {
    throw std::runtime_error(group.name + " outer-map metadata missing");
  }
  std::fill(seen.begin(), seen.end(), 0);
  for (Element image : group.outer_map) {
    if (image >= group.order || seen[image]) {
      throw std::runtime_error(group.name + " outer map is not bijective");
    }
    seen[image] = 1;
  }
  if (group.outer_map[group.identity] != group.identity) {
    throw std::runtime_error(group.name + " outer map moves identity");
  }
  for (int left = 0; left < group.order; ++left) {
    for (int right = 0; right < group.order; ++right) {
      if (group.outer_map[group.multiply(left, right)] !=
          group.multiply(group.outer_map[left], group.outer_map[right])) {
        throw std::runtime_error(group.name + " outer map not homomorphic");
      }
    }
  }
  int exact_order = 0;
  for (int exponent = 1; exponent <= group.outer_order; ++exponent) {
    bool identity_map = true;
    for (int value = 0; value < group.order; ++value) {
      Element image = static_cast<Element>(value);
      for (int step = 0; step < exponent; ++step) {
        image = group.outer_map[image];
      }
      if (image != value) {
        identity_map = false;
        break;
      }
    }
    if (identity_map) {
      exact_order = exponent;
      break;
    }
  }
  if (exact_order != group.outer_order) {
    throw std::runtime_error(group.name + " outer-map order check failed");
  }
}

int permutation_code(
    const std::array<std::uint8_t, 6>& permutation, int degree) {
  int answer = 0;
  int place = 1;
  for (int point = 0; point < degree; ++point) {
    answer += permutation[point] * place;
    place *= degree;
  }
  return answer;
}

std::array<std::uint8_t, 6> compose_permutations(
    const std::array<std::uint8_t, 6>& left,
    const std::array<std::uint8_t, 6>& right,
    int degree) {
  // Product convention: a point is first sent by left, then by right.
  std::array<std::uint8_t, 6> answer{};
  for (int point = 0; point < degree; ++point) {
    answer[point] = right[left[point]];
  }
  return answer;
}

std::array<std::uint8_t, 6> cycle3(
    int degree, int first, int second, int third) {
  std::array<std::uint8_t, 6> answer{};
  for (int point = 0; point < degree; ++point) answer[point] = point;
  answer[first] = second;
  answer[second] = third;
  answer[third] = first;
  return answer;
}

std::array<std::uint8_t, 6> transposition(
    int degree, int first, int second) {
  std::array<std::uint8_t, 6> answer{};
  for (int point = 0; point < degree; ++point) answer[point] = point;
  std::swap(answer[first], answer[second]);
  return answer;
}

bool even_permutation(
    const std::array<std::uint8_t, 6>& permutation, int degree) {
  int inversions = 0;
  for (int first = 0; first < degree; ++first) {
    for (int second = first + 1; second < degree; ++second) {
      inversions += permutation[first] > permutation[second];
    }
  }
  return inversions % 2 == 0;
}

FiniteGroup alternating_group(int degree) {
  FiniteGroup group;
  group.name = "A" + std::to_string(degree);
  group.model = "even_permutations_degree_" + std::to_string(degree);
  group.outer_name =
      degree == 5 ? "full_C2_odd_S5" : "standard_C2_odd_S6";

  std::vector<std::array<std::uint8_t, 6>> elements;
  std::array<std::uint8_t, 6> permutation{};
  for (int point = 0; point < degree; ++point) permutation[point] = point;
  do {
    if (even_permutation(permutation, degree)) {
      elements.push_back(permutation);
    }
  } while (std::next_permutation(
      permutation.begin(), permutation.begin() + degree));

  group.order = static_cast<int>(elements.size());
  int code_count = 1;
  for (int i = 0; i < degree; ++i) code_count *= degree;
  std::vector<int> index(code_count, -1);
  for (int i = 0; i < group.order; ++i) {
    index[permutation_code(elements[i], degree)] = i;
  }
  std::array<std::uint8_t, 6> identity{};
  for (int point = 0; point < degree; ++point) identity[point] = point;
  group.identity = static_cast<Element>(
      index.at(permutation_code(identity, degree)));

  group.product.resize(
      static_cast<std::size_t>(group.order) * group.order);
  for (int left = 0; left < group.order; ++left) {
    for (int right = 0; right < group.order; ++right) {
      const auto product =
          compose_permutations(elements[left], elements[right], degree);
      group.product[
          static_cast<std::size_t>(left) * group.order + right] =
          static_cast<Element>(
              index.at(permutation_code(product, degree)));
    }
  }
  fill_inverses(group);

  for (int point = 2; point < degree; ++point) {
    group.generators.push_back(static_cast<Element>(index.at(
        permutation_code(cycle3(degree, 0, 1, point), degree))));
  }

  const auto odd = transposition(degree, 0, 1);
  group.outer_map.resize(group.order);
  for (int value = 0; value < group.order; ++value) {
    const auto conjugate = compose_permutations(
        compose_permutations(odd, elements[value], degree), odd, degree);
    group.outer_map[value] = static_cast<Element>(
        index.at(permutation_code(conjugate, degree)));
  }
  group.outer_order = 2;
  return group;
}

int gf8_multiply(int left, int right) {
  int answer = 0;
  while (right) {
    if (right & 1) answer ^= left;
    right >>= 1;
    left <<= 1;
    if (left & 8) left ^= 0b1011;  // x^3+x+1
  }
  return answer & 7;
}

int gf8_inverse(int value) {
  if (value == 0) throw std::runtime_error("inverse of zero in GF(8)");
  for (int candidate = 1; candidate < 8; ++candidate) {
    if (gf8_multiply(value, candidate) == 1) return candidate;
  }
  throw std::runtime_error("GF(8) inverse failure");
}

int matrix2_code(int a, int b, int c, int d) {
  return a | (b << 3) | (c << 6) | (d << 9);
}

std::array<int, 4> matrix2_decode(int code) {
  return {code & 7, (code >> 3) & 7, (code >> 6) & 7, (code >> 9) & 7};
}

int matrix2_multiply(int left_code, int right_code) {
  const auto left = matrix2_decode(left_code);
  const auto right = matrix2_decode(right_code);
  return matrix2_code(
      gf8_multiply(left[0], right[0]) ^ gf8_multiply(left[1], right[2]),
      gf8_multiply(left[0], right[1]) ^ gf8_multiply(left[1], right[3]),
      gf8_multiply(left[2], right[0]) ^ gf8_multiply(left[3], right[2]),
      gf8_multiply(left[2], right[1]) ^ gf8_multiply(left[3], right[3]));
}

FiniteGroup psl2_8() {
  FiniteGroup group;
  group.name = "PSL(2,8)";
  group.model = "SL(2,8)_equals_PSL(2,8)";
  group.outer_name = "field_C3";

  std::vector<int> elements;
  std::array<int, 4096> index{};
  index.fill(-1);
  for (int a = 0; a < 8; ++a) {
    for (int b = 0; b < 8; ++b) {
      for (int c = 0; c < 8; ++c) {
        for (int d = 0; d < 8; ++d) {
          const int determinant =
              gf8_multiply(a, d) ^ gf8_multiply(b, c);
          if (determinant != 1) continue;
          const int code = matrix2_code(a, b, c, d);
          index[code] = static_cast<int>(elements.size());
          elements.push_back(code);
        }
      }
    }
  }
  group.order = static_cast<int>(elements.size());
  if (group.order != 504) {
    throw std::runtime_error("SL(2,8) enumeration did not have order 504");
  }
  group.identity =
      static_cast<Element>(index[matrix2_code(1, 0, 0, 1)]);
  group.product.resize(
      static_cast<std::size_t>(group.order) * group.order);
  for (int left = 0; left < group.order; ++left) {
    for (int right = 0; right < group.order; ++right) {
      group.product[
          static_cast<std::size_t>(left) * group.order + right] =
          static_cast<Element>(
              index[matrix2_multiply(elements[left], elements[right])]);
    }
  }
  fill_inverses(group);

  const int primitive = 2;  // x modulo x^3+x+1
  const int primitive_inverse = gf8_inverse(primitive);
  group.generators = {
      static_cast<Element>(index[matrix2_code(1, 1, 0, 1)]),
      static_cast<Element>(index[matrix2_code(1, 0, 1, 1)]),
      static_cast<Element>(index[
          matrix2_code(primitive, 0, 0, primitive_inverse)]),
  };
  group.outer_map.resize(group.order);
  for (int value = 0; value < group.order; ++value) {
    auto matrix = matrix2_decode(elements[value]);
    for (int& entry : matrix) entry = gf8_multiply(entry, entry);
    group.outer_map[value] = static_cast<Element>(
        index[matrix2_code(matrix[0], matrix[1], matrix[2], matrix[3])]);
  }
  group.outer_order = 3;
  return group;
}

using Matrix3 = std::uint16_t;

int matrix3_entry(Matrix3 matrix, int row, int column) {
  return (matrix >> (3 * row + column)) & 1;
}

Matrix3 matrix3_multiply(Matrix3 left, Matrix3 right) {
  Matrix3 answer = 0;
  for (int row = 0; row < 3; ++row) {
    for (int column = 0; column < 3; ++column) {
      int value = 0;
      for (int inner = 0; inner < 3; ++inner) {
        value ^= matrix3_entry(left, row, inner) &
                 matrix3_entry(right, inner, column);
      }
      answer |= static_cast<Matrix3>(value << (3 * row + column));
    }
  }
  return answer;
}

Matrix3 matrix3_transpose(Matrix3 matrix) {
  Matrix3 answer = 0;
  for (int row = 0; row < 3; ++row) {
    for (int column = 0; column < 3; ++column) {
      answer |= static_cast<Matrix3>(
          matrix3_entry(matrix, row, column) <<
          (3 * column + row));
    }
  }
  return answer;
}

bool matrix3_invertible(Matrix3 matrix) {
  std::array<int, 3> rows{};
  for (int row = 0; row < 3; ++row) {
    rows[row] = (matrix >> (3 * row)) & 7;
  }
  int rank = 0;
  for (int column = 0; column < 3; ++column) {
    int pivot = rank;
    while (pivot < 3 && ((rows[pivot] >> column) & 1) == 0) ++pivot;
    if (pivot == 3) continue;
    std::swap(rows[rank], rows[pivot]);
    for (int row = 0; row < 3; ++row) {
      if (row != rank && ((rows[row] >> column) & 1)) {
        rows[row] ^= rows[rank];
      }
    }
    ++rank;
  }
  return rank == 3;
}

FiniteGroup psl2_7() {
  FiniteGroup group;
  group.name = "PSL(2,7)";
  group.model = "GL(3,2)_equals_PSL(2,7)";
  group.outer_name = "graph_C2";

  std::vector<Matrix3> elements;
  std::array<int, 512> index{};
  index.fill(-1);
  for (Matrix3 matrix = 0; matrix < 512; ++matrix) {
    if (!matrix3_invertible(matrix)) continue;
    index[matrix] = static_cast<int>(elements.size());
    elements.push_back(matrix);
  }
  group.order = static_cast<int>(elements.size());
  if (group.order != 168) {
    throw std::runtime_error("GL(3,2) enumeration did not have order 168");
  }
  const Matrix3 identity = (1U << 0) | (1U << 4) | (1U << 8);
  group.identity = static_cast<Element>(index[identity]);
  group.product.resize(
      static_cast<std::size_t>(group.order) * group.order);
  for (int left = 0; left < group.order; ++left) {
    for (int right = 0; right < group.order; ++right) {
      group.product[
          static_cast<std::size_t>(left) * group.order + right] =
          static_cast<Element>(
              index[matrix3_multiply(elements[left], elements[right])]);
    }
  }
  fill_inverses(group);

  const Matrix3 cycle =
      (1U << 1) | (1U << 5) | (1U << 6);
  const Matrix3 transvection = identity | (1U << 1);
  group.generators = {
      static_cast<Element>(index[cycle]),
      static_cast<Element>(index[transvection]),
  };
  group.outer_map.resize(group.order);
  for (int value = 0; value < group.order; ++value) {
    const Matrix3 image =
        matrix3_transpose(elements[group.inverse[value]]);
    group.outer_map[value] = static_cast<Element>(index[image]);
    if (index[image] < 0) {
      throw std::runtime_error("graph automorphism left GL(3,2)");
    }
  }
  group.outer_order = 2;
  return group;
}

struct StateSpace {
  int coordinates = 0;
  Code count = 0;
  std::array<Code, 4> powers{};
  std::vector<Digits> digits;
};

StateSpace make_state_space(int group_order, int k) {
  StateSpace space;
  space.coordinates = k - 1;
  space.powers.fill(1);
  for (int coordinate = 1; coordinate <= space.coordinates; ++coordinate) {
    space.powers[coordinate] =
        space.powers[coordinate - 1] * group_order;
  }
  space.count = space.powers[space.coordinates];
  space.digits.resize(space.count);
  for (Code code = 0; code < space.count; ++code) {
    Code value = code;
    for (int coordinate = 0; coordinate < space.coordinates; ++coordinate) {
      space.digits[code][coordinate] =
          static_cast<Element>(value % group_order);
      value /= group_order;
    }
  }
  return space;
}

Code encode_state(const StateSpace& space, const Digits& digits) {
  Code answer = 0;
  for (int coordinate = 0; coordinate < space.coordinates; ++coordinate) {
    answer += static_cast<Code>(digits[coordinate]) *
              space.powers[coordinate];
  }
  return answer;
}

std::vector<std::vector<int>> alternating_top_generators(int k) {
  std::vector<std::vector<int>> generators;
  for (int point = 2; point < k; ++point) {
    std::vector<int> cycle(k);
    std::iota(cycle.begin(), cycle.end(), 0);
    cycle[0] = 1;
    cycle[1] = point;
    cycle[point] = 0;
    generators.push_back(std::move(cycle));
  }
  return generators;
}

std::vector<int> symmetric_top_generator(int k) {
  std::vector<int> swap(k);
  std::iota(swap.begin(), swap.end(), 0);
  std::swap(swap[0], swap[1]);
  return swap;
}

struct ActionMaps {
  std::vector<std::vector<Code>> base;
  std::vector<Code> symmetric;
  std::vector<Code> outer;
  int checked_maps = 0;
};

void verify_action_map(
    const std::vector<Code>& image, Code origin, const std::string& label) {
  if (image.empty() || image[origin] != origin) {
    throw std::runtime_error(label + " action does not fix origin");
  }
  std::vector<std::uint8_t> seen(image.size(), 0);
  for (Code target : image) {
    if (target >= image.size() || seen[target]) {
      throw std::runtime_error(label + " action is not a permutation");
    }
    seen[target] = 1;
  }
}

ActionMaps build_action_maps(
    const FiniteGroup& group, const StateSpace& space, int k) {
  ActionMaps maps;
  Digits origin_digits{};
  origin_digits.fill(group.identity);
  const Code origin = encode_state(space, origin_digits);

  for (Element conjugator : group.generators) {
    std::vector<Code> image(space.count);
    const Element inverse = group.inverse[conjugator];
    for (Code point = 0; point < space.count; ++point) {
      Digits output{};
      for (int coordinate = 0;
           coordinate < space.coordinates; ++coordinate) {
        output[coordinate] = group.multiply(
            group.multiply(
                inverse, space.digits[point][coordinate]),
            conjugator);
      }
      image[point] = encode_state(space, output);
    }
    verify_action_map(image, origin, group.name + "_inner");
    ++maps.checked_maps;
    maps.base.push_back(std::move(image));
  }

  auto build_top_map = [&](const std::vector<int>& permutation) {
    std::vector<Code> image(space.count);
    for (Code point = 0; point < space.count; ++point) {
      std::array<Element, 4> tuple{};
      tuple[0] = group.identity;
      for (int coordinate = 1; coordinate < k; ++coordinate) {
        tuple[coordinate] =
            space.digits[point][coordinate - 1];
      }
      std::array<Element, 4> permuted{};
      for (int coordinate = 0; coordinate < k; ++coordinate) {
        permuted[coordinate] = tuple[permutation[coordinate]];
      }
      const Element normaliser = group.inverse[permuted[0]];
      Digits output{};
      for (int coordinate = 1; coordinate < k; ++coordinate) {
        output[coordinate - 1] =
            group.multiply(normaliser, permuted[coordinate]);
      }
      image[point] = encode_state(space, output);
    }
    return image;
  };

  for (const auto& permutation : alternating_top_generators(k)) {
    std::vector<Code> image = build_top_map(permutation);
    verify_action_map(image, origin, group.name + "_alternating_top");
    ++maps.checked_maps;
    maps.base.push_back(std::move(image));
  }

  maps.symmetric = build_top_map(symmetric_top_generator(k));
  verify_action_map(maps.symmetric, origin, group.name + "_symmetric_top");
  ++maps.checked_maps;

  maps.outer.resize(space.count);
  for (Code point = 0; point < space.count; ++point) {
    Digits output{};
    for (int coordinate = 0;
         coordinate < space.coordinates; ++coordinate) {
      output[coordinate] =
          group.outer_map[space.digits[point][coordinate]];
    }
    maps.outer[point] = encode_state(space, output);
  }
  verify_action_map(maps.outer, origin, group.name + "_outer");
  ++maps.checked_maps;
  return maps;
}

struct OrbitData {
  std::vector<std::int32_t> id;
  std::vector<std::vector<Code>> points;
};

OrbitData stabiliser_orbits(
    const StateSpace& space,
    const ActionMaps& maps,
    bool symmetric_top,
    bool include_outer,
    int stabiliser_order) {
  std::vector<const std::vector<Code>*> actions;
  actions.reserve(maps.base.size() + 2);
  for (const auto& action : maps.base) actions.push_back(&action);
  if (symmetric_top) actions.push_back(&maps.symmetric);
  if (include_outer) actions.push_back(&maps.outer);

  OrbitData orbit;
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

void print_digits(const StateSpace& space, Code code) {
  std::cout << '[';
  for (int coordinate = 0; coordinate < space.coordinates; ++coordinate) {
    if (coordinate) std::cout << ',';
    std::cout << space.digits[code][coordinate];
  }
  std::cout << ']';
}

void print_index_list(const std::vector<std::uint8_t>& met, bool wanted) {
  std::cout << '[';
  bool first = true;
  for (std::size_t index = 0; index < met.size(); ++index) {
    if (static_cast<bool>(met[index]) != wanted) continue;
    if (!first) std::cout << ',';
    std::cout << index + 1;
    first = false;
  }
  std::cout << ']';
}

struct CaseResult {
  bool pass = false;
  std::uint64_t pairs = 0;
  std::uint64_t membership_tests = 0;
  std::uint64_t reconstruction_checks = 0;
  std::uint64_t missing_pairs = 0;
  std::size_t minimum_met = 0;
};

CaseResult run_case(
    const FiniteGroup& group,
    const StateSpace& space,
    const ActionMaps& maps,
    int k,
    bool symmetric_top,
    bool include_outer,
    bool progress) {
  const auto start = std::chrono::steady_clock::now();
  const std::string top =
      std::string(symmetric_top ? "S" : "A") + std::to_string(k);
  const int top_order =
      factorial(k) / (symmetric_top ? 1 : 2);
  const int stabiliser_order =
      group.order * top_order *
      (include_outer ? group.outer_order : 1);
  const OrbitData orbit = stabiliser_orbits(
      space, maps, symmetric_top, include_outer, stabiliser_order);

  std::vector<int> regular_oids;
  std::vector<std::uint8_t> regular(space.count, 0);
  std::uint64_t regular_points = 0;
  for (int oid = 0; oid < static_cast<int>(orbit.points.size()); ++oid) {
    if (orbit.points[oid].size() !=
        static_cast<std::size_t>(stabiliser_order)) {
      continue;
    }
    regular_oids.push_back(oid);
    regular_points += orbit.points[oid].size();
    for (Code point : orbit.points[oid]) regular[point] = 1;
  }

  std::uint64_t inverse_checks = 0;
  for (int oid : regular_oids) {
    for (Code point : orbit.points[oid]) {
      Digits inverse_digits{};
      for (int coordinate = 0;
           coordinate < space.coordinates; ++coordinate) {
        inverse_digits[coordinate] =
            group.inverse[space.digits[point][coordinate]];
      }
      if (!regular[encode_state(space, inverse_digits)]) {
        throw std::runtime_error(
            group.name + " connection set is not inverse-stable");
      }
      ++inverse_checks;
    }
  }

  std::cout << "CASE|T=" << group.name
            << "|k=" << k
            << "|top=" << top
            << "|outer="
            << (include_outer ? group.outer_name : "none")
            << "|outer_order="
            << (include_outer ? group.outer_order : 1)
            << "|degree=" << space.count
            << "|H_order=" << stabiliser_order
            << "|H_orbits=" << orbit.points.size()
            << "|regular_orbits=" << regular_oids.size()
            << "|regular_points=" << regular_points
            << "|density=" << std::fixed << std::setprecision(9)
            << static_cast<double>(regular_points) / space.count
            << '\n';
  std::cout << "CONVENTION_CHECK|T=" << group.name
            << "|top=" << top
            << "|outer="
            << (include_outer ? group.outer_name : "none")
            << "|coset_representative=(1,x1,...,x"
            << (k - 1) << ")"
            << "|point_product=coordinatewise_left_times_right"
            << "|neighbourhood=R*x"
            << "|membership=lambda*x^-1_in_R"
            << "|R_inverse_checks=" << inverse_checks
            << "|status=PASS\n";

  CaseResult result;
  result.pairs =
      static_cast<std::uint64_t>(orbit.points.size()) *
      regular_oids.size();
  if (regular_oids.empty()) {
    std::cout << "RESULT|starstar=NOT_APPLICABLE|reason=base_size_gt_2\n";
    return result;
  }

  for (std::size_t lambda_index = 0;
       lambda_index < regular_oids.size(); ++lambda_index) {
    const int oid = regular_oids[lambda_index];
    const Code representative = orbit.points[oid][0];
    std::cout << "LAMBDA|index=" << lambda_index + 1
              << "|H_orbit=" << oid + 1
              << "|representative_code=" << representative
              << "|representative=";
    print_digits(space, representative);
    std::cout << "|size=" << orbit.points[oid].size() << '\n';
  }

  std::size_t minimum_met = regular_oids.size() + 1;
  int minimum_target_oid = -1;
  std::vector<std::uint8_t> minimum_met_bitmap;
  for (int target_oid = 0;
       target_oid < static_cast<int>(orbit.points.size()); ++target_oid) {
    const Code target = orbit.points[target_oid][0];
    const Digits& target_digits = space.digits[target];
    Digits target_inverse{};
    for (int coordinate = 0;
         coordinate < space.coordinates; ++coordinate) {
      target_inverse[coordinate] =
          group.inverse[target_digits[coordinate]];
    }

    std::vector<std::uint8_t> met(regular_oids.size(), 0);
    std::uint64_t target_tests = 0;
    for (std::size_t lambda_index = 0;
         lambda_index < regular_oids.size(); ++lambda_index) {
      const auto& lambda = orbit.points[regular_oids[lambda_index]];
      for (Code lambda_point : lambda) {
        ++target_tests;
        ++result.membership_tests;
        Code right_quotient = 0;
        Digits quotient_digits{};
        for (int coordinate = 0;
             coordinate < space.coordinates; ++coordinate) {
          quotient_digits[coordinate] = group.multiply(
              space.digits[lambda_point][coordinate],
              target_inverse[coordinate]);
          right_quotient +=
              static_cast<Code>(quotient_digits[coordinate]) *
              space.powers[coordinate];
        }
        if (!regular[right_quotient]) continue;

        // This exact reconstruction check guards the multiplication direction:
        // r = lambda*x^-1 must satisfy r*x = lambda coordinate by coordinate.
        for (int coordinate = 0;
             coordinate < space.coordinates; ++coordinate) {
          if (group.multiply(
                  quotient_digits[coordinate],
                  target_digits[coordinate]) !=
              space.digits[lambda_point][coordinate]) {
            throw std::runtime_error(
                group.name + " right-quotient reconstruction failed");
          }
        }
        ++result.reconstruction_checks;
        met[lambda_index] = 1;
        break;
      }
    }

    const std::size_t met_count =
        static_cast<std::size_t>(
            std::count(met.begin(), met.end(), std::uint8_t{1}));
    const std::size_t missing_count = met.size() - met_count;
    result.missing_pairs += missing_count;
    if (met_count < minimum_met) {
      minimum_met = met_count;
      minimum_target_oid = target_oid;
      minimum_met_bitmap = met;
    }

    std::cout << "TARGET|H_orbit=" << target_oid + 1
              << "|representative_code=" << target
              << "|representative=";
    print_digits(space, target);
    std::cout << "|orbit_size=" << orbit.points[target_oid].size()
              << "|met_count=" << met_count
              << "|regular_orbits=" << met.size()
              << "|met=";
    if (missing_count == 0) {
      std::cout << "ALL";
    } else {
      print_index_list(met, true);
    }
    std::cout << "|missing=";
    print_index_list(met, false);
    std::cout << "|membership_tests=" << target_tests << '\n';

    if (progress && (target_oid + 1) % 25 == 0) {
      const double seconds = std::chrono::duration<double>(
          std::chrono::steady_clock::now() - start).count();
      std::cerr << "PROGRESS|T=" << group.name
                << "|top=" << top
                << "|outer="
                << (include_outer ? group.outer_name : "none")
                << "|targets=" << target_oid + 1
                << "/" << orbit.points.size()
                << "|wall_seconds=" << std::fixed
                << std::setprecision(3) << seconds << '\n';
    }
  }

  if (minimum_target_oid < 0 || minimum_met > regular_oids.size()) {
    throw std::runtime_error(group.name + " minimum-coverage bookkeeping failed");
  }
  result.minimum_met = minimum_met;
  const Code minimum_target =
      orbit.points[minimum_target_oid][0];
  std::cout << "MINIMUM|m=" << minimum_met
            << "|r=" << regular_oids.size()
            << "|target_H_orbit=" << minimum_target_oid + 1
            << "|target_representative_code=" << minimum_target
            << "|target_representative=";
  print_digits(space, minimum_target);
  std::cout << "|met=";
  if (minimum_met == regular_oids.size()) {
    std::cout << "ALL";
  } else {
    print_index_list(minimum_met_bitmap, true);
  }
  std::cout << "|missing=";
  print_index_list(minimum_met_bitmap, false);
  std::cout << '\n';

  const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  if (result.missing_pairs == 0) {
    result.pass = true;
    std::cout << "RESULT|starstar=PASS"
              << "|target_H_orbits=" << orbit.points.size()
              << "|regular_orbits=" << regular_oids.size()
              << "|minimum_met_regular_orbits=" << minimum_met
              << "|pairs=" << result.pairs
              << "|membership_tests=" << result.membership_tests
              << "|reconstruction_checks="
              << result.reconstruction_checks
              << "|wall_seconds=" << std::fixed
              << std::setprecision(6) << seconds << '\n';
  } else {
    const int missing_lambda = static_cast<int>(
        std::find(
            minimum_met_bitmap.begin(),
            minimum_met_bitmap.end(),
            std::uint8_t{0}) -
        minimum_met_bitmap.begin());
    const int lambda_oid = regular_oids[missing_lambda];
    std::cout << "FAILURE|target_H_orbit="
              << minimum_target_oid + 1
              << "|target_representative_code=" << minimum_target
              << "|target_representative=";
    print_digits(space, minimum_target);
    std::cout << "|minimum_met_regular_orbits=" << minimum_met
              << "|missing_lambda_indices=";
    print_index_list(minimum_met_bitmap, false);
    std::cout << "|first_missing_lambda_index=" << missing_lambda + 1
              << "|lambda_H_orbit=" << lambda_oid + 1
              << "|lambda_representative_code="
              << orbit.points[lambda_oid][0]
              << "|lambda_representative=";
    print_digits(space, orbit.points[lambda_oid][0]);
    std::cout << '\n';
    std::cout << "RESULT|starstar=FAIL"
              << "|target_H_orbits=" << orbit.points.size()
              << "|regular_orbits=" << regular_oids.size()
              << "|minimum_met_regular_orbits=" << minimum_met
              << "|pairs=" << result.pairs
              << "|missing_pairs=" << result.missing_pairs
              << "|membership_tests=" << result.membership_tests
              << "|reconstruction_checks="
              << result.reconstruction_checks
              << "|wall_seconds=" << std::fixed
              << std::setprecision(6) << seconds << '\n';
  }
  return result;
}

struct Arguments {
  bool smoke = false;
  bool progress = false;
};

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--all") {
      arguments.smoke = false;
    } else if (argument == "--smoke") {
      arguments.smoke = true;
    } else if (argument == "--progress") {
      arguments.progress = true;
    } else if (argument == "--help") {
      std::cout
          << "usage: diagonal_starstar_small "
          << "[--all|--smoke] [--progress]\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown argument: " + argument);
    }
  }
  return arguments;
}

struct RunTotals {
  int cases = 0;
  int pass = 0;
  int fail = 0;
  std::uint64_t pairs = 0;
  std::uint64_t membership_tests = 0;
};

void run_group(
    FiniteGroup group,
    int k,
    const Arguments& arguments,
    RunTotals& totals) {
  verify_group(group);
  std::cout << "GROUP|name=" << group.name
            << "|model=" << group.model
            << "|order=" << group.order
            << "|identity_index=" << group.identity
            << "|generators=" << group.generators.size()
            << "|outer=" << group.outer_name
            << "|outer_order=" << group.outer_order
            << "|table_checks=PASS"
            << "|generator_checks=PASS"
            << "|outer_automorphism_checks=PASS\n";

  const StateSpace space = make_state_space(group.order, k);
  const ActionMaps maps = build_action_maps(group, space, k);
  std::cout << "ACTION_CHECK|T=" << group.name
            << "|k=" << k
            << "|degree=" << space.count
            << "|maps_checked=" << maps.checked_maps
            << "|bijections=PASS"
            << "|origin_fixed=PASS\n";

  auto accumulate = [&](bool symmetric_top, bool include_outer) {
    const CaseResult result = run_case(
        group, space, maps, k, symmetric_top, include_outer,
        arguments.progress);
    ++totals.cases;
    totals.pass += result.pass;
    totals.fail += !result.pass;
    totals.pairs += result.pairs;
    totals.membership_tests += result.membership_tests;
  };

  if (arguments.smoke) {
    accumulate(true, true);
    return;
  }
  for (bool symmetric_top : {false, true}) {
    for (bool include_outer : {false, true}) {
      accumulate(symmetric_top, include_outer);
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    std::cout << "ENGINE|name=diagonal_starstar_small"
              << "|schema=1"
              << "|mode=" << (arguments.smoke ? "smoke" : "all")
              << "|exactness=integer_group_tables_and_H_orbits\n";
    std::cout << "MODEL|Omega=Diag(T)_backslash_T^k"
              << "|coordinates=T^(k-1)"
              << "|origin=(1,...,1)"
              << "|R=union_of_regular_H_orbits"
              << "|neighbours_of_x=R*x"
              << "|target_reduction=one_representative_per_H_orbit"
              << "|reason=N(x^h)=N(x)^h_and_Lambda_i^h=Lambda_i\n";

    RunTotals totals;
    run_group(psl2_7(), 3, arguments, totals);
    if (!arguments.smoke) {
      run_group(alternating_group(6), 3, arguments, totals);
      run_group(alternating_group(5), 4, arguments, totals);
      run_group(psl2_8(), 3, arguments, totals);
    }
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

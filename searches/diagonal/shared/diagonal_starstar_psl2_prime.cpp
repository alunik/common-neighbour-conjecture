// Exact (star-star) checker for simple-diagonal groups with
// T = PSL(2,p), k = 3 and top A3 or S3.
//
// This reuses the audited orbit/coverage engine in
// diagonal_starstar_small.cpp and supplies an independent permutation model
// for PSL(2,p) on the projective line.  The outer involution is conjugation by
// an involution in PGL(2,p) \ PSL(2,p).
//
// Usage:
//   diagonal_starstar_psl2_prime 7       # local cross-check
//   diagonal_starstar_psl2_prime 11 13   # next exact cases
//   diagonal_starstar_psl2_prime 17 --progress

#define main diagonal_starstar_small_embedded_main
#include "diagonal_starstar_small.cpp"
#undef main

#include <unordered_map>

namespace {

using ProjectivePermutation = std::vector<std::uint8_t>;

int mod_prime(int value, int prime) {
  value %= prime;
  return value < 0 ? value + prime : value;
}

int inverse_mod_prime(int value, int prime) {
  value = mod_prime(value, prime);
  if (value == 0) {
    throw std::runtime_error("attempted to invert zero modulo prime");
  }
  for (int candidate = 1; candidate < prime; ++candidate) {
    if (value * candidate % prime == 1) return candidate;
  }
  throw std::runtime_error("modulus is not prime");
}

bool is_prime(int value) {
  if (value < 2) return false;
  for (int divisor = 2; divisor * divisor <= value; ++divisor) {
    if (value % divisor == 0) return false;
  }
  return true;
}

bool is_square_mod_prime(int value, int prime) {
  value = mod_prime(value, prime);
  for (int candidate = 0; candidate < prime; ++candidate) {
    if (candidate * candidate % prime == value) return true;
  }
  return false;
}

std::string permutation_key(const ProjectivePermutation& permutation) {
  return std::string(
      reinterpret_cast<const char*>(permutation.data()),
      permutation.size());
}

ProjectivePermutation compose_projective(
    const ProjectivePermutation& left,
    const ProjectivePermutation& right) {
  if (left.size() != right.size()) {
    throw std::runtime_error("projective permutation degree mismatch");
  }
  ProjectivePermutation answer(left.size());
  // Product convention: a point is first sent by left, then by right.
  for (std::size_t point = 0; point < left.size(); ++point) {
    answer[point] = right[left[point]];
  }
  return answer;
}

FiniteGroup psl2_prime(int prime) {
  if (!is_prime(prime) || prime < 5 || prime % 2 == 0 ||
      prime + 1 > 255) {
    throw std::runtime_error(
        "PSL(2,p) constructor requires an odd prime 5 <= p <= 251");
  }

  const int infinity = prime;
  const int degree = prime + 1;
  ProjectivePermutation identity(degree);
  std::iota(identity.begin(), identity.end(), std::uint8_t{0});

  // Images for u:x -> x+1 and s:x -> -1/x.
  ProjectivePermutation unipotent = identity;
  for (int point = 0; point < prime; ++point) {
    unipotent[point] =
        static_cast<std::uint8_t>((point + 1) % prime);
  }
  unipotent[infinity] = static_cast<std::uint8_t>(infinity);

  ProjectivePermutation inversion = identity;
  inversion[0] = static_cast<std::uint8_t>(infinity);
  inversion[infinity] = 0;
  for (int point = 1; point < prime; ++point) {
    inversion[point] = static_cast<std::uint8_t>(
        mod_prime(-inverse_mod_prime(point, prime), prime));
  }

  std::vector<ProjectivePermutation> elements;
  std::unordered_map<std::string, int> index;
  elements.push_back(identity);
  index.emplace(permutation_key(identity), 0);
  std::vector<int> queue{0};
  const std::array<ProjectivePermutation, 2> generators = {
      unipotent, inversion};
  for (std::size_t head = 0; head < queue.size(); ++head) {
    const ProjectivePermutation current = elements[queue[head]];
    for (const auto& generator : generators) {
      ProjectivePermutation image =
          compose_projective(current, generator);
      const std::string key = permutation_key(image);
      auto [position, inserted] =
          index.emplace(key, static_cast<int>(elements.size()));
      if (!inserted) continue;
      elements.push_back(std::move(image));
      queue.push_back(position->second);
    }
  }

  const int expected_order = prime * (prime * prime - 1) / 2;
  if (static_cast<int>(elements.size()) != expected_order) {
    throw std::runtime_error(
        "PSL(2,p) projective generators have unexpected order");
  }

  FiniteGroup group;
  group.name = "PSL(2," + std::to_string(prime) + ")";
  group.model =
      "projective_line_degree_" + std::to_string(degree);
  group.outer_name = "PGL_outer_C2";
  group.order = expected_order;
  group.identity = 0;
  group.product.resize(
      static_cast<std::size_t>(group.order) * group.order);
  for (int left = 0; left < group.order; ++left) {
    for (int right = 0; right < group.order; ++right) {
      const std::string key = permutation_key(
          compose_projective(elements[left], elements[right]));
      const auto position = index.find(key);
      if (position == index.end()) {
        throw std::runtime_error("PSL(2,p) multiplication left group");
      }
      group.product[
          static_cast<std::size_t>(left) * group.order + right] =
          static_cast<Element>(position->second);
    }
  }
  fill_inverses(group);
  group.generators = {
      static_cast<Element>(index.at(permutation_key(unipotent))),
      static_cast<Element>(index.at(permutation_key(inversion))),
  };

  // Choose a so that h:x -> a/x is an involution whose determinant -a is
  // nonsquare.  Thus h lies in PGL(2,p) \ PSL(2,p).
  int scalar = -1;
  for (int candidate = 1; candidate < prime; ++candidate) {
    if (!is_square_mod_prime(-candidate, prime)) {
      scalar = candidate;
      break;
    }
  }
  if (scalar < 0) {
    throw std::runtime_error("failed to construct PGL outer involution");
  }
  ProjectivePermutation outer_conjugator = identity;
  outer_conjugator[0] = static_cast<std::uint8_t>(infinity);
  outer_conjugator[infinity] = 0;
  for (int point = 1; point < prime; ++point) {
    outer_conjugator[point] = static_cast<std::uint8_t>(
        scalar * inverse_mod_prime(point, prime) % prime);
  }
  if (compose_projective(outer_conjugator, outer_conjugator) !=
      identity) {
    throw std::runtime_error("PGL outer conjugator is not involutory");
  }

  group.outer_map.resize(group.order);
  for (int value = 0; value < group.order; ++value) {
    const ProjectivePermutation conjugate = compose_projective(
        compose_projective(outer_conjugator, elements[value]),
        outer_conjugator);
    const auto position = index.find(permutation_key(conjugate));
    if (position == index.end()) {
      throw std::runtime_error("PGL conjugation did not normalize PSL");
    }
    group.outer_map[value] =
        static_cast<Element>(position->second);
  }
  group.outer_order = 2;
  return group;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    Arguments arguments;
    std::vector<int> primes;
    for (int index = 1; index < argc; ++index) {
      const std::string argument = argv[index];
      if (argument == "--progress") {
        arguments.progress = true;
      } else if (argument == "--help") {
        std::cout
            << "usage: diagonal_starstar_psl2_prime "
            << "p [p ...] [--progress]\n";
        return 0;
      } else {
        primes.push_back(std::stoi(argument));
      }
    }
    if (primes.empty()) primes = {11, 13};

    std::cout << "ENGINE|name=diagonal_starstar_psl2_prime"
              << "|schema=1"
              << "|exactness=integer_permutation_tables_and_H_orbits\n";
    RunTotals totals;
    for (int prime : primes) {
      run_group(psl2_prime(prime), 3, arguments, totals);
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

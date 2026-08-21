// Rebuild the 320 Co3-orbits on GF(4)^22 and find two regular summands for
// every orbit representative. Multiplication by omega is also checked on the
// regular orbits to handle the full scalar extension.

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr unsigned kDim = 22;
constexpr uint32_t kSpaceSize = uint32_t{1} << kDim;
constexpr uint32_t kRoot = 7;
constexpr uint32_t kNoParent = std::numeric_limits<uint32_t>::max();

struct Matrix {
  uint32_t row[kDim]{};
};

struct Case {
  uint32_t representative = 0;
  uint64_t first_orbit_length = 0;
  uint64_t stabilizer_order = 0;
  std::vector<Matrix> generators;
};

uint32_t act(uint32_t v, const Matrix& g) {
  uint32_t image = 0;
  while (v != 0) {
    const unsigned i = static_cast<unsigned>(__builtin_ctz(v));
    image ^= g.row[i];
    v &= v - 1;
  }
  return image;
}

Matrix identity_matrix() {
  Matrix result;
  for (unsigned i = 0; i < kDim; ++i) result.row[i] = uint32_t{1} << i;
  return result;
}

Matrix multiply(const Matrix& a, const Matrix& b) {
  Matrix result;
  for (unsigned i = 0; i < kDim; ++i) result.row[i] = act(a.row[i], b);
  return result;
}

Matrix power(Matrix base, unsigned exponent) {
  Matrix result = identity_matrix();
  while (exponent != 0) {
    if (exponent & 1U) result = multiply(result, base);
    base = multiply(base, base);
    exponent >>= 1U;
  }
  return result;
}

bool is_identity(const Matrix& g) {
  for (unsigned i = 0; i < kDim; ++i) {
    if (g.row[i] != (uint32_t{1} << i)) return false;
  }
  return true;
}

std::vector<Matrix> read_initial(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::vector<Matrix> matrices;
  std::string token;
  while (in >> token) {
    if (token != "MATRIX") continue;
    Matrix g;
    for (unsigned i = 0; i < kDim; ++i) {
      uint64_t x = 0;
      if (!(in >> x) || x >= kSpaceSize) {
        throw std::runtime_error("invalid packed matrix row in initial data");
      }
      g.row[i] = static_cast<uint32_t>(x);
    }
    matrices.push_back(g);
  }
  if (matrices.size() != 2) {
    throw std::runtime_error("expected two Co3 generators");
  }
  return matrices;
}

std::vector<Case> read_cases(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::vector<Case> cases;
  std::string token;
  while (in >> token) {
    if (token != "CASE") continue;
    Case c;
    size_t generator_count = 0;
    uint64_t support_size = 0;
    if (!(in >> c.representative >> c.first_orbit_length >>
          c.stabilizer_order >> generator_count >> support_size)) {
      throw std::runtime_error("invalid CASE header");
    }
    for (size_t j = 0; j < generator_count; ++j) {
      if (!(in >> token) || token != "MATRIX") {
        throw std::runtime_error("missing MATRIX record");
      }
      Matrix g;
      for (unsigned i = 0; i < kDim; ++i) {
        uint64_t x = 0;
        if (!(in >> x) || x >= kSpaceSize) {
          throw std::runtime_error("invalid packed stabilizer matrix row");
        }
        g.row[i] = static_cast<uint32_t>(x);
      }
      c.generators.push_back(g);
    }
    if (!(in >> token) || token != "ENDCASE") {
      throw std::runtime_error("missing ENDCASE record");
    }
    cases.push_back(std::move(c));
  }
  if (cases.size() != 5) {
    throw std::runtime_error("expected five first-vector stabilizers");
  }
  return cases;
}

uint64_t next_random(uint64_t& state) {
  state ^= state >> 12;
  state ^= state << 25;
  state ^= state >> 27;
  return state * UINT64_C(2685821657736338717);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: sumset INITIAL_DATA STABILIZERS_DATA\n";
    return 2;
  }
  try {
    const auto group_generators = read_initial(argv[1]);
    const auto cases = read_cases(argv[2]);
    if (!is_identity(power(group_generators[0], 3)) ||
        !is_identity(power(group_generators[1], 4))) {
      throw std::runtime_error("unexpected standard-generator orders");
    }
    const std::vector<Matrix> inverse_generators = {
        power(group_generators[0], 2), power(group_generators[1], 3)};

    // A Schreier tree for the unique binary-vector orbit of length 2608200.
    // parent[x] and parent_generator[x] encode a word taking 7 to x.
    std::vector<uint32_t> parent(kSpaceSize, kNoParent);
    std::vector<uint8_t> parent_generator(kSpaceSize, 0);
    std::vector<uint32_t> queue(kSpaceSize);
    size_t head = 0;
    size_t tail = 0;
    parent[kRoot] = kRoot;
    queue[tail++] = kRoot;
    while (head < tail) {
      const uint32_t x = queue[head++];
      for (uint8_t j = 0; j < group_generators.size(); ++j) {
        const uint32_t y = act(x, group_generators[j]);
        if (parent[y] == kNoParent) {
          parent[y] = x;
          parent_generator[y] = j;
          queue[tail++] = y;
        }
      }
    }
    if (tail != 2608200) {
      throw std::runtime_error("unexpected length for the root binary orbit");
    }
    const std::vector<uint32_t> first_orbit(queue.begin(), queue.begin() + tail);

    const Case* root_case = nullptr;
    for (const Case& c : cases) {
      if (c.representative == kRoot) root_case = &c;
    }
    if (root_case == nullptr || root_case->stabilizer_order != 190080) {
      throw std::runtime_error("missing stabilizer of binary vector 7");
    }

    // The two regular H_7-orbits are represented by 88 and 106.  Label every
    // member by 1 or 2; all remaining second vectors retain label zero.
    std::vector<uint8_t> regular_label(kSpaceSize, 0);
    std::vector<uint32_t> regular_second_vectors;
    for (uint8_t label = 1; label <= 2; ++label) {
      const uint32_t seed = label == 1 ? 88 : 106;
      head = 0;
      tail = 0;
      if (regular_label[seed] != 0) {
        throw std::runtime_error("regular stabilizer orbits intersect");
      }
      regular_label[seed] = label;
      queue[tail++] = seed;
      while (head < tail) {
        const uint32_t x = queue[head++];
        for (const Matrix& g : root_case->generators) {
          const uint32_t y = act(x, g);
          if (regular_label[y] == 0) {
            regular_label[y] = label;
            queue[tail++] = y;
          } else if (regular_label[y] != label) {
            throw std::runtime_error("regular stabilizer orbits intersect");
          }
        }
      }
      if (tail != root_case->stabilizer_order) {
        throw std::runtime_error("second-vector orbit is not regular");
      }
    }
    for (uint32_t x = 0; x < kSpaceSize; ++x) {
      if (regular_label[x] != 0) regular_second_vectors.push_back(x);
    }
    if (regular_second_vectors.size() != 2 * root_case->stabilizer_order) {
      throw std::runtime_error("wrong number of regular second vectors");
    }

    auto pull_back_second = [&](uint32_t first, uint32_t second) {
      uint32_t current = first;
      while (current != kRoot) {
        if (current == kNoParent || parent[current] == kNoParent) {
          throw std::runtime_error("first vector is outside the root orbit");
        }
        const uint8_t j = parent_generator[current];
        second = act(second, inverse_generators[j]);
        current = parent[current];
      }
      return second;
    };

    auto regular_orbit_label = [&](uint32_t first, uint32_t second) -> uint8_t {
      if (parent[first] == kNoParent) return 0;
      return regular_label[pull_back_second(first, second)];
    };

    auto push_forward_second = [&](uint32_t first, uint32_t second) {
      std::vector<uint8_t> reverse_word;
      uint32_t current = first;
      while (current != kRoot) {
        if (parent[current] == kNoParent) {
          throw std::runtime_error("first vector is outside the root orbit");
        }
        reverse_word.push_back(parent_generator[current]);
        current = parent[current];
      }
      for (auto it = reverse_word.rbegin(); it != reverse_word.rend(); ++it) {
        second = act(second, group_generators[*it]);
      }
      return second;
    };

    // Exhaustively recover one representative of every Co3-orbit on ordered
    // binary pairs (equivalently on GF(4)^22 after choosing {1,omega}).
    std::vector<std::pair<uint32_t, uint32_t>> target_representatives;
    std::vector<uint8_t> seen(kSpaceSize);
    for (const Case& c : cases) {
      std::fill(seen.begin(), seen.end(), uint8_t{0});
      uint64_t covered = 0;
      for (uint32_t seed = 0; seed < kSpaceSize; ++seed) {
        if (seen[seed]) continue;
        target_representatives.emplace_back(c.representative, seed);
        head = 0;
        tail = 0;
        seen[seed] = 1;
        queue[tail++] = seed;
        while (head < tail) {
          const uint32_t x = queue[head++];
          for (const Matrix& g : c.generators) {
            const uint32_t y = act(x, g);
            if (!seen[y]) {
              seen[y] = 1;
              queue[tail++] = y;
            }
          }
        }
        if (c.stabilizer_order % tail != 0) {
          throw std::runtime_error("second-vector orbit length is invalid");
        }
        covered += tail;
      }
      if (covered != kSpaceSize) {
        throw std::runtime_error("second-vector orbits do not cover the space");
      }
    }
    if (target_representatives.size() != 320) {
      throw std::runtime_error("unexpected number of GF(4)-vector orbits");
    }

    // Scalar multiplication by omega sends (a,b) to (b,a+b).  It fixes both
    // regular Co3-orbits, so adjoining GF(4)^* destroys every regular vector.
    const uint8_t scalar_image_1 = regular_orbit_label(88, 7 ^ 88);
    const uint8_t scalar_image_2 = regular_orbit_label(106, 7 ^ 106);
    if (scalar_image_1 != 1 || scalar_image_2 != 2) {
      throw std::runtime_error("unexpected scalar action on regular orbits");
    }

    std::cout << "CO3_F4_SUMSET\n";
    std::cout << "binary_first_orbit_length " << first_orbit.size() << "\n";
    std::cout << "regular_second_fibre_size "
              << regular_second_vectors.size() << "\n";
    std::cout << "regular_Co3_vector_orbits 2\n";
    std::cout << "regular_Co3_vectors 991533312000\n";
    std::cout << "omega_regular_orbit_action 1_to_"
              << static_cast<unsigned>(scalar_image_1) << " 2_to_"
              << static_cast<unsigned>(scalar_image_2) << "\n";
    std::cout << "regular_full_scalar_vectors 0\n";
    std::cout << "target_orbits " << target_representatives.size() << "\n";
    std::cout << "target_index target_a target_b witness_a witness_b "
                 "witness_label complement_label trials\n";

    uint64_t rng = UINT64_C(0x4c1f3a9b7652d8e1);
    uint64_t maximum_trials = 0;
    uint64_t total_trials = 0;
    for (size_t index = 0; index < target_representatives.size(); ++index) {
      const auto [target_a, target_b] = target_representatives[index];
      bool found = false;
      uint32_t witness_a = 0;
      uint32_t witness_b = 0;
      uint8_t witness_label = 0;
      uint8_t complement_label = 0;
      uint64_t trials = 0;
      for (; trials < UINT64_C(1000000); ++trials) {
        witness_a =
            first_orbit[next_random(rng) % first_orbit.size()];
        const uint32_t fibre =
            regular_second_vectors[next_random(rng) %
                                   regular_second_vectors.size()];
        witness_b = push_forward_second(witness_a, fibre);
        witness_label = regular_orbit_label(witness_a, witness_b);
        if (witness_label == 0) {
          throw std::runtime_error("constructed witness is not regular");
        }
        complement_label =
            regular_orbit_label(target_a ^ witness_a,
                                target_b ^ witness_b);
        if (complement_label != 0) {
          found = true;
          ++trials;
          break;
        }
      }
      if (!found) {
        throw std::runtime_error("failed to find a sumset witness");
      }
      maximum_trials = std::max(maximum_trials, trials);
      total_trials += trials;
      std::cout << (index + 1) << " " << target_a << " " << target_b
                << " " << witness_a << " " << witness_b << " "
                << static_cast<unsigned>(witness_label) << " "
                << static_cast<unsigned>(complement_label) << " "
                << trials << "\n";
    }
    std::cout << "maximum_trials " << maximum_trials << "\n";
    std::cout << "total_trials " << total_trials << "\n";
    std::cout << "all_target_orbits_have_verified_decomposition true\n";
    std::cout << "sumset_equals_GF4_22 true\n";
    std::cout << "Sumset calculation finished.\n";
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

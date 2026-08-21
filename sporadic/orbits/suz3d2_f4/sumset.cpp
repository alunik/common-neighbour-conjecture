// Find regular-plus-regular decompositions for every 3.Suz.2-orbit on
// GF(4)^24, both before and after adjoining the nontrivial field scalars.

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

constexpr unsigned kDim = 24;
constexpr std::uint32_t kSpaceSize = std::uint32_t{1} << kDim;
constexpr std::uint32_t kMask = kSpaceSize - 1;
constexpr std::uint64_t kGroupOrder = UINT64_C(2690072985600);
constexpr std::uint32_t kNoParent =
    std::numeric_limits<std::uint32_t>::max();

struct Matrix {
  std::uint32_t row[kDim]{};
};

struct Case {
  std::uint32_t representative = 0;
  std::uint64_t first_orbit_length = 0;
  std::uint64_t stabilizer_order = 0;
  std::vector<Matrix> generators;
};

struct RegularOrbit {
  std::uint8_t label = 0;
  std::size_t case_index = 0;
  std::uint32_t second_representative = 0;
};

std::uint32_t Act(std::uint32_t vector, const Matrix& matrix) {
  std::uint32_t image = 0;
  while (vector != 0) {
    const unsigned position =
        static_cast<unsigned>(__builtin_ctz(vector));
    image ^= matrix.row[position];
    vector &= vector - 1;
  }
  return image;
}

Matrix IdentityMatrix() {
  Matrix result;
  for (unsigned position = 0; position < kDim; ++position) {
    result.row[position] = std::uint32_t{1} << position;
  }
  return result;
}

Matrix Multiply(const Matrix& left, const Matrix& right) {
  Matrix result;
  for (unsigned position = 0; position < kDim; ++position) {
    result.row[position] = Act(left.row[position], right);
  }
  return result;
}

Matrix Power(Matrix base, unsigned exponent) {
  Matrix result = IdentityMatrix();
  while (exponent != 0) {
    if ((exponent & 1U) != 0) {
      result = Multiply(result, base);
    }
    base = Multiply(base, base);
    exponent >>= 1U;
  }
  return result;
}

bool IsIdentity(const Matrix& matrix) {
  for (unsigned position = 0; position < kDim; ++position) {
    if (matrix.row[position] != (std::uint32_t{1} << position)) {
      return false;
    }
  }
  return true;
}

std::vector<Matrix> ReadInitial(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open " + path);
  }
  std::vector<Matrix> matrices;
  std::string token;
  while (input >> token) {
    if (token != "MATRIX") {
      continue;
    }
    Matrix matrix;
    for (unsigned position = 0; position < kDim; ++position) {
      std::uint64_t row = 0;
      if (!(input >> row) || row >= kSpaceSize) {
        throw std::runtime_error("invalid packed initial matrix row");
      }
      matrix.row[position] = static_cast<std::uint32_t>(row);
    }
    matrices.push_back(matrix);
  }
  if (matrices.size() != 2) {
    throw std::runtime_error("expected two standard generators");
  }
  return matrices;
}

std::vector<Case> ReadCases(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open " + path);
  }
  std::vector<Case> cases;
  std::string token;
  while (input >> token) {
    if (token != "CASE") {
      continue;
    }
    Case item;
    std::size_t generator_count = 0;
    std::uint64_t support_size = 0;
    if (!(input >> item.representative >> item.first_orbit_length >>
          item.stabilizer_order >> generator_count >> support_size)) {
      throw std::runtime_error("invalid CASE header");
    }
    if (generator_count == 0) {
      throw std::runtime_error("empty stabilizer generating set");
    }
    for (std::size_t index = 0; index < generator_count; ++index) {
      if (!(input >> token) || token != "MATRIX") {
        throw std::runtime_error("missing MATRIX record");
      }
      Matrix matrix;
      for (unsigned position = 0; position < kDim; ++position) {
        std::uint64_t row = 0;
        if (!(input >> row) || row >= kSpaceSize) {
          throw std::runtime_error("invalid packed stabilizer matrix row");
        }
        matrix.row[position] = static_cast<std::uint32_t>(row);
      }
      item.generators.push_back(matrix);
    }
    if (!(input >> token) || token != "ENDCASE") {
      throw std::runtime_error("missing ENDCASE record");
    }
    cases.push_back(std::move(item));
  }
  if (cases.size() != 5) {
    throw std::runtime_error("expected five first-vector stabilizers");
  }
  return cases;
}

std::uint64_t NextRandom(std::uint64_t& state) {
  state ^= state >> 12;
  state ^= state << 25;
  state ^= state >> 27;
  return state * UINT64_C(2685821657736338717);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr
        << "usage: sumset INITIAL_DATA STABILIZERS_DATA "
           "WITNESSES_TSV\n";
    return 2;
  }
  try {
    const std::vector<Matrix> group_generators = ReadInitial(argv[1]);
    const std::vector<Case> cases = ReadCases(argv[2]);
    if (!IsIdentity(Power(group_generators[0], 2)) ||
        !IsIdentity(Power(group_generators[1], 3))) {
      throw std::runtime_error("unexpected standard-generator orders");
    }
    const std::vector<Matrix> inverse_generators = {
        group_generators[0], Power(group_generators[1], 2)};

    // Build a Schreier forest for all five binary-vector orbits.  For every
    // x, parent[x] and parent_generator[x] record the last edge of a word from
    // the case representative to x.
    std::vector<std::uint32_t> parent(kSpaceSize, kNoParent);
    std::vector<std::uint8_t> parent_generator(kSpaceSize, 0);
    std::vector<std::vector<std::uint32_t>> first_orbits(cases.size());
    std::vector<std::uint32_t> queue(kSpaceSize);
    std::size_t total_first_vectors = 0;
    for (std::size_t case_index = 0; case_index < cases.size();
         ++case_index) {
      const Case& item = cases[case_index];
      if (parent[item.representative] != kNoParent) {
        throw std::runtime_error("binary first-orbit roots intersect");
      }
      std::size_t head = 0;
      std::size_t tail = 0;
      parent[item.representative] = item.representative;
      queue[tail++] = item.representative;
      while (head < tail) {
        const std::uint32_t point = queue[head++];
        for (std::uint8_t generator_index = 0;
             generator_index < group_generators.size();
             ++generator_index) {
          const std::uint32_t image =
              Act(point, group_generators[generator_index]);
          if (parent[image] == kNoParent) {
            parent[image] = point;
            parent_generator[image] = generator_index;
            queue[tail++] = image;
          }
        }
      }
      if (tail != item.first_orbit_length) {
        throw std::runtime_error("unexpected binary first-orbit length");
      }
      first_orbits[case_index].assign(queue.begin(), queue.begin() + tail);
      total_first_vectors += tail;
    }
    if (total_first_vectors != kSpaceSize ||
        std::find(parent.begin(), parent.end(), kNoParent) != parent.end()) {
      throw std::runtime_error(
          "binary first-orbit forest does not cover the space");
    }

    auto CaseIndexFromRoot = [&](std::uint32_t root) -> std::size_t {
      for (std::size_t case_index = 0; case_index < cases.size();
           ++case_index) {
        if (cases[case_index].representative == root) {
          return case_index;
        }
      }
      throw std::runtime_error("unknown binary-orbit root");
    };

    auto NormalizeSecond =
        [&](std::uint32_t first,
            std::uint32_t second) -> std::pair<std::size_t, std::uint32_t> {
      std::uint32_t current = first;
      while (parent[current] != current) {
        if (parent[current] == kNoParent) {
          throw std::runtime_error("first vector is outside Schreier forest");
        }
        const std::uint8_t generator_index = parent_generator[current];
        second = Act(second, inverse_generators[generator_index]);
        current = parent[current];
      }
      return {CaseIndexFromRoot(current), second};
    };

    auto PushForwardSecond =
        [&](std::uint32_t first, std::uint32_t root_second) {
      std::vector<std::uint8_t> reverse_word;
      std::uint32_t current = first;
      while (parent[current] != current) {
        if (parent[current] == kNoParent) {
          throw std::runtime_error("first vector is outside Schreier forest");
        }
        reverse_word.push_back(parent_generator[current]);
        current = parent[current];
      }
      for (auto position = reverse_word.rbegin();
           position != reverse_word.rend(); ++position) {
        root_second = Act(root_second, group_generators[*position]);
      }
      return root_second;
    };

    // Enumerate every stabilizer orbit on the second component.  This both
    // constructs a complete G-orbit transversal on GF(4)^24 and labels every
    // regular fibre orbit for exact membership tests.
    std::vector<std::vector<std::uint8_t>> regular_fibre_label(
        cases.size(), std::vector<std::uint8_t>(kSpaceSize, 0));
    std::vector<std::uint8_t> seen(kSpaceSize);
    std::vector<std::pair<std::uint32_t, std::uint32_t>>
        target_representatives;
    std::vector<RegularOrbit> regular_orbits;
    std::uint64_t weighted_pair_coverage = 0;
    for (std::size_t case_index = 0; case_index < cases.size();
         ++case_index) {
      const Case& item = cases[case_index];
      std::fill(seen.begin(), seen.end(), std::uint8_t{0});
      std::uint64_t covered = 0;
      for (std::uint32_t seed = 0; seed < kSpaceSize; ++seed) {
        if (seen[seed] != 0) {
          continue;
        }
        target_representatives.emplace_back(item.representative, seed);
        std::size_t head = 0;
        std::size_t tail = 0;
        seen[seed] = 1;
        queue[tail++] = seed;
        while (head < tail) {
          const std::uint32_t point = queue[head++];
          for (const Matrix& generator : item.generators) {
            const std::uint32_t image = Act(point, generator);
            if (seen[image] == 0) {
              seen[image] = 1;
              queue[tail++] = image;
            }
          }
        }
        if (item.stabilizer_order % tail != 0) {
          throw std::runtime_error(
              "second orbit length does not divide stabilizer order");
        }
        covered += tail;
        weighted_pair_coverage += item.first_orbit_length * tail;
        if (tail == item.stabilizer_order) {
          if (regular_orbits.size() >= 254) {
            throw std::runtime_error("too many regular orbits for labels");
          }
          const std::uint8_t label =
              static_cast<std::uint8_t>(regular_orbits.size() + 1);
          regular_orbits.push_back(
              RegularOrbit{label, case_index, seed});
          for (std::size_t position = 0; position < tail; ++position) {
            if (regular_fibre_label[case_index][queue[position]] != 0) {
              throw std::runtime_error("regular fibre orbits intersect");
            }
            regular_fibre_label[case_index][queue[position]] = label;
          }
        }
      }
      if (covered != kSpaceSize) {
        throw std::runtime_error(
            "second-component stabilizer orbits do not cover the space");
      }
    }
    if (weighted_pair_coverage !=
        std::uint64_t{kSpaceSize} * kSpaceSize) {
      throw std::runtime_error(
          "ordered-pair orbits do not cover GF(4)^24");
    }
    if (regular_orbits.empty()) {
      throw std::runtime_error("3.Suz.2 has no regular vectors");
    }

    auto RegularOrbitLabel =
        [&](std::uint32_t first, std::uint32_t second) -> std::uint8_t {
      const auto [case_index, normalized_second] =
          NormalizeSecond(first, second);
      return regular_fibre_label[case_index][normalized_second];
    };

    // Multiplication by omega in the basis {1,omega}, omega^2=omega+1, is
    // (a,b) -> (b,a+b).  Its permutation on the regular G-orbits determines
    // exactly which vectors remain regular after adjoining external GF(4)^*.
    std::vector<std::uint8_t> omega_image(regular_orbits.size() + 1, 0);
    std::vector<std::uint8_t> omega_preimages(
        regular_orbits.size() + 1, 0);
    for (const RegularOrbit& orbit : regular_orbits) {
      const std::uint32_t first =
          cases[orbit.case_index].representative;
      const std::uint32_t second = orbit.second_representative;
      const std::uint8_t image =
          RegularOrbitLabel(second, first ^ second);
      if (image == 0 || image > regular_orbits.size()) {
        throw std::runtime_error(
            "external scalar does not permute regular G-orbits");
      }
      omega_image[orbit.label] = image;
      ++omega_preimages[image];
    }
    std::vector<std::uint8_t> g_labels;
    std::vector<std::uint8_t> hmax_labels;
    for (std::uint8_t label = 1; label <= regular_orbits.size(); ++label) {
      if (omega_preimages[label] != 1 ||
          omega_image[omega_image[omega_image[label]]] != label) {
        throw std::runtime_error(
            "external scalar labels do not form an order-three permutation");
      }
      g_labels.push_back(label);
      if (omega_image[label] != label) {
        hmax_labels.push_back(label);
      }
    }
    if (hmax_labels.size() % 3 != 0) {
      throw std::runtime_error(
          "external scalar extension has an incomplete regular orbit cycle");
    }

    std::ofstream witnesses(argv[3]);
    if (!witnesses) {
      throw std::runtime_error("cannot create witness table");
    }
    witnesses << "group\ttarget_index\ttarget_a\ttarget_b\twitness_a"
                 "\twitness_b\tcomplement_a\tcomplement_b\twitness_label"
                 "\tcomplement_label\ttrials\n";

    std::uint64_t rng = UINT64_C(0x7410d3e269abc85f);
    auto RunSumset =
        [&](const std::string& group_name,
            const std::vector<std::uint8_t>& allowed_labels) {
      std::vector<std::uint8_t> allowed(
          regular_orbits.size() + 1, 0);
      for (const std::uint8_t label : allowed_labels) {
        allowed[label] = 1;
      }
      std::uint64_t maximum_trials = 0;
      std::uint64_t total_trials = 0;
      for (std::size_t target_index = 0;
           target_index < target_representatives.size(); ++target_index) {
        const auto [target_a, target_b] =
            target_representatives[target_index];
        bool found = false;
        std::uint32_t witness_a = 0;
        std::uint32_t witness_b = 0;
        std::uint8_t witness_label = 0;
        std::uint8_t complement_label = 0;
        std::uint64_t trials = 0;
        for (; trials < UINT64_C(10000000); ++trials) {
          const std::uint8_t chosen_label =
              allowed_labels[
                  NextRandom(rng) % allowed_labels.size()];
          const RegularOrbit& chosen_orbit =
              regular_orbits[chosen_label - 1];
          const std::size_t case_index = chosen_orbit.case_index;
          const auto& first_orbit = first_orbits[case_index];
          witness_a =
              first_orbit[NextRandom(rng) % first_orbit.size()];
          std::uint32_t root_second = 0;
          do {
            root_second =
                static_cast<std::uint32_t>(NextRandom(rng)) & kMask;
          } while (
              regular_fibre_label[case_index][root_second] !=
              chosen_label);
          witness_b = PushForwardSecond(witness_a, root_second);
          witness_label = RegularOrbitLabel(witness_a, witness_b);
          if (witness_label != chosen_label ||
              allowed[witness_label] == 0) {
            throw std::runtime_error(
                "constructed witness is not in the requested regular set");
          }
          complement_label =
              RegularOrbitLabel(target_a ^ witness_a,
                                target_b ^ witness_b);
          if (complement_label != 0 &&
              allowed[complement_label] != 0) {
            found = true;
            ++trials;
            break;
          }
        }
        if (!found) {
          throw std::runtime_error(
              group_name + " sumset witness search exhausted");
        }
        const std::uint32_t complement_a = target_a ^ witness_a;
        const std::uint32_t complement_b = target_b ^ witness_b;
        if ((witness_a ^ complement_a) != target_a ||
            (witness_b ^ complement_b) != target_b ||
            RegularOrbitLabel(complement_a, complement_b) !=
                complement_label) {
          throw std::runtime_error("sumset witness verification failed");
        }
        maximum_trials = std::max(maximum_trials, trials);
        total_trials += trials;
        witnesses << group_name << '\t' << (target_index + 1) << '\t'
                  << target_a << '\t' << target_b << '\t' << witness_a
                  << '\t' << witness_b << '\t' << complement_a << '\t'
                  << complement_b << '\t'
                  << static_cast<unsigned>(witness_label) << '\t'
                  << static_cast<unsigned>(complement_label) << '\t'
                  << trials << '\n';
      }
      std::cout << group_name << "_allowed_G_orbit_labels";
      for (const std::uint8_t label : allowed_labels) {
        std::cout << " " << static_cast<unsigned>(label);
      }
      std::cout << "\n";
      std::cout << group_name << "_maximum_trials " << maximum_trials
                << "\n";
      std::cout << group_name << "_total_trials " << total_trials
                << "\n";
      std::cout << group_name
                << "_all_target_orbits_have_verified_decomposition true\n";
      std::cout << group_name << "_sumset_equals_GF4_24 true\n";
    };

    std::cout << "SUZ3D2_F4_SUMSET\n";
    std::cout << "binary_first_orbits " << cases.size() << "\n";
    std::cout << "binary_first_vectors_covered " << total_first_vectors
              << "\n";
    std::cout << "GF4_vector_orbits "
              << target_representatives.size() << "\n";
    std::cout << "regular_G_vector_orbits " << regular_orbits.size()
              << "\n";
    std::cout << "regular_G_vectors "
              << regular_orbits.size() * kGroupOrder << "\n";
    std::cout << "omega_regular_orbit_action";
    for (std::uint8_t label = 1; label <= regular_orbits.size(); ++label) {
      std::cout << " " << static_cast<unsigned>(label) << "->"
                << static_cast<unsigned>(omega_image[label]);
    }
    std::cout << "\n";
    std::cout << "regular_Hmax_G_orbit_components "
              << hmax_labels.size() << "\n";
    std::cout << "regular_Hmax_vector_orbits "
              << hmax_labels.size() / 3 << "\n";
    std::cout << "regular_Hmax_vectors "
              << hmax_labels.size() * kGroupOrder << "\n";
    std::cout << "Hmax_has_regular_vectors "
              << (hmax_labels.empty() ? "false" : "true") << "\n";
    std::cout << "target_orbits "
              << target_representatives.size() << "\n";
    RunSumset("G", g_labels);
    if (hmax_labels.empty()) {
      std::cout << "Hmax_affine_base_size_greater_than_two true\n";
      std::cout << "Hmax_sumset_test_not_applicable true\n";
    } else {
      RunSumset("Hmax", hmax_labels);
    }
    std::cout << "Sumset calculation finished.\n";
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
  return 0;
}

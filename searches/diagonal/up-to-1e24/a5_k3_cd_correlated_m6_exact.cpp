// Exact common-neighbour witnesses for the residual correlated A5,k=3,m=6
// quotient classes.  Magma supplies one representative of every quotient
// orbit by a double-coset calculation; this engine independently rebuilds
// the quotient group and the 55-colour/77-target incidence relation.

#define A5_CORRELATED_M3_EMBEDDED
#include "../up-to-1e18/a5_k3_cd_correlated_m3_exact.cpp"
#undef A5_CORRELATED_M3_EMBEDDED

#include <bit>
#include <atomic>
#include <charconv>
#include <cstdlib>
#include <limits>
#include <memory>
#include <thread>

namespace {

using Tuple6 = std::array<std::uint8_t, 6>;

struct W6 {
  std::array<std::uint8_t, 6> components{};
  Tuple6 top{};
};

struct Case6 {
  int number = 0;
  int m = 0;
  int order = 0;
  std::uint64_t regular_points = 0;
  std::vector<W6> generators;
  std::vector<Tuple6> x_reps;
  std::vector<Tuple6> y_reps;
};

std::string field_value(const std::string& line, const std::string& key) {
  const std::string needle = "|" + key + "=";
  const auto begin = line.find(needle);
  if (begin == std::string::npos) return {};
  const auto value_begin = begin + needle.size();
  const auto end = line.find('|', value_begin);
  return line.substr(value_begin, end == std::string::npos
                                      ? std::string::npos : end - value_begin);
}

std::uint64_t unsigned_field(const std::string& line, const std::string& key) {
  const std::string value = field_value(line, key);
  if (value.empty()) throw std::runtime_error("missing field " + key);
  std::uint64_t answer = 0;
  const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), answer);
  if (error != std::errc{} || end != value.data() + value.size()) {
    throw std::runtime_error("bad integer field " + key);
  }
  return answer;
}

Tuple6 tuple6_from_numbers(const std::string& text, int width, int offset = 0) {
  const auto values = parse_numbers(text);
  if (values.size() != static_cast<std::size_t>(width + offset)) {
    throw std::runtime_error("bad tuple descriptor");
  }
  Tuple6 answer{0,1,2,3,4,5};
  for (int i = 0; i < width; ++i) answer[i] = values[i + offset] - 1;
  return answer;
}

std::vector<Case6> parse_representatives(
    const std::string& path, const CoreData& core) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open orbit-representative transcript");
  std::map<int, Case6> cases;
  std::string line;
  int width = 0;
  while (std::getline(input, line)) {
    if (line.starts_with("ORBIT_REP_SHAPE|")) {
      width = unsigned_field(line, "m");
      if (width < 3 || width > 6) throw std::runtime_error("bad compound width");
    } else if (line.starts_with("ORBIT_REP_CASE|")) {
      if (width == 0) throw std::runtime_error("case before shape");
      Case6 item;
      item.number = unsigned_field(line, "case");
      item.m = width;
      item.order = unsigned_field(line, "M");
      item.regular_points = unsigned_field(line, "regular");
      cases[item.number] = item;
    } else if (line.starts_with("ORBIT_REP_GENERATOR|")) {
      const int number = unsigned_field(line, "case");
      if (!cases.contains(number)) throw std::runtime_error("generator before case");
      const auto top_begin = line.find("|top=");
      const auto component_begin = line.find("|components=");
      if (top_begin == std::string::npos || component_begin == std::string::npos) {
        throw std::runtime_error("bad generator row");
      }
      W6 generator;
      generator.top = tuple6_from_numbers(
          line.substr(top_begin, component_begin - top_begin), width);
      generator.components.fill(core.s3_index.at(SmallPerm{0,1,2}));
      std::string components = line.substr(component_begin + 12);
      if (components.front() != '[' || components.back() != ']') {
        throw std::runtime_error("bad components bracket");
      }
      components = components.substr(1, components.size() - 2);
      std::stringstream stream(components);
      std::string component;
      int coordinate = 0;
      while (std::getline(stream, component, ';')) {
        const auto values = parse_numbers(component);
        if (values.size() != 4 || coordinate >= width) {
          throw std::runtime_error("bad component descriptor");
        }
        SmallPerm permutation{};
        for (int i = 0; i < 3; ++i) permutation[i] = values[i + 1] - 1;
        generator.components[coordinate++] =
            values[0] * 6 + core.s3_index.at(permutation);
      }
      if (coordinate != width) throw std::runtime_error("bad component count");
      cases[number].generators.push_back(generator);
    } else if (line.starts_with("X_REGULAR_REP|")) {
      const int number = unsigned_field(line, "case");
      cases.at(number).x_reps.push_back(
          tuple6_from_numbers(field_value(line, "tuple"), width));
    } else if (line.starts_with("Y_TARGET_REP|")) {
      const int number = unsigned_field(line, "case");
      cases.at(number).y_reps.push_back(
          tuple6_from_numbers(field_value(line, "tuple"), width));
    } else if (line.starts_with("ORBIT_REP_CASE_COMPLETE|")) {
      const int number = unsigned_field(line, "case");
      const auto& item = cases.at(number);
      if (item.x_reps.size() != unsigned_field(line, "x") ||
          item.y_reps.size() != unsigned_field(line, "y")) {
        throw std::runtime_error("representative census mismatch");
      }
    }
  }
  std::vector<Case6> answer;
  for (auto& [number, item] : cases) {
    if (item.generators.empty() || item.x_reps.empty() || item.y_reps.empty()) {
      throw std::runtime_error("incomplete residual case");
    }
    answer.push_back(std::move(item));
  }
  return answer;
}

Tuple6 compose6(const Tuple6& left, const Tuple6& right, int width) {
  Tuple6 answer{};
  for (int i = 0; i < 6; ++i) answer[i] = i;
  for (int i = 0; i < width; ++i) answer[i] = right[left[i]];
  return answer;
}

W6 multiply6(
    const W6& left, const W6& right, const CoreData& core, int width) {
  W6 answer;
  answer.top = compose6(left.top, right.top, width);
  answer.components.fill(core.s3_index.at(SmallPerm{0,1,2}));
  for (int i = 0; i < width; ++i) {
    answer.components[i] = r_product(
        left.components[i], right.components[left.top[i]], core);
  }
  return answer;
}

std::uint64_t key6(const W6& element, int width) {
  std::uint64_t key = width;
  for (int i = 0; i < width; ++i) key = key * 6 + element.top[i];
  for (int i = 0; i < width; ++i) key = key * 12 + element.components[i];
  return key;
}

std::vector<W6> closure6(const Case6& item, const CoreData& core) {
  W6 identity;
  identity.top = {0,1,2,3,4,5};
  identity.components.fill(core.s3_index.at(SmallPerm{0,1,2}));
  std::vector<W6> elements = {identity};
  std::unordered_map<std::uint64_t, int> seen;
  seen.reserve(static_cast<std::size_t>(item.order) * 5 / 4);
  seen.emplace(key6(identity, item.m), 0);
  for (std::size_t head = 0; head < elements.size(); ++head) {
    for (const auto& generator : item.generators) {
      const W6 product = multiply6(elements[head], generator, core, item.m);
      const auto [it, inserted] = seen.emplace(key6(product, item.m), elements.size());
      if (inserted) elements.push_back(product);
    }
  }
  if (elements.size() != static_cast<std::size_t>(item.order)) {
    throw std::runtime_error("quotient closure order mismatch");
  }
  return elements;
}

Tuple6 act6(
    const Tuple6& input, const W6& element, const CoreData& core, int width) {
  Tuple6 output{};
  for (int i = 0; i < width; ++i) {
    output[element.top[i]] = core.x_action[element.components[i]][input[i]];
  }
  return output;
}

constexpr std::array<std::size_t, 6> powers = {
    1, 55, 55*55, 55*55*55, 55ULL*55*55*55,
    55ULL*55*55*55*55};

constexpr std::size_t triple_codes = 55ULL * 55 * 55;
constexpr std::uint64_t private_dense_build_threshold = 400'000'000ULL;
constexpr std::uint64_t bidirectional_sparse_threshold = 1'000'000'000ULL;
// The largest residual full-frontier case has 9,825,753,600 regular points.
// A forward-only 32-bit CSR for it occupies about 39.3 GB and avoids the much
// slower dense Cartesian-box scan over hundreds of thousands of targets.
constexpr std::uint64_t sparse_edge_threshold = 10'000'000'000ULL;

struct SparseRegularIndex {
  std::vector<std::uint64_t> forward_offsets;
  std::vector<std::uint64_t> reverse_offsets;
  std::vector<std::uint32_t> forward_edges;
  std::vector<std::uint32_t> reverse_edges;
};

struct RegularSet6 {
  std::vector<std::uint64_t> masks;
  std::unique_ptr<SparseRegularIndex> sparse;
};

std::vector<std::uint64_t> regular_prefix_masks(
    const Case6& item, const CoreData& core, const std::vector<W6>& group,
    int thread_count) {
  const bool private_dense_build =
      item.regular_points > private_dense_build_threshold && thread_count > 1;
  const int builder_threads = private_dense_build
      ? std::min(thread_count, 8) : thread_count;
  std::vector<std::uint64_t> masks;
  std::vector<std::vector<std::uint64_t>> private_masks;
  if (private_dense_build) {
    private_masks.resize(builder_threads);
    for (auto& private_mask : private_masks) {
      private_mask.resize(powers[item.m - 1], 0);
    }
  } else {
    masks.resize(powers[item.m - 1], 0);
  }
  std::atomic<std::size_t> next_representative{0};
  auto worker = [&](int thread) {
    while (true) {
      const std::size_t number = next_representative.fetch_add(
          1, std::memory_order_relaxed);
      if (number >= item.x_reps.size()) break;
      const auto& representative = item.x_reps[number];
      for (const auto& element : group) {
        const Tuple6 image = act6(representative, element, core, item.m);
        std::size_t prefix = 0;
        for (int i = 0; i < item.m - 1; ++i) {
          prefix += image[i] * powers[i];
        }
        const std::uint64_t bit = std::uint64_t{1} << image[item.m - 1];
        if (private_dense_build) {
          private_masks[thread][prefix] |= bit;
        } else if (thread_count == 1) {
          masks[prefix] |= bit;
        } else {
          __atomic_fetch_or(&masks[prefix], bit, __ATOMIC_RELAXED);
        }
      }
    }
  };
  std::vector<std::thread> threads;
  threads.reserve(std::max(thread_count, builder_threads));
  for (int thread = 0; thread < builder_threads; ++thread) {
    threads.emplace_back(worker, thread);
  }
  for (auto& thread : threads) thread.join();
  if (private_dense_build) {
    masks = std::move(private_masks[0]);
  }
  std::vector<std::uint64_t> partial_counts(thread_count, 0);
  threads.clear();
  for (int thread = 0; thread < thread_count; ++thread) {
    threads.emplace_back([&, thread]() {
      const std::size_t begin = masks.size() * thread / thread_count;
      const std::size_t end = masks.size() * (thread + 1) / thread_count;
      std::uint64_t count = 0;
      for (std::size_t index = begin; index < end; ++index) {
        if (private_dense_build) {
          for (int source = 1; source < builder_threads; ++source) {
            masks[index] |= private_masks[source][index];
          }
        }
        count += std::popcount(masks[index]);
      }
      partial_counts[thread] = count;
    });
  }
  for (auto& thread : threads) thread.join();
  std::uint64_t count = 0;
  for (std::uint64_t partial : partial_counts) count += partial;
  if (count != item.regular_points ||
      count != item.x_reps.size() * static_cast<std::uint64_t>(item.order)) {
    throw std::runtime_error("regular-set census mismatch: got " +
        std::to_string(count) + " expected " +
        std::to_string(item.regular_points));
  }
  return masks;
}

RegularSet6 build_regular_set(
    const Case6& item, const CoreData& core, const std::vector<W6>& group,
    int thread_count) {
  RegularSet6 answer;
  answer.masks = regular_prefix_masks(item, core, group, thread_count);
  if (item.regular_points > sparse_edge_threshold) return answer;

  auto sparse = std::make_unique<SparseRegularIndex>();
  const bool bidirectional =
      item.regular_points <= bidirectional_sparse_threshold;
  std::vector<std::uint64_t> forward_degrees(triple_codes, 0);
  std::vector<std::uint64_t> reverse_degrees(
      bidirectional ? triple_codes : 0, 0);
  if (!bidirectional) {
    std::vector<std::uint64_t> partial_counts(thread_count, 0);
    std::vector<std::thread> threads;
    threads.reserve(thread_count);
    for (int thread = 0; thread < thread_count; ++thread) {
      threads.emplace_back([&, thread]() {
        const std::size_t begin = triple_codes * thread / thread_count;
        const std::size_t end = triple_codes * (thread + 1) / thread_count;
        std::uint64_t local_count = 0;
        for (std::size_t left = begin; left < end; ++left) {
          std::uint64_t degree = 0;
          for (std::size_t upper = 0; upper < 55 * 55; ++upper) {
            degree += std::popcount(
                answer.masks[left + triple_codes * upper]);
          }
          forward_degrees[left] = degree;
          local_count += degree;
        }
        partial_counts[thread] = local_count;
      });
    }
    for (auto& thread : threads) thread.join();
    std::uint64_t count = 0;
    for (std::uint64_t partial : partial_counts) count += partial;
    if (count != item.regular_points) {
      throw std::runtime_error("forward sparse regular-set census mismatch");
    }
    sparse->forward_offsets.resize(triple_codes + 1, 0);
    for (std::size_t left = 0; left < triple_codes; ++left) {
      sparse->forward_offsets[left + 1] =
          sparse->forward_offsets[left] + forward_degrees[left];
    }
    sparse->forward_edges.resize(item.regular_points);
    std::atomic<bool> row_census_mismatch(false);
    threads.clear();
    for (int thread = 0; thread < thread_count; ++thread) {
      threads.emplace_back([&, thread]() {
        const std::size_t begin = triple_codes * thread / thread_count;
        const std::size_t end = triple_codes * (thread + 1) / thread_count;
        for (std::size_t left = begin; left < end; ++left) {
          std::uint64_t cursor = sparse->forward_offsets[left];
          for (std::uint32_t upper = 0; upper < 55 * 55; ++upper) {
            std::uint64_t bits =
                answer.masks[left + triple_codes * upper];
            const std::uint32_t x3 = upper % 55;
            const std::uint32_t x4 = upper / 55;
            while (bits) {
              const std::uint32_t x5 = std::countr_zero(bits);
              bits &= bits - 1;
              sparse->forward_edges[cursor++] =
                  x3 + 55 * x4 + 55 * 55 * x5;
            }
          }
          if (cursor != sparse->forward_offsets[left + 1]) {
            row_census_mismatch.store(true, std::memory_order_relaxed);
          }
        }
      });
    }
    for (auto& thread : threads) thread.join();
    if (row_census_mismatch.load(std::memory_order_relaxed)) {
      throw std::runtime_error("forward sparse row census mismatch");
    }
    answer.masks.clear();
    answer.masks.shrink_to_fit();
    answer.sparse = std::move(sparse);
    return answer;
  }
  std::uint64_t count = 0;
  for (std::size_t prefix = 0; prefix < answer.masks.size(); ++prefix) {
    std::uint64_t bits = answer.masks[prefix];
    if (!bits) continue;
    const std::uint32_t left = prefix % triple_codes;
    const std::size_t upper = prefix / triple_codes;
    const std::uint32_t x3 = upper % 55;
    const std::uint32_t x4 = upper / 55;
    while (bits) {
      const std::uint32_t x5 = std::countr_zero(bits);
      bits &= bits - 1;
      const std::uint32_t right = x3 + 55 * x4 + 55 * 55 * x5;
      ++forward_degrees[left];
      if (bidirectional) ++reverse_degrees[right];
      ++count;
    }
  }
  if (count != item.regular_points) {
    throw std::runtime_error("sparse regular-set census mismatch");
  }
  auto make_offsets = [](const std::vector<std::uint64_t>& degrees) {
    std::vector<std::uint64_t> offsets(degrees.size() + 1, 0);
    for (std::size_t index = 0; index < degrees.size(); ++index) {
      offsets[index + 1] = offsets[index] + degrees[index];
    }
    return offsets;
  };
  sparse->forward_offsets = make_offsets(forward_degrees);
  if (bidirectional) {
    sparse->reverse_offsets = make_offsets(reverse_degrees);
  }
  sparse->forward_edges.resize(item.regular_points);
  if (bidirectional) sparse->reverse_edges.resize(item.regular_points);
  std::vector<std::uint64_t> forward_cursor = sparse->forward_offsets;
  std::vector<std::uint64_t> reverse_cursor = sparse->reverse_offsets;
  for (std::size_t prefix = 0; prefix < answer.masks.size(); ++prefix) {
    std::uint64_t bits = answer.masks[prefix];
    if (!bits) continue;
    const std::uint32_t left = prefix % triple_codes;
    const std::size_t upper = prefix / triple_codes;
    const std::uint32_t x3 = upper % 55;
    const std::uint32_t x4 = upper / 55;
    while (bits) {
      const std::uint32_t x5 = std::countr_zero(bits);
      bits &= bits - 1;
      const std::uint32_t right = x3 + 55 * x4 + 55 * 55 * x5;
      sparse->forward_edges[forward_cursor[left]++] = right;
      if (bidirectional) {
        sparse->reverse_edges[reverse_cursor[right]++] = left;
      }
    }
  }
  answer.masks.clear();
  answer.masks.shrink_to_fit();
  answer.sparse = std::move(sparse);
  return answer;
}

bool sparse_box_search(
    const std::array<std::uint64_t,6>& allowed,
    const SparseRegularIndex& sparse, Tuple6* witness,
    std::uint64_t* probes) {
  const std::uint64_t forward_work =
      std::popcount(allowed[0]) * std::popcount(allowed[1]) *
      std::popcount(allowed[2]);
  const std::uint64_t reverse_work = sparse.reverse_edges.empty()
      ? std::numeric_limits<std::uint64_t>::max() :
      std::popcount(allowed[3]) * std::popcount(allowed[4]) *
      std::popcount(allowed[5]);
  if (forward_work <= reverse_work) {
    std::uint64_t first = allowed[0];
    while (first) {
      const std::uint32_t x0 = std::countr_zero(first);
      first &= first - 1;
      std::uint64_t second = allowed[1];
      while (second) {
        const std::uint32_t x1 = std::countr_zero(second);
        second &= second - 1;
        std::uint64_t third = allowed[2];
        while (third) {
          const std::uint32_t x2 = std::countr_zero(third);
          third &= third - 1;
          const std::uint32_t left = x0 + 55 * x1 + 55 * 55 * x2;
          for (std::uint64_t edge = sparse.forward_offsets[left];
               edge < sparse.forward_offsets[left + 1]; ++edge) {
            ++*probes;
            const std::uint32_t right = sparse.forward_edges[edge];
            const std::uint32_t x3 = right % 55;
            const std::uint32_t x4 = (right / 55) % 55;
            const std::uint32_t x5 = right / (55 * 55);
            if (((allowed[3] >> x3) & 1) &&
                ((allowed[4] >> x4) & 1) &&
                ((allowed[5] >> x5) & 1)) {
              *witness = Tuple6{
                  static_cast<std::uint8_t>(x0),
                  static_cast<std::uint8_t>(x1),
                  static_cast<std::uint8_t>(x2),
                  static_cast<std::uint8_t>(x3),
                  static_cast<std::uint8_t>(x4),
                  static_cast<std::uint8_t>(x5)};
              return true;
            }
          }
        }
      }
    }
  } else {
    std::uint64_t fourth = allowed[3];
    while (fourth) {
      const std::uint32_t x3 = std::countr_zero(fourth);
      fourth &= fourth - 1;
      std::uint64_t fifth = allowed[4];
      while (fifth) {
        const std::uint32_t x4 = std::countr_zero(fifth);
        fifth &= fifth - 1;
        std::uint64_t sixth = allowed[5];
        while (sixth) {
          const std::uint32_t x5 = std::countr_zero(sixth);
          sixth &= sixth - 1;
          const std::uint32_t right = x3 + 55 * x4 + 55 * 55 * x5;
          for (std::uint64_t edge = sparse.reverse_offsets[right];
               edge < sparse.reverse_offsets[right + 1]; ++edge) {
            ++*probes;
            const std::uint32_t left = sparse.reverse_edges[edge];
            const std::uint32_t x0 = left % 55;
            const std::uint32_t x1 = (left / 55) % 55;
            const std::uint32_t x2 = left / (55 * 55);
            if (((allowed[0] >> x0) & 1) &&
                ((allowed[1] >> x1) & 1) &&
                ((allowed[2] >> x2) & 1)) {
              *witness = Tuple6{
                  static_cast<std::uint8_t>(x0),
                  static_cast<std::uint8_t>(x1),
                  static_cast<std::uint8_t>(x2),
                  static_cast<std::uint8_t>(x3),
                  static_cast<std::uint8_t>(x4),
                  static_cast<std::uint8_t>(x5)};
              return true;
            }
          }
        }
      }
    }
  }
  return false;
}

bool box_search_recursive(
    int depth, int width, const std::array<int,5>& order,
    const std::array<std::uint64_t,6>& allowed,
    std::size_t prefix, const std::vector<std::uint64_t>& regular,
    Tuple6* witness, std::uint64_t* probes) {
  if (depth == width - 1) {
    ++*probes;
    const std::uint64_t possible = regular[prefix] & allowed[width - 1];
    if (!possible) return false;
    witness->at(width - 1) = std::countr_zero(possible);
    return true;
  }
  const int coordinate = order[depth];
  std::uint64_t values = allowed[coordinate];
  while (values) {
    const int value = std::countr_zero(values);
    values &= values - 1;
    witness->at(coordinate) = value;
    if (box_search_recursive(depth + 1, width, order, allowed,
                             prefix + value * powers[coordinate], regular,
                             witness, probes)) return true;
  }
  return false;
}

bool target_witness(
    const Tuple6& target, const Case6& item, const CoreData& core,
    const std::vector<W6>& group, const RegularSet6& regular,
    Tuple6* left, Tuple6* right, std::size_t* left_index,
    std::uint64_t* probes, std::uint64_t* expanded_left_tests) {
  auto try_left = [&](const Tuple6& u, std::size_t index) {
    std::array<std::uint64_t,6> allowed{};
    for (int i = 0; i < item.m; ++i) allowed[i] = core.incidence[target[i]][u[i]];
    std::array<int,5> order = {0,1,2,3,4};
    std::sort(order.begin(), order.begin() + item.m - 1, [&](int a, int b) {
      return std::popcount(allowed[a]) < std::popcount(allowed[b]);
    });
    Tuple6 v{};
    const bool found = regular.sparse
        ? sparse_box_search(allowed, *regular.sparse, &v, probes)
        : box_search_recursive(
              0, item.m, order, allowed, 0, regular.masks, &v, probes);
    if (found) {
      *left = u;
      *right = v;
      *left_index = index;
      return true;
    }
    return false;
  };
  for (std::size_t index = 0; index < item.x_reps.size(); ++index) {
    if (try_left(item.x_reps[index], index)) return true;
  }
  // A target representative need not be aligned with the selected regular
  // source-orbit representative.  The fallback is exhaustive: act every
  // source representative by every quotient element.  It is normally entered
  // only for a tiny number of target orbits.
  for (std::size_t index = 0; index < item.x_reps.size(); ++index) {
    for (std::size_t element = 1; element < group.size(); ++element) {
      ++*expanded_left_tests;
      if (try_left(act6(item.x_reps[index], group[element], core, item.m), index)) {
        return true;
      }
    }
  }
  return false;
}

std::uint64_t mix_witness(
    std::uint64_t hash, const Tuple6& target, const Tuple6& left,
    const Tuple6& right, std::size_t left_index, int width) {
  auto mix = [&](std::uint64_t value) {
    hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
  };
  mix(left_index);
  for (int i = 0; i < width; ++i) { mix(target[i]); mix(left[i]); mix(right[i]); }
  return hash;
}

struct Witness6 {
  Tuple6 left{};
  Tuple6 right{};
  std::size_t left_index = 0;
  bool pass = false;
};

void audit6(
    const std::string& mode, const Case6& item, const CoreData& core,
    int thread_count) {
  const auto group = closure6(item, core);
  const auto regular = build_regular_set(
      item, core, group, thread_count);
  std::vector<Witness6> witnesses(item.y_reps.size());
  std::vector<std::uint64_t> thread_probes(thread_count, 0);
  std::vector<std::uint64_t> thread_expanded(thread_count, 0);
  std::atomic<std::size_t> next_target{0};
  auto worker = [&](int thread) {
    while (true) {
      constexpr std::size_t chunk = 64;
      const std::size_t begin = next_target.fetch_add(
          chunk, std::memory_order_relaxed);
      if (begin >= item.y_reps.size()) break;
      const std::size_t end = std::min(begin + chunk, item.y_reps.size());
      for (std::size_t target_number = begin; target_number < end;
           ++target_number) {
        Witness6& witness = witnesses[target_number];
        witness.pass = target_witness(
            item.y_reps[target_number], item, core, group, regular,
            &witness.left, &witness.right, &witness.left_index,
            &thread_probes[thread], &thread_expanded[thread]);
      }
    }
  };
  std::vector<std::thread> threads;
  threads.reserve(thread_count);
  for (int thread = 0; thread < thread_count; ++thread) {
    threads.emplace_back(worker, thread);
  }
  for (auto& thread : threads) thread.join();

  std::uint64_t probes = 0;
  std::uint64_t expanded_left_tests = 0;
  for (int thread = 0; thread < thread_count; ++thread) {
    probes += thread_probes[thread];
    expanded_left_tests += thread_expanded[thread];
  }
  std::uint64_t digest = 0xcbf29ce484222325ULL;
  std::size_t maximum_left_index = 0;
  for (std::size_t target_number = 0;
       target_number < item.y_reps.size(); ++target_number) {
    const Witness6& witness = witnesses[target_number];
    if (!witness.pass) {
      std::cout << "CORRELATED_M6_DEFECT|mode=" << mode
                << "|case=" << item.number
                << "|target_orbit=" << target_number + 1 << "\n";
      throw std::runtime_error("target orbit has no certified witness");
    }
    maximum_left_index = std::max(maximum_left_index, witness.left_index);
    digest = mix_witness(
        digest, item.y_reps[target_number], witness.left, witness.right,
        witness.left_index, item.m);
  }
  std::cout << "CORRELATED_M6_RESULT|mode=" << mode
            << "|case=" << item.number
            << "|m=" << item.m
            << "|quotient_order=" << item.order
            << "|regular_points=" << item.regular_points
            << "|regular_orbits=" << item.x_reps.size()
            << "|target_orbits=" << item.y_reps.size()
            << "|maximum_left_rep=" << maximum_left_index + 1
            << "|box_probes=" << probes
            << "|expanded_left_tests=" << expanded_left_tests
            << "|witness_digest=" << digest
            << "|status=PASS\n" << std::flush;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 3) {
      throw std::runtime_error("usage: engine FULL_TRANSCRIPT S3_TRANSCRIPT");
    }
    std::cout << "ENGINE|name=a5_k3_cd_correlated_m6_exact|schema=1\n";
    const CoreData core = build_core_data();
    const auto full = std::string(argv[1]) == "-"
        ? std::vector<Case6>{} : parse_representatives(argv[1], core);
    const auto s3 = std::string(argv[2]) == "-"
        ? std::vector<Case6>{} : parse_representatives(argv[2], core);
    if (full.empty() && s3.empty()) {
      throw std::runtime_error("unexpected residual-case census");
    }
    int thread_count = 1;
    if (const char* value = std::getenv("BG_M6_THREADS")) {
      const std::string text = value;
      const auto [end, error] = std::from_chars(
          text.data(), text.data() + text.size(), thread_count);
      if (error != std::errc{} || end != text.data() + text.size() ||
          thread_count < 1 || thread_count > 64) {
        throw std::runtime_error("BG_M6_THREADS must be an integer in 1..64");
      }
    }
    for (const auto& item : full) audit6("full", item, core, thread_count);
    for (const auto& item : s3) audit6("s3", item, core, thread_count);
    std::cout << "CORRELATED_M6_COMPLETE|cases=" << full.size() + s3.size()
              << "|failures=0|status=ALL_PASS\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ERROR|" << error.what() << '\n';
    return 1;
  }
}

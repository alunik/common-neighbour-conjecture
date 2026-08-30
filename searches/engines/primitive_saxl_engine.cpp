// primitive_saxl_engine.cpp
//
// Exact Saxl-graph engine for a finite transitive permutation action.
// Magma supplies the union R of the regular suborbits of a point stabiliser
// and, only when needed, generators for the ambient permutation group.  The
// engine transports R along a Schreier tree.  For stabiliser-orbit target
// representatives it tests regular points one at a time and stops at the
// first common-neighbour witness; a negative target still exhausts R.
//
// Build:
//   clang++ -O3 -march=native -DNDEBUG -std=c++17
//       src/primitive_saxl_engine.cpp -o primitive_saxl_engine
//
// Input (multiple action blocks may follow the header):
//   PRIMITIVE_SAXL_V1
//   action
//   label Primitive_55_1
//   degree 55
//   stabilizer_order 12
//   classification graph
//   regular_orbits 2
//   regular_count 24
//   regular 2 3 ...
//   gens 2
//   <each permutation either as degree one-based images, or as
//    "packed_gen width data" using base-94 little-endian image codes>
//   end
//
// classification is one of base1, order_obstruction, compute,
// no_regular_orbit, complete, density, graph.  Production exporters use
// compute and provide hgens plus gens; C++ then finds the regular suborbits.
// The derived forms are retained for compact regression fixtures.  JSONL
// results are written to stdout.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;

#ifndef PRIMITIVE_SAXL_DENSE_LIMIT_BYTES
#define PRIMITIVE_SAXL_DENSE_LIMIT_BYTES (256ULL * 1024ULL * 1024ULL)
#endif

namespace {

constexpr u32 UNSEEN = std::numeric_limits<u32>::max();

struct Action {
    std::string label;
    u32 degree = 0;
    std::string stabilizer_order;
    std::string classification;
    u32 regular_orbits = 0;
    u32 regular_count = 0;
    std::vector<u32> regular;
    std::vector<u32> orbit_representatives;
    std::vector<u32> orbit_sizes;
    std::vector<std::vector<u32>> stabilizer_generators;
    std::vector<std::vector<u32>> generators;
    bool stabilizer_orbit_enumeration_complete = true;
};

std::string jsonEscape(const std::string& text) {
    std::ostringstream out;
    for (unsigned char c : text) {
        switch (c) {
            case '\\': out << "\\\\"; break;
            case '"': out << "\\\""; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (c < 0x20) {
                    out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                        << static_cast<unsigned>(c) << std::dec;
                } else {
                    out << static_cast<char>(c);
                }
        }
    }
    return out.str();
}

bool decimalString(const std::string& value) {
    return !value.empty() && value != "0" &&
           std::all_of(value.begin(), value.end(), [](unsigned char c) {
               return c >= '0' && c <= '9';
           });
}

u32 oneBasedPoint(const std::string& text, u32 degree,
                  const std::string& description) {
    if (text.empty() ||
        !std::all_of(text.begin(), text.end(), [](unsigned char c) {
            return c >= '0' && c <= '9';
        }))
        throw std::runtime_error("bad " + description);
    u64 value = 0;
    for (unsigned char c : text) {
        value = 10 * value + static_cast<unsigned>(c - '0');
        if (value > degree) throw std::runtime_error("bad " + description);
    }
    if (value < 1) throw std::runtime_error("bad " + description);
    return static_cast<u32>(value - 1);
}

void readGenerator(std::istream& input, std::vector<u32>& generator,
                   u32 degree, const std::string& description) {
    std::string first;
    if (!(input >> first)) throw std::runtime_error("missing " + description);
    if (first != "packed_gen") {
        generator[0] = oneBasedPoint(first, degree, description + " image");
        for (u32 index = 1; index < degree; ++index) {
            std::string image;
            if (!(input >> image)) throw std::runtime_error("missing " + description);
            generator[index] = oneBasedPoint(image, degree, description + " image");
        }
        return;
    }

    u32 width = 0;
    std::string packed;
    if (!(input >> width >> packed) || width < 1 || width > 5 ||
        packed.size() != static_cast<std::size_t>(degree) * width)
        throw std::runtime_error("bad packed " + description);
    for (u32 point = 0; point < degree; ++point) {
        u64 image = 0;
        u64 place = 1;
        for (u32 digit_index = 0; digit_index < width; ++digit_index) {
            const unsigned char code = static_cast<unsigned char>(
                packed[static_cast<std::size_t>(point) * width + digit_index]);
            if (code < 33 || code > 126)
                throw std::runtime_error("bad packed " + description + " digit");
            image += static_cast<u64>(code - 33) * place;
            place *= 94;
        }
        if (image >= degree)
            throw std::runtime_error("bad packed " + description + " image");
        generator[point] = static_cast<u32>(image);
    }
}

void readActions(std::istream& input, const std::function<void(Action&)>& consume) {
    std::string token;
    if (!(input >> token) || token != "PRIMITIVE_SAXL_V1")
        throw std::runtime_error("input must start with PRIMITIVE_SAXL_V1");

    while (input >> token) {
        if (token != "action")
            throw std::runtime_error("expected 'action', got '" + token + "'");
        Action action;
        bool have_label = false, have_degree = false, have_order = false;
        bool have_classification = false, have_regular_orbits = false;
        bool have_regular_count = false, have_regular = false;
        bool have_orbit_count = false, have_orbit_representatives = false;
        bool have_hgens = false, have_gens = false;
        bool ended = false;
        while (input >> token) {
            if (token == "end") {
                ended = true;
                break;
            }
            if (token == "label") {
                if (!(input >> action.label) || action.label.empty())
                    throw std::runtime_error("missing action label");
                have_label = true;
            } else if (token == "degree") {
                u64 value = 0;
                if (!(input >> value) || value < 2 || value > std::numeric_limits<u32>::max())
                    throw std::runtime_error("bad action degree");
                action.degree = static_cast<u32>(value);
                have_degree = true;
            } else if (token == "stabilizer_order") {
                if (!(input >> action.stabilizer_order) ||
                    !decimalString(action.stabilizer_order))
                    throw std::runtime_error("bad stabilizer order");
                have_order = true;
            } else if (token == "classification") {
                if (!(input >> action.classification))
                    throw std::runtime_error("missing classification");
                have_classification = true;
            } else if (token == "regular_orbits") {
                u64 value = 0;
                if (!(input >> value) || value > std::numeric_limits<u32>::max())
                    throw std::runtime_error("bad regular-orbit count");
                action.regular_orbits = static_cast<u32>(value);
                have_regular_orbits = true;
            } else if (token == "regular_count") {
                u64 value = 0;
                if (!(input >> value) || value > std::numeric_limits<u32>::max())
                    throw std::runtime_error("bad regular-point count");
                action.regular_count = static_cast<u32>(value);
                have_regular_count = true;
            } else if (token == "regular") {
                if (!have_degree || !have_regular_count)
                    throw std::runtime_error("degree and regular_count must precede regular");
                action.regular.resize(action.regular_count);
                for (u32& point : action.regular) {
                    u64 value = 0;
                    if (!(input >> value) || value < 1 || value > action.degree)
                        throw std::runtime_error("bad regular point");
                    point = static_cast<u32>(value - 1);
                }
                have_regular = true;
            } else if (token == "orbit_representatives_count") {
                u64 count = 0;
                if (!(input >> count) || count > action.degree)
                    throw std::runtime_error("bad orbit-representative count");
                action.orbit_representatives.resize(static_cast<std::size_t>(count));
                action.orbit_sizes.resize(static_cast<std::size_t>(count));
                have_orbit_count = true;
            } else if (token == "orbit_representatives") {
                if (!have_degree || !have_orbit_count)
                    throw std::runtime_error(
                        "degree and orbit_representatives_count must precede representatives");
                for (std::size_t index = 0;
                     index < action.orbit_representatives.size(); ++index) {
                    std::string point;
                    u64 size = 0;
                    if (!(input >> point >> size) || size < 1 || size > action.degree)
                        throw std::runtime_error("bad orbit representative");
                    action.orbit_representatives[index] =
                        oneBasedPoint(point, action.degree, "orbit representative");
                    action.orbit_sizes[index] = static_cast<u32>(size);
                }
                have_orbit_representatives = true;
            } else if (token == "hgens") {
                if (!have_degree) throw std::runtime_error("degree must precede hgens");
                u64 count = 0;
                if (!(input >> count) || count > std::numeric_limits<u32>::max())
                    throw std::runtime_error("bad stabilizer-generator count");
                action.stabilizer_generators.assign(
                    static_cast<std::size_t>(count), std::vector<u32>(action.degree));
                for (auto& generator : action.stabilizer_generators) {
                    readGenerator(input, generator, action.degree,
                                  "stabilizer permutation");
                }
                have_hgens = true;
            } else if (token == "gens") {
                if (!have_degree) throw std::runtime_error("degree must precede gens");
                u64 count = 0;
                if (!(input >> count) || count > std::numeric_limits<u32>::max())
                    throw std::runtime_error("bad generator count");
                action.generators.assign(static_cast<std::size_t>(count),
                                         std::vector<u32>(action.degree));
                for (auto& generator : action.generators) {
                    readGenerator(input, generator, action.degree, "permutation");
                }
                have_gens = true;
            } else {
                throw std::runtime_error("unknown input field '" + token + "'");
            }
        }
        if (!ended) throw std::runtime_error("unterminated action block");
        if (!(have_label && have_degree && have_order && have_classification &&
              have_regular_orbits && have_regular_count))
            throw std::runtime_error("incomplete action block");

        const bool graph = action.classification == "graph";
        const bool compute = action.classification == "compute";
        if (have_orbit_count != have_orbit_representatives)
            throw std::runtime_error("incomplete orbit-representative data");
        if ((graph && !(have_regular && have_gens)) ||
            (compute && !(have_hgens && have_gens)))
            throw std::runtime_error("graph/compute block is missing permutation data");
        if (!graph && !compute && (have_regular || have_hgens || have_gens))
            throw std::runtime_error("terminal block contains unnecessary permutation data");
        static const std::vector<std::string> valid = {
            "base1", "order_obstruction", "no_regular_orbit",
            "complete", "density", "graph", "compute"
        };
        if (std::find(valid.begin(), valid.end(), action.classification) == valid.end())
            throw std::runtime_error("unknown classification '" + action.classification + "'");
        consume(action);
    }
}

class BitMatrix {
  public:
    explicit BitMatrix(u32 n)
        : n_(n), words_((static_cast<u64>(n) + 63) / 64),
          data_(static_cast<std::size_t>(n) * words_, 0) {}

    u64* row(u32 vertex) { return data_.data() + static_cast<std::size_t>(vertex) * words_; }
    const u64* row(u32 vertex) const {
        return data_.data() + static_cast<std::size_t>(vertex) * words_;
    }
    std::size_t words() const { return words_; }

    void set(u32 vertex, u32 neighbour) {
        row(vertex)[neighbour >> 6] |= u64{1} << (neighbour & 63);
    }
    bool test(u32 vertex, u32 neighbour) const {
        return ((row(vertex)[neighbour >> 6] >> (neighbour & 63)) & 1U) != 0;
    }
    bool intersects(u32 left, u32 right) const {
        const u64* a = row(left);
        const u64* b = row(right);
        for (std::size_t word = 0; word < words_; ++word)
            if ((a[word] & b[word]) != 0) return true;
        return false;
    }

    template <typename Function>
    void forEachSetBit(u32 vertex, Function function) const {
        const u64* bits = row(vertex);
        for (std::size_t word = 0; word < words_; ++word) {
            u64 value = bits[word];
            while (value != 0) {
                const unsigned offset = static_cast<unsigned>(__builtin_ctzll(value));
                const u32 point = static_cast<u32>(64 * word + offset);
                if (point < n_) function(point);
                value &= value - 1;
            }
        }
    }

  private:
    u32 n_;
    std::size_t words_;
    std::vector<u64> data_;
};

void validatePermutation(const std::vector<u32>& permutation, u32 degree) {
    if (permutation.size() != degree) throw std::runtime_error("wrong permutation degree");
    std::vector<unsigned char> seen(degree, 0);
    for (u32 image : permutation) {
        if (image >= degree || seen[image]) throw std::runtime_error("generator is not a permutation");
        seen[image] = 1;
    }
}

void verifyTransitiveGenerators(
    const std::vector<std::vector<u32>>& generators, u32 degree) {
    if (generators.empty()) throw std::runtime_error("ambient group has no generators");
    std::vector<unsigned char> reached(degree, 0);
    std::vector<u32> queue;
    queue.reserve(degree);
    reached[0] = 1;
    queue.push_back(0);
    for (std::size_t head = 0; head < queue.size(); ++head) {
        for (const auto& generator : generators) {
            const u32 image = generator[queue[head]];
            if (reached[image]) continue;
            reached[image] = 1;
            queue.push_back(image);
        }
    }
    if (queue.size() != degree)
        throw std::runtime_error("ambient generators are not transitive");
}

u32 smallDecimal(const std::string& text, u32 maximum) {
    u64 value = 0;
    for (char c : text) {
        value = 10 * value + static_cast<unsigned>(c - '0');
        if (value > maximum) throw std::runtime_error("stabilizer order exceeds compute bound");
    }
    return static_cast<u32>(value);
}

void computeRegularSuborbits(Action& action) {
    const u32 n = action.degree;
    const u32 order = smallDecimal(action.stabilizer_order, n - 1);
    if (order <= 1) throw std::runtime_error("compute block has trivial stabilizer");
    if (action.stabilizer_generators.empty())
        throw std::runtime_error("compute block has no stabilizer generators");
    for (const auto& generator : action.stabilizer_generators) {
        validatePermutation(generator, n);
        if (generator[0] != 0) throw std::runtime_error("hgen does not fix the base point");
    }
    for (const auto& generator : action.generators) validatePermutation(generator, n);

    std::vector<unsigned char> seen(n, 0);
    std::vector<u32> queue;
    for (u32 seed = 0; seed < n; ++seed) {
        if (seen[seed]) continue;
        queue.clear();
        seen[seed] = 1;
        queue.push_back(seed);
        for (std::size_t head = 0; head < queue.size(); ++head) {
            const u32 point = queue[head];
            for (const auto& generator : action.stabilizer_generators) {
                const u32 image = generator[point];
                if (!seen[image]) {
                    seen[image] = 1;
                    queue.push_back(image);
                }
            }
        }
        action.orbit_representatives.push_back(seed);
        action.orbit_sizes.push_back(static_cast<u32>(queue.size()));
        if (order % queue.size() != 0)
            throw std::runtime_error("stabilizer orbit length does not divide its order");
        if (queue.size() == order) {
            ++action.regular_orbits;
            action.regular.insert(action.regular.end(), queue.begin(), queue.end());
            if (action.regular.size() == n - 1) {
                action.regular_count = static_cast<u32>(action.regular.size());
                action.classification = "complete";
                verifyTransitiveGenerators(action.generators, n);
                return;
            }
            if (2 * static_cast<u64>(action.regular.size()) > n) {
                // The strict density inequality alone certifies a common
                // neighbour for every pair.  Retain the regular points and
                // orbit receipts already found, but do not enumerate the
                // remaining stabiliser orbits merely to obtain an unused
                // exact valency.
                action.regular_count = static_cast<u32>(action.regular.size());
                action.classification = "density";
                action.stabilizer_orbit_enumeration_complete = false;
                verifyTransitiveGenerators(action.generators, n);
                return;
            }
        }
    }
    action.regular_count = static_cast<u32>(action.regular.size());
    if (action.regular_count == 0) {
        action.classification = "no_regular_orbit";
        verifyTransitiveGenerators(action.generators, n);
    }
    else if (action.regular_count == n - 1) action.classification = "complete";
    else if (2 * static_cast<u64>(action.regular_count) > n) action.classification = "density";
    else action.classification = "graph";
}

struct GraphResult {
    bool common_neighbour_all_pairs = true;
    u32 common_neighbour_failures = 0;
    bool diameter_at_most_2 = true;
    u32 distance_ge3_witness = UNSEEN;
    bool diameter_exact = true;
    u32 diameter = 0;
    std::vector<u32> distance_layers;
};

GraphResult analyseGraphDense(const Action& action) {
    const u32 n = action.degree;
    const u32 k = action.regular_count;
    if (action.generators.empty()) throw std::runtime_error("graph block has no generators");
    if (k == 0 || 2 * static_cast<u64>(k) > n)
        throw std::runtime_error("graph block violates sparse-graph precondition");
    for (const auto& generator : action.generators) validatePermutation(generator, n);

    std::vector<unsigned char> in_regular(n, 0);
    for (u32 point : action.regular) {
        if (point == 0) throw std::runtime_error("base point occurs in regular suborbits");
        if (in_regular[point]) throw std::runtime_error("duplicate regular point");
        in_regular[point] = 1;
    }

    BitMatrix adjacency(n);
    for (u32 point : action.regular) adjacency.set(0, point);

    // A Schreier tree transports the one known adjacency row.  A set bit x
    // in row(v) becomes x^g in row(v^g).  Each row is constructed once.
    std::vector<u32> queue;
    queue.reserve(n);
    std::vector<unsigned char> reached(n, 0);
    reached[0] = 1;
    queue.push_back(0);
    for (std::size_t head = 0; head < queue.size(); ++head) {
        const u32 vertex = queue[head];
        for (const auto& generator : action.generators) {
            const u32 image_vertex = generator[vertex];
            if (reached[image_vertex]) continue;
            reached[image_vertex] = 1;
            adjacency.forEachSetBit(vertex, [&](u32 neighbour) {
                adjacency.set(image_vertex, generator[neighbour]);
            });
            queue.push_back(image_vertex);
        }
    }
    if (queue.size() != n) throw std::runtime_error("ambient generators are not transitive");

    // Symmetry is a compact end-to-end check on the exported regular set,
    // permutation orientation, and Schreier transport.
    for (u32 vertex = 0; vertex < n; ++vertex) {
        if (adjacency.test(vertex, vertex)) throw std::runtime_error("Saxl graph has a loop");
        u32 row_count = 0;
        adjacency.forEachSetBit(vertex, [&](u32 neighbour) {
            ++row_count;
            if (!adjacency.test(neighbour, vertex))
                throw std::runtime_error("transported adjacency is not symmetric");
        });
        if (row_count != k) throw std::runtime_error("transported row has wrong valency");
    }

    GraphResult result;
    for (u32 vertex = 1; vertex < n; ++vertex) {
        const bool common = adjacency.intersects(0, vertex);
        if (!common) {
            result.common_neighbour_all_pairs = false;
            ++result.common_neighbour_failures;
            if (!adjacency.test(0, vertex) && result.distance_ge3_witness == UNSEEN) {
                result.diameter_at_most_2 = false;
                result.distance_ge3_witness = vertex;
            }
        }
    }

    if (result.diameter_at_most_2) {
        result.diameter = k == n - 1 ? 1 : 2;
        result.distance_layers = {1, k, n - 1 - k};
        if (result.diameter == 1) result.distance_layers.resize(2);
        return result;
    }

    // Only exceptional actions pay for the full BFS.  Since each adjacency
    // row is already available, every edge is examined at most twice.
    std::vector<u32> distance(n, UNSEEN);
    queue.clear();
    distance[0] = 0;
    queue.push_back(0);
    result.distance_layers = {1};
    for (std::size_t head = 0; head < queue.size(); ++head) {
        const u32 vertex = queue[head];
        adjacency.forEachSetBit(vertex, [&](u32 neighbour) {
            if (distance[neighbour] != UNSEEN) return;
            distance[neighbour] = distance[vertex] + 1;
            if (result.distance_layers.size() <= distance[neighbour])
                result.distance_layers.resize(distance[neighbour] + 1, 0);
            ++result.distance_layers[distance[neighbour]];
            queue.push_back(neighbour);
        });
    }
    if (queue.size() != n) throw std::runtime_error("nonempty primitive orbital union is disconnected");
    result.diameter = static_cast<u32>(result.distance_layers.size() - 1);
    return result;
}

GraphResult analyseGraphStreaming(const Action& action) {
    const u32 n = action.degree;
    const u32 k = action.regular_count;
    if (action.generators.empty()) throw std::runtime_error("graph block has no generators");
    if (k == 0 || 2 * static_cast<u64>(k) > n)
        throw std::runtime_error("graph block violates sparse-graph precondition");
    for (const auto& generator : action.generators) validatePermutation(generator, n);

    std::vector<unsigned char> in_regular(n, 0);
    for (u32 point : action.regular) {
        if (point == 0) throw std::runtime_error("base point occurs in regular suborbits");
        if (in_regular[point]) throw std::runtime_error("duplicate regular point");
        in_regular[point] = 1;
    }

    // Build only the Schreier tree, not all n adjacency rows.  The tree is
    // then traversed depth first while a single copy of R is transported in
    // place.  Applying the inverse generator on return restores the parent
    // row.  This uses O(n + |R|) working memory instead of O(n^2) bits.
    const u32 generator_count = static_cast<u32>(action.generators.size());
    std::vector<std::vector<u32>> inverses(
        generator_count, std::vector<u32>(n));
    for (u32 index = 0; index < generator_count; ++index) {
        for (u32 point = 0; point < n; ++point)
            inverses[index][action.generators[index][point]] = point;
    }

    std::vector<u32> first_child(n, UNSEEN), next_sibling(n, UNSEEN);
    std::vector<u32> parent_generator(n, UNSEEN), queue;
    std::vector<unsigned char> reached(n, 0);
    queue.reserve(n);
    reached[0] = 1;
    queue.push_back(0);
    for (std::size_t head = 0; head < queue.size(); ++head) {
        const u32 vertex = queue[head];
        for (u32 index = 0; index < generator_count; ++index) {
            const u32 child = action.generators[index][vertex];
            if (reached[child]) continue;
            reached[child] = 1;
            parent_generator[child] = index;
            next_sibling[child] = first_child[vertex];
            first_child[vertex] = child;
            queue.push_back(child);
        }
    }
    if (queue.size() != n) throw std::runtime_error("ambient generators are not transitive");

    struct Frame { u32 vertex; u32 next_child; };
    std::vector<Frame> stack;
    stack.reserve(n);
    stack.push_back({0, first_child[0]});
    std::vector<u32> transported = action.regular;
    GraphResult result;

    auto inspect_row = [&](u32 vertex) {
        bool common = false;
        bool contains_base = false;
        bool contains_vertex = false;
        for (u32 point : transported) {
            common = common || in_regular[point];
            contains_base = contains_base || point == 0;
            contains_vertex = contains_vertex || point == vertex;
        }
        if (contains_vertex) throw std::runtime_error("Saxl graph has a loop");
        if (contains_base != static_cast<bool>(in_regular[vertex]))
            throw std::runtime_error("transported adjacency is not symmetric");
        if (!common && vertex != 0) {
            result.common_neighbour_all_pairs = false;
            ++result.common_neighbour_failures;
            if (!in_regular[vertex] && result.distance_ge3_witness == UNSEEN) {
                result.diameter_at_most_2 = false;
                result.distance_ge3_witness = vertex;
            }
        }
    };
    inspect_row(0);

    while (!stack.empty()) {
        Frame& frame = stack.back();
        if (frame.next_child != UNSEEN) {
            const u32 child = frame.next_child;
            frame.next_child = next_sibling[child];
            const u32 index = parent_generator[child];
            for (u32& point : transported) point = action.generators[index][point];
            inspect_row(child);
            stack.push_back({child, first_child[child]});
        } else {
            const u32 vertex = frame.vertex;
            stack.pop_back();
            if (vertex != 0) {
                const u32 index = parent_generator[vertex];
                for (u32& point : transported) point = inverses[index][point];
            }
        }
    }
    if (transported != action.regular)
        throw std::runtime_error("Schreier traversal failed to restore regular set");

    if (result.diameter_at_most_2) {
        result.diameter = k == n - 1 ? 1 : 2;
        result.distance_layers = {1, k, n - 1 - k};
        if (result.diameter == 1) result.distance_layers.resize(2);
    } else {
        result.diameter_exact = false;
    }
    return result;
}

GraphResult analyseGraphByOrbitRepresentatives(const Action& action) {
    const u32 n = action.degree;
    const u32 k = action.regular_count;
    if (action.generators.empty() || action.orbit_representatives.empty())
        throw std::runtime_error("representative graph block is incomplete");
    if (k == 0 || 2 * static_cast<u64>(k) > n)
        throw std::runtime_error("graph block violates sparse-graph precondition");
    for (const auto& generator : action.generators) validatePermutation(generator, n);

    std::vector<unsigned char> in_regular(n, 0), seen_representative(n, 0);
    for (u32 point : action.regular) {
        if (point == 0) throw std::runtime_error("base point occurs in regular suborbits");
        if (in_regular[point]) throw std::runtime_error("duplicate regular point");
        in_regular[point] = 1;
    }
    u64 orbit_size_sum = 0;
    bool have_base = false;
    for (std::size_t index = 0; index < action.orbit_representatives.size(); ++index) {
        const u32 representative = action.orbit_representatives[index];
        if (representative >= n || seen_representative[representative])
            throw std::runtime_error("bad/duplicate orbit representative");
        seen_representative[representative] = 1;
        have_base = have_base || representative == 0;
        orbit_size_sum += action.orbit_sizes[index];
    }
    if (!have_base || orbit_size_sum != n)
        throw std::runtime_error("orbit representatives do not partition the degree");

    // Add inverse generators so the Schreier tree is shallow, then record a
    // certified transporter word from the base point to every vertex.
    std::vector<std::vector<u32>> transports = action.generators;
    for (const auto& generator : action.generators) {
        std::vector<u32> inverse(n);
        for (u32 point = 0; point < n; ++point) inverse[generator[point]] = point;
        transports.push_back(std::move(inverse));
    }
    std::vector<u32> parent(n, UNSEEN), parent_transport(n, UNSEEN), queue;
    std::vector<unsigned char> reached(n, 0);
    queue.reserve(n);
    reached[0] = 1;
    queue.push_back(0);
    for (std::size_t head = 0; head < queue.size(); ++head) {
        const u32 vertex = queue[head];
        for (u32 index = 0; index < transports.size(); ++index) {
            const u32 child = transports[index][vertex];
            if (reached[child]) continue;
            reached[child] = 1;
            parent[child] = vertex;
            parent_transport[child] = index;
            queue.push_back(child);
        }
    }
    if (queue.size() != n) throw std::runtime_error("ambient generators are not transitive");

    GraphResult result;
    const u32 generator_count = static_cast<u32>(action.generators.size());
    auto inverse_transport = [generator_count](u32 index) {
        return index < generator_count ? index + generator_count
                                       : index - generator_count;
    };
    std::vector<u32> path;
    for (std::size_t orbit_index = 0;
         orbit_index < action.orbit_representatives.size(); ++orbit_index) {
        const u32 representative = action.orbit_representatives[orbit_index];
        path.clear();
        for (u32 vertex = representative; vertex != 0; vertex = parent[vertex]) {
            if (parent[vertex] == UNSEEN || parent_transport[vertex] == UNSEEN)
                throw std::runtime_error("broken Schreier transporter");
            path.push_back(parent_transport[vertex]);
        }

        // Replay the transporter on the base point and its inverse on the
        // base point.  The latter gives the unique point whose transported
        // row contains the base, providing the old symmetry check without
        // scanning the whole transported regular set.
        u32 image_of_base = 0;
        for (auto iterator = path.rbegin(); iterator != path.rend(); ++iterator) {
            image_of_base = transports[*iterator][image_of_base];
        }
        if (image_of_base != representative)
            throw std::runtime_error("Schreier transporter has wrong endpoint");
        u32 preimage_of_base = 0;
        for (u32 index : path)
            preimage_of_base = transports[inverse_transport(index)][preimage_of_base];
        if (static_cast<bool>(in_regular[preimage_of_base]) !=
            static_cast<bool>(in_regular[representative]))
            throw std::runtime_error("transported adjacency is not symmetric");

        bool common = false;
        for (u32 regular_point : action.regular) {
            u32 transported_point = regular_point;
            for (auto iterator = path.rbegin(); iterator != path.rend(); ++iterator)
                transported_point = transports[*iterator][transported_point];
            if (in_regular[transported_point]) {
                common = true;
                break;
            }
        }
        if (!common) {
            result.common_neighbour_all_pairs = false;
            result.common_neighbour_failures += action.orbit_sizes[orbit_index];
            if (!in_regular[representative] && representative != 0 &&
                result.distance_ge3_witness == UNSEEN) {
                result.diameter_at_most_2 = false;
                result.distance_ge3_witness = representative;
            }
        }
    }

    if (result.diameter_at_most_2) {
        result.diameter = k == n - 1 ? 1 : 2;
        result.distance_layers = {1, k, n - 1 - k};
        if (result.diameter == 1) result.distance_layers.resize(2);
    } else {
        result.diameter_exact = false;
    }
    return result;
}

GraphResult analyseGraph(const Action& action) {
    // Dense storage is retained for modest exceptional graphs because it
    // gives the exact diameter by BFS.  Larger graphs use the linear-memory
    // screen; any distance-at-least-three witness is then rerun separately.
    const u64 words = (static_cast<u64>(action.degree) + 63) / 64;
    const u64 dense_bytes = static_cast<u64>(action.degree) * words * 8;
    constexpr u64 DENSE_LIMIT_BYTES = PRIMITIVE_SAXL_DENSE_LIMIT_BYTES;
    if (dense_bytes <= DENSE_LIMIT_BYTES) return analyseGraphDense(action);
    if (!action.orbit_representatives.empty())
        return analyseGraphByOrbitRepresentatives(action);
    return analyseGraphStreaming(action);
}

void printLayers(const std::vector<u32>& layers) {
    std::cout << '[';
    for (std::size_t i = 0; i < layers.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << layers[i];
    }
    std::cout << ']';
}

void emitResult(Action& action) {
    const auto started = std::chrono::steady_clock::now();
    if (action.classification == "compute") computeRegularSuborbits(action);
    const u32 n = action.degree;
    const u32 k = action.regular_count;

    std::string status;
    std::string base_size;
    bool common_neighbour = false;
    u32 common_failures = 0;
    bool diameter_at_most_2 = false;
    bool has_diameter = false;
    u32 diameter = 0;
    u32 witness = UNSEEN;
    std::vector<u32> layers;

    if (action.classification == "base1") {
        if (k != 0) throw std::runtime_error("base1 block has regular points");
        status = "trivial_stabiliser";
        base_size = "1";
    } else if (action.classification == "order_obstruction") {
        if (k != 0) throw std::runtime_error("order-obstruction block has regular points");
        status = "order_obstruction";
        base_size = ">2";
    } else if (action.classification == "no_regular_orbit") {
        if (k != 0) throw std::runtime_error("no-regular-orbit block has regular points");
        status = "no_regular_orbit";
        base_size = ">2";
    } else if (action.classification == "complete") {
        if (k != n - 1 || action.regular_orbits == 0)
            throw std::runtime_error("bad complete-graph certificate");
        status = "complete_certificate";
        base_size = "2";
        common_neighbour = n >= 3;
        common_failures = n == 2 ? 1 : 0;
        diameter_at_most_2 = true;
        has_diameter = true;
        diameter = 1;
        layers = {1, n - 1};
    } else if (action.classification == "density") {
        if (2 * static_cast<u64>(k) <= n || k >= n || action.regular_orbits == 0)
            throw std::runtime_error("bad density certificate");
        status = "density_certificate";
        base_size = "2";
        common_neighbour = true;
        diameter_at_most_2 = true;
        // In compute mode k can be only the regular mass discovered before
        // the strict-density shortcut fired.  This proves the common-
        // neighbour property, but it does not distinguish a complete Saxl
        // graph (diameter 1) from a noncomplete graph of diameter 2 and it
        // does not give the exact distance-layer sizes.
        has_diameter = false;
    } else {
        if (action.regular_orbits == 0 || k == 0)
            throw std::runtime_error("graph block has no regular suborbit");
        const GraphResult graph = analyseGraph(action);
        status = graph.diameter_at_most_2 ? "graph_certificate" : "diameter_at_least_3";
        base_size = "2";
        common_neighbour = graph.common_neighbour_all_pairs;
        common_failures = graph.common_neighbour_failures;
        diameter_at_most_2 = graph.diameter_at_most_2;
        has_diameter = graph.diameter_exact;
        diameter = graph.diameter;
        witness = graph.distance_ge3_witness;
        layers = graph.distance_layers;
    }

    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout << "{\"schema\":\"PRIMITIVE_SAXL_RESULT_V1\""
              << ",\"label\":\"" << jsonEscape(action.label) << "\""
              << ",\"degree\":" << n
              << ",\"stabilizer_order_decimal\":\"" << action.stabilizer_order << "\""
              << ",\"classification\":\"" << action.classification << "\""
              << ",\"status\":\"" << status << "\""
              << ",\"base_size\":\"" << base_size << "\""
              << ",\"regular_orbits\":" << action.regular_orbits
              << ",\"regular_points\":" << k
              << ",\"stabilizer_orbit_representatives\":"
              << action.orbit_representatives.size()
              << ",\"stabilizer_orbit_enumeration_complete\":"
              << (action.stabilizer_orbit_enumeration_complete ? "true" : "false")
              << ",\"stabilizer_generators\":" << action.stabilizer_generators.size()
              << ",\"generators\":" << action.generators.size()
              << ",\"common_neighbour_all_pairs\":";
    if (base_size == "2") std::cout << (common_neighbour ? "true" : "false");
    else std::cout << "null";
    std::cout << ",\"common_neighbour_failures_from_basepoint\":";
    if (base_size == "2") std::cout << common_failures;
    else std::cout << "null";
    std::cout << ",\"diameter_at_most_2\":";
    if (base_size == "2") std::cout << (diameter_at_most_2 ? "true" : "false");
    else std::cout << "null";
    std::cout << ",\"diameter\":";
    if (has_diameter) std::cout << diameter;
    else std::cout << "null";
    std::cout << ",\"diameter_lower_bound\":";
    if (base_size == "2") {
        if (has_diameter) std::cout << diameter;
        else if (diameter_at_most_2) std::cout << 1;
        else std::cout << 3;
    } else std::cout << "null";
    std::cout << ",\"distance_layers\":";
    if (has_diameter) printLayers(layers);
    else std::cout << "null";
    std::cout << ",\"distance_ge3_witness\":";
    if (witness == UNSEEN) std::cout << "null";
    else std::cout << witness + 1;
    std::cout << ",\"seconds\":" << std::fixed << std::setprecision(6) << seconds
              << "}\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc > 2) throw std::runtime_error("usage: primitive_saxl_engine [input|-]");
        if (argc == 2 && std::string(argv[1]) != "-") {
            std::ifstream input(argv[1]);
            if (!input) throw std::runtime_error("cannot open input file");
            readActions(input, [](Action& action) { emitResult(action); });
        } else {
            readActions(std::cin, [](Action& action) { emitResult(action); });
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "primitive_saxl_engine: " << error.what() << '\n';
        return 1;
    }
}

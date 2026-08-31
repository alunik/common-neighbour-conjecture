// affine_saxl_engine.cpp
//
// Exact batch engine for affine groups G = V : H, where V = F_p^n and H is
// supplied by prime-field matrix generators together with its exact order.
//
// Mathematics used by the engine:
//   * v is a regular vector iff |v^H| = |H|;
//   * the Saxl graph is Cay(V,R), where R is the union of regular H-orbits;
//   * R+X is H-invariant when X is H-invariant, so a sumset can be computed
//     from one representative of one operand and every vector of the other;
//   * graph distance is computed by an exact orbit-compressed BFS.
//
// Build from the repository root:
//   clang++ -O3 -march=native -DNDEBUG -std=c++17 \
//       searches/engines/affine_saxl_engine.cpp -o affine_saxl_engine
//
// Input (several ACTION blocks may follow one header):
//   AFFINE_SAXL_V1
//   action
//   label <whitespace-free-label>
//   p <prime>
//   n <dimension>
//   order <exact-order-of-H>
//   orientation row|column
//   gens <number-of-generators>
//   <generator matrices, row-major, n*n integers each>
//   end
//
// Output is one JSON object per action on stdout. Diagnostics go to stderr.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <unordered_map>
#include <vector>

using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

namespace {

constexpr u32 UNASSIGNED = std::numeric_limits<u32>::max();

struct Options {
    u64 max_vectors = 1ULL << 24;
    u64 add_table_mib = 256;
    u64 pair_budget = 500000000ULL;
    int max_diameter = 12;
    bool diameter_ge3_only = false;
    bool orbit_fingerprint_only = false;
    bool dual_orbit_fingerprint = false;
    bool two_space_fingerprint = false;
};

struct Action {
    std::string label;
    u32 p = 0;
    u32 n = 0;
    u64 order = 0;
    std::string order_decimal;
    bool order_fits_u64 = false;
    bool row_action = true;
    std::vector<std::vector<u32>> gens;
};

std::string jsonEscape(const std::string& s) {
    std::ostringstream out;
    for (unsigned char c : s) {
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

u64 parseU64(const std::string& text, const char* option) {
    std::size_t used = 0;
    u64 value = 0;
    try {
        value = std::stoull(text, &used);
    } catch (...) {
        throw std::runtime_error(std::string("bad value for ") + option + ": " + text);
    }
    if (used != text.size())
        throw std::runtime_error(std::string("bad value for ") + option + ": " + text);
    return value;
}

bool isPrime(u32 p) {
    if (p < 2) return false;
    if ((p & 1U) == 0) return p == 2;
    for (u32 d = 3; static_cast<u64>(d) * d <= p; d += 2)
        if (p % d == 0) return false;
    return true;
}

u64 mulMod(u64 a, u64 b, u64 p) {
    return static_cast<u64>((static_cast<__uint128_t>(a) * b) % p);
}

u32 powMod(u32 a, u32 e, u32 p) {
    u64 result = 1, base = a;
    while (e) {
        if (e & 1U) result = mulMod(result, base, p);
        base = mulMod(base, base, p);
        e >>= 1U;
    }
    return static_cast<u32>(result);
}

bool invertible(const std::vector<u32>& matrix, u32 n, u32 p) {
    std::vector<u32> a = matrix;
    u32 rank = 0;
    for (u32 col = 0; col < n && rank < n; ++col) {
        u32 pivot = rank;
        while (pivot < n && a[static_cast<std::size_t>(pivot) * n + col] == 0) ++pivot;
        if (pivot == n) continue;
        if (pivot != rank) {
            for (u32 j = 0; j < n; ++j)
                std::swap(a[static_cast<std::size_t>(pivot) * n + j],
                          a[static_cast<std::size_t>(rank) * n + j]);
        }
        const u32 inv = powMod(a[static_cast<std::size_t>(rank) * n + col], p - 2, p);
        for (u32 j = col; j < n; ++j)
            a[static_cast<std::size_t>(rank) * n + j] =
                static_cast<u32>(mulMod(a[static_cast<std::size_t>(rank) * n + j], inv, p));
        for (u32 i = 0; i < n; ++i) {
            if (i == rank) continue;
            const u32 factor = a[static_cast<std::size_t>(i) * n + col];
            if (factor == 0) continue;
            for (u32 j = col; j < n; ++j) {
                const u64 sub = mulMod(factor, a[static_cast<std::size_t>(rank) * n + j], p);
                const u64 cur = a[static_cast<std::size_t>(i) * n + j];
                a[static_cast<std::size_t>(i) * n + j] =
                    static_cast<u32>((cur + p - sub) % p);
            }
        }
        ++rank;
    }
    return rank == n;
}

// Return A^{-T}.  For a row-vector action this is the dual representation.
// GL-conjugate natural modules have GL-conjugate duals, so the dual orbit-size
// multiset is another exact negative invariant.
std::vector<u32> inverseTranspose(const std::vector<u32>& matrix, u32 n, u32 p) {
    std::vector<u32> augmented(static_cast<std::size_t>(n) * 2 * n, 0);
    const u32 width = 2 * n;
    for (u32 row = 0; row < n; ++row) {
        for (u32 column = 0; column < n; ++column)
            augmented[static_cast<std::size_t>(row) * width + column] =
                matrix[static_cast<std::size_t>(row) * n + column];
        augmented[static_cast<std::size_t>(row) * width + n + row] = 1;
    }
    for (u32 column = 0; column < n; ++column) {
        u32 pivot = column;
        while (pivot < n &&
               augmented[static_cast<std::size_t>(pivot) * width + column] == 0)
            ++pivot;
        if (pivot == n) throw std::runtime_error("singular generator in dual action");
        if (pivot != column) {
            for (u32 j = 0; j < width; ++j)
                std::swap(augmented[static_cast<std::size_t>(pivot) * width + j],
                          augmented[static_cast<std::size_t>(column) * width + j]);
        }
        const u32 inverse = powMod(
            augmented[static_cast<std::size_t>(column) * width + column], p - 2, p);
        for (u32 j = 0; j < width; ++j)
            augmented[static_cast<std::size_t>(column) * width + j] =
                static_cast<u32>(mulMod(
                    augmented[static_cast<std::size_t>(column) * width + j], inverse, p));
        for (u32 row = 0; row < n; ++row) {
            if (row == column) continue;
            const u32 factor = augmented[static_cast<std::size_t>(row) * width + column];
            if (factor == 0) continue;
            for (u32 j = 0; j < width; ++j) {
                const u32 subtract = static_cast<u32>(mulMod(
                    factor, augmented[static_cast<std::size_t>(column) * width + j], p));
                u32& entry = augmented[static_cast<std::size_t>(row) * width + j];
                entry = (entry + p - subtract) % p;
            }
        }
    }
    std::vector<u32> dual(static_cast<std::size_t>(n) * n);
    for (u32 row = 0; row < n; ++row)
        for (u32 column = 0; column < n; ++column)
            dual[static_cast<std::size_t>(row) * n + column] =
                augmented[static_cast<std::size_t>(column) * width + n + row];
    return dual;
}

std::vector<Action> readActions(std::istream& in) {
    std::string token;
    if (!(in >> token) || token != "AFFINE_SAXL_V1")
        throw std::runtime_error("input must start with AFFINE_SAXL_V1");

    std::vector<Action> actions;
    while (in >> token) {
        if (token != "action") throw std::runtime_error("expected 'action', got '" + token + "'");
        Action a;
        bool have_label = false, have_p = false, have_n = false;
        bool have_order = false, have_orientation = false, have_gens = false;
        bool ended = false;
        while (in >> token) {
            if (token == "end") {
                ended = true;
                break;
            }
            if (token == "label") {
                if (!(in >> a.label)) throw std::runtime_error("missing label value");
                have_label = true;
            } else if (token == "p") {
                u64 x;
                if (!(in >> x) || x > std::numeric_limits<u32>::max())
                    throw std::runtime_error("bad field characteristic");
                a.p = static_cast<u32>(x);
                have_p = true;
            } else if (token == "n") {
                u64 x;
                if (!(in >> x) || x == 0 || x > std::numeric_limits<u32>::max())
                    throw std::runtime_error("bad dimension");
                a.n = static_cast<u32>(x);
                have_n = true;
            } else if (token == "order") {
                std::string decimal;
                if (!(in >> decimal) || decimal.empty() ||
                    !std::all_of(decimal.begin(), decimal.end(), [](unsigned char c) {
                        return c >= '0' && c <= '9';
                    }))
                    throw std::runtime_error("bad group order");
                const std::size_t first_nonzero = decimal.find_first_not_of('0');
                if (first_nonzero == std::string::npos) throw std::runtime_error("bad group order");
                a.order_decimal = decimal.substr(first_nonzero);
                a.order = 0;
                a.order_fits_u64 = true;
                for (char c : a.order_decimal) {
                    const u32 digit = static_cast<u32>(c - '0');
                    if (a.order > (std::numeric_limits<u64>::max() - digit) / 10) {
                        a.order_fits_u64 = false;
                        break;
                    }
                    a.order = 10 * a.order + digit;
                }
                have_order = true;
            } else if (token == "orientation") {
                std::string orientation;
                if (!(in >> orientation)) throw std::runtime_error("missing orientation");
                if (orientation == "row") a.row_action = true;
                else if (orientation == "column") a.row_action = false;
                else throw std::runtime_error("orientation must be row or column");
                have_orientation = true;
            } else if (token == "gens") {
                if (!have_p || !have_n)
                    throw std::runtime_error("p and n must precede gens");
                u64 count;
                if (!(in >> count) || count > std::numeric_limits<u32>::max())
                    throw std::runtime_error("bad generator count");
                a.gens.assign(static_cast<std::size_t>(count),
                              std::vector<u32>(static_cast<std::size_t>(a.n) * a.n));
                for (auto& matrix : a.gens) {
                    for (u32& entry : matrix) {
                        std::int64_t x;
                        if (!(in >> x)) throw std::runtime_error("truncated generator matrix");
                        const std::int64_t reduced =
                            ((x % static_cast<std::int64_t>(a.p)) + a.p) % a.p;
                        entry = static_cast<u32>(reduced);
                    }
                }
                have_gens = true;
            } else {
                throw std::runtime_error("unknown input field '" + token + "'");
            }
        }
        if (!ended) throw std::runtime_error("unterminated action block");
        if (!(have_label && have_p && have_n && have_order && have_orientation && have_gens))
            throw std::runtime_error("incomplete action block");
        actions.push_back(std::move(a));
    }
    return actions;
}

struct HalfAddTable {
    u64 q = 0;
    std::vector<u16> small;
    std::vector<u32> large;

    bool enabled() const { return !small.empty() || !large.empty(); }

    u32 get(u32 a, u32 b) const {
        const std::size_t pos = static_cast<std::size_t>(a) * static_cast<std::size_t>(q) + b;
        return small.empty() ? large[pos] : small[pos];
    }

    void build(u64 size, u32 p) {
        q = size;
        const u64 cells = q * q;
        if (q <= std::numeric_limits<u16>::max()) small.assign(static_cast<std::size_t>(cells), 0);
        else large.assign(static_cast<std::size_t>(cells), 0);
        for (u64 a = 0; a < q; ++a) {
            for (u64 b = 0; b < q; ++b) {
                if (a == 0 && b == 0) continue;
                const u64 parent = (a / p) * q + (b / p);
                const u32 parent_value = small.empty() ? large[static_cast<std::size_t>(parent)]
                                                       : small[static_cast<std::size_t>(parent)];
                const u32 value = static_cast<u32>(((a % p) + (b % p)) % p +
                                                   static_cast<u64>(p) * parent_value);
                const std::size_t pos = static_cast<std::size_t>(a * q + b);
                if (small.empty()) large[pos] = value;
                else small[pos] = static_cast<u16>(value);
            }
        }
    }
};

struct PackedVectorOps {
    u32 p;
    u32 n;
    std::vector<u64> powers;
    u64 vectors;
    u32 low_dimensions;
    u64 low_size;
    u64 high_size;
    bool xor_addition;
    bool table_addition = false;
    HalfAddTable low_add;
    HalfAddTable high_add;

    PackedVectorOps(u32 characteristic, u32 dimension, u64 table_mib)
        : p(characteristic), n(dimension), powers(dimension + 1, 1), vectors(1),
          low_dimensions((dimension + 1) / 2), low_size(1), high_size(1),
          xor_addition(characteristic == 2) {
        for (u32 i = 1; i <= n; ++i) powers[i] = powers[i - 1] * p;
        vectors = powers[n];
        low_size = powers[low_dimensions];
        high_size = powers[n - low_dimensions];
        if (xor_addition) return;

        const __uint128_t low_cells = static_cast<__uint128_t>(low_size) * low_size;
        const __uint128_t high_cells = static_cast<__uint128_t>(high_size) * high_size;
        const u64 low_bytes = low_size <= std::numeric_limits<u16>::max() ? 2 : 4;
        const u64 high_bytes = high_size <= std::numeric_limits<u16>::max() ? 2 : 4;
        const __uint128_t required = low_cells * low_bytes + high_cells * high_bytes;
        const __uint128_t allowed = static_cast<__uint128_t>(table_mib) << 20;
        if (required <= allowed &&
            low_cells <= std::numeric_limits<std::size_t>::max() &&
            high_cells <= std::numeric_limits<std::size_t>::max()) {
            low_add.build(low_size, p);
            high_add.build(high_size, p);
            table_addition = true;
        }
    }

    u32 add(u32 a, u32 b) const {
        if (a == 0) return b;
        if (b == 0) return a;
        if (xor_addition) return a ^ b;
        if (table_addition) {
            const u32 low = low_add.get(static_cast<u32>(a % low_size),
                                       static_cast<u32>(b % low_size));
            const u32 high = high_add.get(static_cast<u32>(a / low_size),
                                         static_cast<u32>(b / low_size));
            return static_cast<u32>(low + low_size * static_cast<u64>(high));
        }
        u64 x = a, y = b, result = 0, place = 1;
        for (u32 i = 0; i < n; ++i) {
            result += ((x % p + y % p) % p) * place;
            x /= p;
            y /= p;
            place *= p;
        }
        return static_cast<u32>(result);
    }

    std::string digits(u32 value) const {
        std::ostringstream out;
        out << '[';
        for (u32 i = 0; i < n; ++i) {
            if (i) out << ',';
            out << value % p;
            value /= p;
        }
        out << ']';
        return out.str();
    }
};

struct GeneratorTable {
    std::vector<u32> low;
    std::vector<u32> high;
};

u32 applyGeneratorToPacked(u32 value, const PackedVectorOps& ops,
                           const GeneratorTable& table) {
    const u32 low = static_cast<u32>(value % ops.low_size);
    const u32 high = static_cast<u32>(value / ops.low_size);
    return ops.add(table.low[low], table.high[high]);
}

u32 applyMatrixToPacked(u32 value, bool high_half, const Action& a,
                        const PackedVectorOps& ops, const std::vector<u32>& matrix) {
    std::vector<u32> input(a.n, 0);
    const u32 offset = high_half ? ops.low_dimensions : 0;
    const u32 count = high_half ? a.n - ops.low_dimensions : ops.low_dimensions;
    for (u32 i = 0; i < count; ++i) {
        input[offset + i] = value % a.p;
        value /= a.p;
    }

    u64 packed = 0;
    for (u32 out = 0; out < a.n; ++out) {
        __uint128_t sum = 0;
        if (a.row_action) {
            for (u32 in = 0; in < a.n; ++in)
                sum += static_cast<u64>(input[in]) * matrix[static_cast<std::size_t>(in) * a.n + out];
        } else {
            for (u32 in = 0; in < a.n; ++in)
                sum += static_cast<u64>(matrix[static_cast<std::size_t>(out) * a.n + in]) * input[in];
        }
        packed += static_cast<u64>(sum % a.p) * ops.powers[out];
    }
    return static_cast<u32>(packed);
}

std::vector<GeneratorTable> makeGeneratorTables(const Action& a, const PackedVectorOps& ops) {
    std::vector<GeneratorTable> tables(a.gens.size());
    for (std::size_t g = 0; g < a.gens.size(); ++g) {
        tables[g].low.resize(static_cast<std::size_t>(ops.low_size));
        tables[g].high.resize(static_cast<std::size_t>(ops.high_size));
        for (u64 v = 0; v < ops.low_size; ++v)
            tables[g].low[static_cast<std::size_t>(v)] =
                applyMatrixToPacked(static_cast<u32>(v), false, a, ops, a.gens[g]);
        for (u64 v = 0; v < ops.high_size; ++v)
            tables[g].high[static_cast<std::size_t>(v)] =
                applyMatrixToPacked(static_cast<u32>(v), true, a, ops, a.gens[g]);
    }
    return tables;
}

struct OrbitData {
    std::vector<u32> id;
    std::vector<u32> representative;
    std::vector<u32> size;
    bool complete = true;
    u64 regular_orbits_seen = 0;
    u64 regular_vectors_seen = 0;
};

OrbitData enumerateOrbits(const Action& a, const PackedVectorOps& ops,
                          const std::vector<GeneratorTable>& tables,
                          bool stop_after_density_certificate,
                          bool validate_orbit_orders = true) {
    OrbitData result;
    result.id.assign(static_cast<std::size_t>(ops.vectors), UNASSIGNED);
    std::vector<u32> queue;
    for (u64 seed64 = 0; seed64 < ops.vectors; ++seed64) {
        const u32 seed = static_cast<u32>(seed64);
        if (result.id[seed] != UNASSIGNED) continue;
        if (result.representative.size() == UNASSIGNED)
            throw std::runtime_error("too many H-orbits for 32-bit orbit identifiers");
        const u32 orbit = static_cast<u32>(result.representative.size());
        result.representative.push_back(seed);
        result.size.push_back(0);
        result.id[seed] = orbit;
        queue.clear();
        queue.push_back(seed);
        for (std::size_t head = 0; head < queue.size(); ++head) {
            const u32 v = queue[head];
            ++result.size[orbit];
            const u32 low = static_cast<u32>(v % ops.low_size);
            const u32 high = static_cast<u32>(v / ops.low_size);
            for (const auto& table : tables) {
                const u32 image = ops.add(table.low[low], table.high[high]);
                if (result.id[image] == UNASSIGNED) {
                    result.id[image] = orbit;
                    queue.push_back(image);
                }
            }
        }
        if (validate_orbit_orders &&
            (result.size[orbit] > a.order || a.order % result.size[orbit] != 0)) {
            std::ostringstream message;
            message << "orbit size " << result.size[orbit]
                    << " does not divide stated group order " << a.order;
            throw std::runtime_error(message.str());
        }
        if (result.size[orbit] == a.order) {
            ++result.regular_orbits_seen;
            result.regular_vectors_seen += result.size[orbit];
            if (stop_after_density_certificate &&
                2 * result.regular_vectors_seen > ops.vectors) {
                result.complete = false;
                return result;
            }
        }
    }
    return result;
}

// The multiset of orbit lengths on V is invariant under GL-conjugacy.  This
// inexpensive digest is used only to partition exact Magma conjugacy tests:
// a hash collision can cause an extra test, but can never fuse two groups.
std::pair<u64, u64> orbitSizeFingerprint(const OrbitData& orbits) {
    std::vector<u32> sizes = orbits.size;
    std::sort(sizes.begin(), sizes.end());
    u64 first = 1469598103934665603ULL;
    u64 second = 0x9e3779b97f4a7c15ULL;
    for (u32 size : sizes) {
        first ^= static_cast<u64>(size);
        first *= 1099511628211ULL;
        second ^= static_cast<u64>(size) + 0x9e3779b97f4a7c15ULL +
                  (second << 6U) + (second >> 2U);
    }
    first ^= static_cast<u64>(sizes.size());
    first *= 1099511628211ULL;
    second ^= static_cast<u64>(sizes.size()) * 0xbf58476d1ce4e5b9ULL;
    return {first, second};
}

struct ProjectiveOrbitFingerprint {
    u32 orbits = 0;
    u64 first = 0;
    u64 second = 0;
};

// Recover the H-orbits on projective points from the already enumerated
// vector orbits.  Scalar multiplication commutes with H, so the projective
// orbit containing v is exactly the union of the vector orbits containing
// c*v for c in F_p^*.  This avoids a second orbit enumeration and is much
// cheaper than Magma's generic OrbitsOfSpaces prefilter.
ProjectiveOrbitFingerprint projectiveOrbitSizeFingerprint(
        const OrbitData& vector_orbits, const PackedVectorOps& ops) {
    const std::size_t count = vector_orbits.size.size();
    std::vector<u32> parent(count);
    std::iota(parent.begin(), parent.end(), 0U);
    auto find = [&parent](u32 value) {
        u32 root = value;
        while (parent[root] != root) root = parent[root];
        while (parent[value] != value) {
            const u32 next = parent[value];
            parent[value] = root;
            value = next;
        }
        return root;
    };
    auto unite = [&parent, &find](u32 left, u32 right) {
        left = find(left);
        right = find(right);
        if (left != right) parent[right] = left;
    };
    auto scale = [&ops](u32 value, u32 scalar) {
        u64 result = 0;
        u64 place = 1;
        for (u32 coordinate = 0; coordinate < ops.n; ++coordinate) {
            result += static_cast<u64>((value % ops.p) * scalar % ops.p) * place;
            value /= ops.p;
            place *= ops.p;
        }
        return static_cast<u32>(result);
    };

    for (u32 orbit = 0; orbit < count; ++orbit) {
        const u32 representative = vector_orbits.representative[orbit];
        if (representative == 0) continue;
        for (u32 scalar = 2; scalar < ops.p; ++scalar) {
            unite(orbit, vector_orbits.id[scale(representative, scalar)]);
        }
    }

    std::vector<u64> root_vector_sizes(count, 0);
    for (u32 orbit = 0; orbit < count; ++orbit) {
        if (vector_orbits.representative[orbit] == 0) continue;
        root_vector_sizes[find(orbit)] += vector_orbits.size[orbit];
    }
    std::vector<u32> line_sizes;
    for (u64 vectors : root_vector_sizes) {
        if (vectors == 0) continue;
        if (vectors % (ops.p - 1) != 0)
            throw std::runtime_error("projective orbit size is not integral");
        const u64 lines = vectors / (ops.p - 1);
        if (lines > std::numeric_limits<u32>::max())
            throw std::runtime_error("projective orbit size exceeds 32 bits");
        line_sizes.push_back(static_cast<u32>(lines));
    }
    std::sort(line_sizes.begin(), line_sizes.end());
    u64 first = 1469598103934665603ULL;
    u64 second = 0x9e3779b97f4a7c15ULL;
    for (u32 size : line_sizes) {
        first ^= static_cast<u64>(size);
        first *= 1099511628211ULL;
        second ^= static_cast<u64>(size) + 0x9e3779b97f4a7c15ULL +
                  (second << 6U) + (second >> 2U);
    }
    first ^= static_cast<u64>(line_sizes.size());
    first *= 1099511628211ULL;
    second ^= static_cast<u64>(line_sizes.size()) * 0xbf58476d1ce4e5b9ULL;
    return {static_cast<u32>(line_sizes.size()), first, second};
}

struct TwoSpaceDomain {
    u32 p;
    std::vector<std::pair<u32, u32>> bases;
    std::unordered_map<u64, u32> index;

    explicit TwoSpaceDomain(u32 characteristic) : p(characteristic) {
        for (u32 first_pivot = 0; first_pivot < 4; ++first_pivot) {
            for (u32 second_pivot = first_pivot + 1; second_pivot < 4;
                 ++second_pivot) {
                std::vector<std::pair<u32, u32>> free_positions;
                for (u32 column = first_pivot + 1; column < 4; ++column)
                    if (column != second_pivot)
                        free_positions.emplace_back(0, column);
                for (u32 column = second_pivot + 1; column < 4; ++column)
                    free_positions.emplace_back(1, column);
                u64 assignments = 1;
                for (std::size_t i = 0; i < free_positions.size(); ++i)
                    assignments *= p;
                for (u64 code = 0; code < assignments; ++code) {
                    u32 rows[2][4] = {};
                    rows[0][first_pivot] = 1;
                    rows[1][second_pivot] = 1;
                    u64 value = code;
                    for (const auto& position : free_positions) {
                        rows[position.first][position.second] =
                            static_cast<u32>(value % p);
                        value /= p;
                    }
                    u32 packed[2] = {0, 0};
                    u64 key = 0, place = 1;
                    for (u32 row = 0; row < 2; ++row) {
                        u64 vector_place = 1;
                        for (u32 column = 0; column < 4; ++column) {
                            packed[row] += static_cast<u32>(rows[row][column] * vector_place);
                            vector_place *= p;
                            key += static_cast<u64>(rows[row][column]) * place;
                            place *= p;
                        }
                    }
                    const u32 position = static_cast<u32>(bases.size());
                    bases.emplace_back(packed[0], packed[1]);
                    if (!index.emplace(key, position).second)
                        throw std::runtime_error("duplicate two-space RREF key");
                }
            }
        }
        const u64 q2 = static_cast<u64>(p) * p;
        const u64 expected = (q2 + 1) * (q2 + p + 1);
        if (bases.size() != expected)
            throw std::runtime_error("two-space domain size mismatch");
    }

    u64 rrefKey(u32 first, u32 second) const {
        u32 rows[2][4];
        for (u32 column = 0; column < 4; ++column) {
            rows[0][column] = first % p;
            first /= p;
            rows[1][column] = second % p;
            second /= p;
        }
        u32 rank = 0;
        for (u32 column = 0; column < 4 && rank < 2; ++column) {
            u32 pivot = rank;
            while (pivot < 2 && rows[pivot][column] == 0) ++pivot;
            if (pivot == 2) continue;
            if (pivot != rank)
                for (u32 j = 0; j < 4; ++j)
                    std::swap(rows[pivot][j], rows[rank][j]);
            const u32 inverse = powMod(rows[rank][column], p - 2, p);
            for (u32 j = column; j < 4; ++j)
                rows[rank][j] = static_cast<u32>(mulMod(rows[rank][j], inverse, p));
            for (u32 row = 0; row < 2; ++row) {
                if (row == rank) continue;
                const u32 factor = rows[row][column];
                if (factor == 0) continue;
                for (u32 j = column; j < 4; ++j) {
                    const u32 subtract = static_cast<u32>(mulMod(
                        factor, rows[rank][j], p));
                    rows[row][j] = (rows[row][j] + p - subtract) % p;
                }
            }
            ++rank;
        }
        if (rank != 2) throw std::runtime_error("singular two-space image");
        u64 key = 0, place = 1;
        for (u32 row = 0; row < 2; ++row)
            for (u32 column = 0; column < 4; ++column) {
                key += static_cast<u64>(rows[row][column]) * place;
                place *= p;
            }
        return key;
    }
};

const TwoSpaceDomain& cachedTwoSpaceDomain(u32 p) {
    static std::map<u32, std::unique_ptr<TwoSpaceDomain>> cache;
    auto found = cache.find(p);
    if (found == cache.end())
        found = cache.emplace(p, std::make_unique<TwoSpaceDomain>(p)).first;
    return *found->second;
}

ProjectiveOrbitFingerprint twoSpaceOrbitSizeFingerprint(
        const std::vector<GeneratorTable>& tables, const PackedVectorOps& ops) {
    if (ops.n != 4) throw std::runtime_error("two-space fingerprint requires n=4");
    const TwoSpaceDomain& domain = cachedTwoSpaceDomain(ops.p);
    std::vector<u32> orbit_id(domain.bases.size(), UNASSIGNED);
    std::vector<u32> orbit_sizes;
    std::vector<u32> queue;
    for (u32 seed = 0; seed < domain.bases.size(); ++seed) {
        if (orbit_id[seed] != UNASSIGNED) continue;
        const u32 orbit = static_cast<u32>(orbit_sizes.size());
        orbit_sizes.push_back(0);
        orbit_id[seed] = orbit;
        queue.clear();
        queue.push_back(seed);
        for (std::size_t head = 0; head < queue.size(); ++head) {
            const u32 space = queue[head];
            ++orbit_sizes[orbit];
            const auto& basis = domain.bases[space];
            for (const auto& table : tables) {
                const u32 first = applyGeneratorToPacked(basis.first, ops, table);
                const u32 second = applyGeneratorToPacked(basis.second, ops, table);
                const u64 key = domain.rrefKey(first, second);
                const auto found = domain.index.find(key);
                if (found == domain.index.end())
                    throw std::runtime_error("two-space image absent from domain");
                const u32 image = found->second;
                if (orbit_id[image] == UNASSIGNED) {
                    orbit_id[image] = orbit;
                    queue.push_back(image);
                }
            }
        }
    }
    std::sort(orbit_sizes.begin(), orbit_sizes.end());
    u64 first = 1469598103934665603ULL;
    u64 second = 0x9e3779b97f4a7c15ULL;
    for (u32 size : orbit_sizes) {
        first ^= static_cast<u64>(size);
        first *= 1099511628211ULL;
        second ^= static_cast<u64>(size) + 0x9e3779b97f4a7c15ULL +
                  (second << 6U) + (second >> 2U);
    }
    first ^= static_cast<u64>(orbit_sizes.size());
    first *= 1099511628211ULL;
    second ^= static_cast<u64>(orbit_sizes.size()) * 0xbf58476d1ce4e5b9ULL;
    return {static_cast<u32>(orbit_sizes.size()), first, second};
}

u64 flagSize(const std::vector<std::uint8_t>& flags, const OrbitData& orbits) {
    u64 size = 0;
    for (std::size_t i = 0; i < flags.size(); ++i)
        if (flags[i]) size += orbits.size[i];
    return size;
}

std::vector<u32> flaggedRepresentatives(const std::vector<std::uint8_t>& flags,
                                        const OrbitData& orbits) {
    std::vector<u32> result;
    for (std::size_t i = 0; i < flags.size(); ++i)
        if (flags[i]) result.push_back(orbits.representative[i]);
    return result;
}

std::vector<u32> flaggedVectors(const std::vector<std::uint8_t>& flags,
                                const OrbitData& orbits) {
    std::vector<u32> result;
    result.reserve(static_cast<std::size_t>(flagSize(flags, orbits)));
    for (u32 v = 0; v < orbits.id.size(); ++v)
        if (flags[orbits.id[v]]) result.push_back(v);
    return result;
}

struct TwoSumWitnessResult {
    bool completed = true;
    bool full = true;
    u64 probes = 0;
    u64 target_orbits_checked = 0;
    u32 first_missing = 0;
};

// R is H-invariant and R = -R.  Hence R + R is H-invariant, and an orbit
// representative x lies in R + R exactly when x + r lies in R for some
// r in R.  Testing membership and stopping at the first witness is much
// cheaper than materialising every pair sum when only the predicate
// R + R = V is required.
TwoSumWitnessResult twoSumFullByOrbitWitness(
        const std::vector<std::uint8_t>& regular,
        const OrbitData& orbits, const PackedVectorOps& ops,
        u64 pair_budget) {
    TwoSumWitnessResult result;
    const std::vector<u32> regular_vectors = flaggedVectors(regular, orbits);
    for (u32 target : orbits.representative) {
        ++result.target_orbits_checked;
        bool covered = false;
        for (u32 r : regular_vectors) {
            if (pair_budget != 0 && result.probes == pair_budget) {
                result.completed = false;
                return result;
            }
            ++result.probes;
            if (regular[orbits.id[ops.add(target, r)]]) {
                covered = true;
                break;
            }
        }
        if (!covered) {
            result.full = false;
            result.first_missing = target;
            return result;
        }
    }
    return result;
}

bool sumOrbitFlags(const std::vector<std::uint8_t>& left_flags,
                   const std::vector<std::uint8_t>& right_flags,
                   const OrbitData& orbits, const PackedVectorOps& ops,
                   u64 pair_budget, std::vector<std::uint8_t>& output,
                   u64& pair_count) {
    const std::vector<u32> left_reps = flaggedRepresentatives(left_flags, orbits);
    const std::vector<u32> right_reps = flaggedRepresentatives(right_flags, orbits);
    const u64 left_size = flagSize(left_flags, orbits);
    const u64 right_size = flagSize(right_flags, orbits);
    const __uint128_t left_cost = static_cast<__uint128_t>(left_reps.size()) * right_size;
    const __uint128_t right_cost = static_cast<__uint128_t>(right_reps.size()) * left_size;
    const bool representatives_from_left = left_cost <= right_cost;
    const __uint128_t cost = std::min(left_cost, right_cost);
    if (pair_budget != 0 && cost > pair_budget) {
        pair_count = cost > std::numeric_limits<u64>::max()
                         ? std::numeric_limits<u64>::max()
                         : static_cast<u64>(cost);
        return false;
    }

    const std::vector<u32>& reps = representatives_from_left ? left_reps : right_reps;
    const std::vector<u32> all =
        flaggedVectors(representatives_from_left ? right_flags : left_flags, orbits);
    output.assign(orbits.representative.size(), 0);
    pair_count = static_cast<u64>(cost);
    for (u32 representative : reps)
        for (u32 v : all)
            output[orbits.id[ops.add(representative, v)]] = 1;
    return true;
}

void printPrefix(const Action& a, const PackedVectorOps& ops, double seconds) {
    std::cout << "{\"schema\":\"AFFINE_SAXL_RESULT_V1\""
              << ",\"label\":\"" << jsonEscape(a.label) << "\""
              << ",\"p\":" << a.p
              << ",\"n\":" << a.n
              << ",\"degree\":" << ops.vectors
              << ",\"group_order\":";
    if (a.order_fits_u64) std::cout << a.order;
    else std::cout << "null";
    std::cout << ",\"group_order_decimal\":\"" << a.order_decimal << "\""
              << ",\"generators\":" << a.gens.size()
              << ",\"orientation\":\"" << (a.row_action ? "row" : "column") << "\""
              << ",\"add_mode\":\""
              << (ops.xor_addition ? "xor" : (ops.table_addition ? "split_table" : "digitwise"))
              << "\""
              << ",\"seconds\":" << std::fixed << std::setprecision(6) << seconds;
}

void emitEarlyResult(const Action& a, const PackedVectorOps& ops, const std::string& status,
                     const std::string& base_size, const std::string& reason,
                     std::chrono::steady_clock::time_point started) {
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    printPrefix(a, ops, seconds);
    std::cout << ",\"status\":\"" << status << "\""
              << ",\"base_size\":\"" << base_size << "\""
              << ",\"reason\":\"" << jsonEscape(reason) << "\""
              << ",\"regular_orbits\":0,\"regular_vectors\":0"
              << ",\"two_sum_full\":null,\"diameter\":null}\n";
}

void processAction(const Action& a, const Options& options) {
    const auto started = std::chrono::steady_clock::now();
    if (!isPrime(a.p)) throw std::runtime_error("p is not prime");
    if (a.n > 63) throw std::runtime_error("dimension exceeds implementation limit 63");
    for (std::size_t i = 0; i < a.gens.size(); ++i)
        if (!invertible(a.gens[i], a.n, a.p))
            throw std::runtime_error("generator " + std::to_string(i + 1) + " is singular");

    u64 vectors = 1;
    for (u32 i = 0; i < a.n; ++i) {
        if (vectors > std::numeric_limits<u32>::max() / a.p)
            throw std::runtime_error("p^n exceeds the 32-bit packed-vector limit");
        vectors *= a.p;
    }
    if (vectors > options.max_vectors)
        throw std::runtime_error("p^n exceeds --max-vectors");
    const PackedVectorOps ops(a.p, a.n, options.add_table_mib);

    if (options.orbit_fingerprint_only) {
        if (a.gens.empty() && a.order_decimal != "1")
            throw std::runtime_error("nontrivial H has no generators");
        const auto tables = makeGeneratorTables(a, ops);
        const OrbitData orbits = enumerateOrbits(a, ops, tables, false, false);
        const auto fingerprint = orbitSizeFingerprint(orbits);
        const auto projective = projectiveOrbitSizeFingerprint(orbits, ops);
        ProjectiveOrbitFingerprint dual_fingerprint;
        std::pair<u64, u64> dual_vector_fingerprint{0, 0};
        u32 dual_orbits_count = 0;
        if (options.dual_orbit_fingerprint) {
            Action dual = a;
            for (std::vector<u32>& generator : dual.gens)
                generator = inverseTranspose(generator, dual.n, dual.p);
            const auto dual_tables = makeGeneratorTables(dual, ops);
            const OrbitData dual_orbits = enumerateOrbits(dual, ops, dual_tables, false, false);
            dual_orbits_count = static_cast<u32>(dual_orbits.size.size());
            dual_vector_fingerprint = orbitSizeFingerprint(dual_orbits);
            dual_fingerprint = projectiveOrbitSizeFingerprint(dual_orbits, ops);
        }
        const bool have_two_spaces = options.two_space_fingerprint && a.n == 4;
        const ProjectiveOrbitFingerprint two_spaces = have_two_spaces
            ? twoSpaceOrbitSizeFingerprint(tables, ops)
            : ProjectiveOrbitFingerprint{};
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        printPrefix(a, ops, seconds);
        std::cout << ",\"status\":\"orbit_fingerprint\""
                  << ",\"orbits\":" << orbits.size.size()
                  << ",\"orbit_size_hash_1\":\"" << std::hex
                  << std::setw(16) << std::setfill('0') << fingerprint.first
                  << "\",\"orbit_size_hash_2\":\""
                  << std::setw(16) << fingerprint.second << std::dec
                  << std::setfill(' ') << "\""
                  << ",\"projective_orbits\":" << projective.orbits
                  << ",\"projective_orbit_size_hash_1\":\"" << std::hex
                  << std::setw(16) << std::setfill('0') << projective.first
                  << "\",\"projective_orbit_size_hash_2\":\""
                  << std::setw(16) << projective.second << std::dec
                  << std::setfill(' ') << "\""
                  << ",\"dual_orbits\":";
        if (options.dual_orbit_fingerprint) {
            std::cout << dual_orbits_count
                      << ",\"dual_orbit_size_hash_1\":\"" << std::hex
                      << std::setw(16) << std::setfill('0') << dual_vector_fingerprint.first
                      << "\",\"dual_orbit_size_hash_2\":\""
                      << std::setw(16) << dual_vector_fingerprint.second << std::dec
                      << std::setfill(' ') << "\""
                      << ",\"dual_projective_orbits\":" << dual_fingerprint.orbits
                      << ",\"dual_projective_orbit_size_hash_1\":\"" << std::hex
                      << std::setw(16) << std::setfill('0') << dual_fingerprint.first
                      << "\",\"dual_projective_orbit_size_hash_2\":\""
                      << std::setw(16) << dual_fingerprint.second << std::dec
                      << std::setfill(' ') << "\"";
        } else {
            std::cout << "null,\"dual_orbit_size_hash_1\":null"
                      << ",\"dual_orbit_size_hash_2\":null"
                      << ",\"dual_projective_orbits\":null"
                      << ",\"dual_projective_orbit_size_hash_1\":null"
                      << ",\"dual_projective_orbit_size_hash_2\":null";
        }
        std::cout
                  << ",\"two_space_orbits\":";
        if (have_two_spaces) {
            std::cout << two_spaces.orbits
                      << ",\"two_space_orbit_size_hash_1\":\"" << std::hex
                      << std::setw(16) << std::setfill('0') << two_spaces.first
                      << "\",\"two_space_orbit_size_hash_2\":\""
                      << std::setw(16) << two_spaces.second << std::dec
                      << std::setfill(' ') << "\"";
        } else {
            std::cout << "null,\"two_space_orbit_size_hash_1\":null"
                      << ",\"two_space_orbit_size_hash_2\":null";
        }
        std::cout << "}\n";
        return;
    }

    if (a.order_fits_u64 && a.order == 1) {
        emitEarlyResult(a, ops, "trivial_stabiliser", "1",
                        "H is trivial, so the affine action has base size 1", started);
        return;
    }
    if (!a.order_fits_u64 || a.order > vectors - 1) {
        emitEarlyResult(a, ops, "order_obstruction", ">2",
                        "a regular H-orbit cannot fit in V minus {0}", started);
        return;
    }
    // Every nonidentity linear transformation fixes at most p^(n-1)
    // vectors.  Thus 2(|H|-1) < p implies that the union of all
    // nonidentity fixed-point spaces has size < |V|/2.  More than half of V
    // is consequently regular, and any two translates of the regular set
    // meet, proving R + R = V without enumerating a single vector.
    if (options.diameter_ge3_only && 2 * (a.order - 1) < a.p) {
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        printPrefix(a, ops, seconds);
        std::cout << ",\"status\":\"fixed_space_union_bound\",\"base_size\":\"2\""
                  << ",\"reason\":\"2(|H|-1) < p, so more than half of V is regular\""
                  << ",\"regular_orbits\":null,\"regular_vectors_lower_bound\":"
                  << (vectors / 2 + 1)
                  << ",\"two_sum_full\":true,\"diameter\":null"
                  << ",\"diameter_upper_bound\":2,\"diameter_at_least_3\":false}\n";
        return;
    }
    if (a.gens.empty()) throw std::runtime_error("nontrivial H has no generators");

    std::cerr << "[affine-saxl] " << a.label << ": building generator tables for "
              << vectors << " vectors\n";
    const auto tables = makeGeneratorTables(a, ops);
    OrbitData orbits = enumerateOrbits(a, ops, tables, options.diameter_ge3_only);

    if (!orbits.complete) {
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        printPrefix(a, ops, seconds);
        std::cout << ",\"status\":\"density_certificate\",\"base_size\":\"2\""
                  << ",\"orbits_partial\":" << orbits.size.size()
                  << ",\"regular_orbits_lower_bound\":" << orbits.regular_orbits_seen
                  << ",\"regular_vectors_lower_bound\":" << orbits.regular_vectors_seen
                  << ",\"two_sum_full\":true,\"diameter\":null"
                  << ",\"diameter_upper_bound\":2,\"diameter_at_least_3\":false}\n";
        return;
    }

    std::vector<std::uint8_t> regular(orbits.representative.size(), 0);
    u64 regular_orbits = 0, regular_vectors = 0;
    for (std::size_t i = 0; i < orbits.size.size(); ++i) {
        if (orbits.size[i] == a.order) {
            regular[i] = 1;
            ++regular_orbits;
            regular_vectors += orbits.size[i];
        }
    }
    if (regular_vectors == 0) {
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        printPrefix(a, ops, seconds);
        std::cout << ",\"status\":\"no_regular_vector\",\"base_size\":\">2\""
                  << ",\"orbits\":" << orbits.size.size()
                  << ",\"regular_orbits\":0,\"regular_vectors\":0"
                  << ",\"two_sum_full\":null,\"diameter\":null}\n";
        return;
    }

    std::vector<std::uint8_t> reached(orbits.representative.size(), 0);
    std::vector<std::uint8_t> frontier = regular;
    reached[orbits.id[0]] = 1;
    for (std::size_t i = 0; i < regular.size(); ++i)
        if (regular[i]) reached[i] = 1;
    u64 reached_size = flagSize(reached, orbits);
    std::vector<u64> layer_sizes = {1, regular_vectors};
    std::vector<u64> pair_counts;
    bool two_sum_known = false, two_sum_full = false;
    u64 two_sum_size = 0;
    u32 first_two_sum_missing = 0;
    bool have_first_two_sum_missing = false;
    u32 first_distance_ge3 = 0;
    bool have_first_distance_ge3 = false;
    int diameter = regular_vectors == vectors - 1 ? 1 : -1;
    std::string status = "exact";

    if (2 * regular_vectors > vectors) {
        two_sum_known = true;
        two_sum_full = true;
        two_sum_size = vectors;
        if (diameter < 0) {
            diameter = 2;
            layer_sizes.push_back(vectors - reached_size);
            reached_size = vectors;
        }
    }

    if (options.diameter_ge3_only && diameter < 0) {
        const TwoSumWitnessResult check = twoSumFullByOrbitWitness(
            regular, orbits, ops, options.pair_budget);
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        printPrefix(a, ops, seconds);
        if (!check.completed) {
            std::cout << ",\"status\":\"pair_budget_exceeded\",\"base_size\":\"2\""
                      << ",\"regular_orbits\":" << regular_orbits
                      << ",\"regular_vectors\":" << regular_vectors
                      << ",\"two_sum_full\":null,\"diameter\":null"
                      << ",\"diameter_at_least_3\":false"
                      << ",\"target_orbits_checked\":" << check.target_orbits_checked
                      << ",\"membership_probes\":" << check.probes << "}\n";
        } else if (check.full) {
            std::cout << ",\"status\":\"two_sum_certificate\",\"base_size\":\"2\""
                      << ",\"regular_orbits\":" << regular_orbits
                      << ",\"regular_vectors\":" << regular_vectors
                      << ",\"two_sum_full\":true,\"diameter\":null"
                      << ",\"diameter_upper_bound\":2,\"diameter_at_least_3\":false"
                      << ",\"target_orbits_checked\":" << check.target_orbits_checked
                      << ",\"membership_probes\":" << check.probes << "}\n";
        } else {
            std::cout << ",\"status\":\"diameter_at_least_3\",\"base_size\":\"2\""
                      << ",\"regular_orbits\":" << regular_orbits
                      << ",\"regular_vectors\":" << regular_vectors
                      << ",\"two_sum_full\":false,\"diameter\":null"
                      << ",\"diameter_at_least_3\":true"
                      << ",\"first_two_sum_missing\":" << check.first_missing
                      << ",\"first_two_sum_missing_vector\":"
                      << ops.digits(check.first_missing)
                      << ",\"first_distance_ge3\":" << check.first_missing
                      << ",\"first_distance_ge3_vector\":"
                      << ops.digits(check.first_missing)
                      << ",\"target_orbits_checked\":" << check.target_orbits_checked
                      << ",\"membership_probes\":" << check.probes << "}\n";
        }
        return;
    }

    int level = 1;
    while (diameter < 0 && reached_size < vectors && level < options.max_diameter) {
        std::vector<std::uint8_t> sums;
        u64 pairs = 0;
        if (!sumOrbitFlags(regular, frontier, orbits, ops, options.pair_budget, sums, pairs)) {
            status = "pair_budget_exceeded";
            pair_counts.push_back(pairs);
            break;
        }
        pair_counts.push_back(pairs);
        if (level == 1) {
            two_sum_known = true;
            two_sum_size = flagSize(sums, orbits);
            two_sum_full = two_sum_size == vectors;
            if (!two_sum_full) {
                for (std::size_t i = 0; i < sums.size(); ++i) {
                    if (!sums[i]) {
                        first_two_sum_missing = orbits.representative[i];
                        have_first_two_sum_missing = true;
                        break;
                    }
                }
            }
        }

        std::vector<std::uint8_t> next(orbits.representative.size(), 0);
        u64 next_size = 0;
        for (std::size_t i = 0; i < sums.size(); ++i) {
            if (sums[i] && !reached[i]) {
                next[i] = 1;
                reached[i] = 1;
                next_size += orbits.size[i];
            }
        }
        ++level;
        layer_sizes.push_back(next_size);
        reached_size += next_size;
        if (reached_size == vectors) {
            diameter = level;
            break;
        }
        if (options.diameter_ge3_only && level == 2) {
            status = "diameter_at_least_3";
            for (std::size_t i = 0; i < reached.size(); ++i) {
                if (!reached[i]) {
                    first_distance_ge3 = orbits.representative[i];
                    have_first_distance_ge3 = true;
                    break;
                }
            }
            break;
        }
        if (next_size == 0) {
            status = "disconnected";
            break;
        }
        frontier.swap(next);
    }
    if (diameter < 0 && status == "exact") status = "diameter_cap_reached";

    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    printPrefix(a, ops, seconds);
    std::cout << ",\"status\":\"" << status << "\",\"base_size\":\"2\""
              << ",\"orbits\":" << orbits.size.size()
              << ",\"regular_orbits\":" << regular_orbits
              << ",\"regular_vectors\":" << regular_vectors
              << ",\"regular_density\":" << std::setprecision(12)
              << static_cast<double>(regular_vectors) / static_cast<double>(vectors)
              << ",\"two_sum_full\":";
    if (two_sum_known) std::cout << (two_sum_full ? "true" : "false");
    else std::cout << "null";
    std::cout << ",\"two_sum_size\":";
    if (two_sum_known) std::cout << two_sum_size;
    else std::cout << "null";
    std::cout << ",\"first_two_sum_missing\":";
    if (have_first_two_sum_missing) std::cout << first_two_sum_missing;
    else std::cout << "null";
    std::cout << ",\"first_two_sum_missing_vector\":";
    if (have_first_two_sum_missing) std::cout << ops.digits(first_two_sum_missing);
    else std::cout << "null";
    std::cout << ",\"diameter\":";
    if (diameter >= 0) std::cout << diameter;
    else std::cout << "null";
    std::cout << ",\"diameter_at_least_3\":"
              << ((status == "diameter_at_least_3" || diameter >= 3) ? "true" : "false")
              << ",\"first_distance_ge3\":";
    if (have_first_distance_ge3) std::cout << first_distance_ge3;
    else std::cout << "null";
    std::cout << ",\"first_distance_ge3_vector\":";
    if (have_first_distance_ge3) std::cout << ops.digits(first_distance_ge3);
    else std::cout << "null";
    std::cout << ",\"reached_vectors\":" << reached_size << ",\"layer_sizes\":[";
    for (std::size_t i = 0; i < layer_sizes.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << layer_sizes[i];
    }
    std::cout << "],\"sum_pairs\":[";
    for (std::size_t i = 0; i < pair_counts.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << pair_counts[i];
    }
    std::cout << "]}\n";
}

void usage(const char* program) {
    std::cerr << "usage: " << program
              << " [input|-] [--max-vectors=N] [--add-table-mib=N]"
                 " [--pair-budget=N] [--max-diameter=N] [--diameter-ge3-only]"
                 " [--orbit-fingerprint-only] [--dual-orbit-fingerprint]"
                 " [--two-space-fingerprint]\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        Options options;
        std::string path = "-";
        bool have_path = false;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg.rfind("--max-vectors=", 0) == 0) {
                options.max_vectors = parseU64(arg.substr(14), "--max-vectors");
            } else if (arg.rfind("--add-table-mib=", 0) == 0) {
                options.add_table_mib = parseU64(arg.substr(16), "--add-table-mib");
            } else if (arg.rfind("--pair-budget=", 0) == 0) {
                options.pair_budget = parseU64(arg.substr(14), "--pair-budget");
            } else if (arg.rfind("--max-diameter=", 0) == 0) {
                const u64 value = parseU64(arg.substr(15), "--max-diameter");
                if (value == 0 || value > std::numeric_limits<int>::max())
                    throw std::runtime_error("--max-diameter must be positive");
                options.max_diameter = static_cast<int>(value);
            } else if (arg == "--diameter-ge3-only") {
                options.diameter_ge3_only = true;
            } else if (arg == "--orbit-fingerprint-only") {
                options.orbit_fingerprint_only = true;
            } else if (arg == "--dual-orbit-fingerprint") {
                options.dual_orbit_fingerprint = true;
            } else if (arg == "--two-space-fingerprint") {
                options.two_space_fingerprint = true;
            } else if (arg == "--help" || arg == "-h") {
                usage(argv[0]);
                return 0;
            } else if (!have_path) {
                path = arg;
                have_path = true;
            } else {
                throw std::runtime_error("unexpected argument: " + arg);
            }
        }

        if (options.two_space_fingerprint && !options.orbit_fingerprint_only)
            throw std::runtime_error(
                "--two-space-fingerprint requires --orbit-fingerprint-only");
        if (options.dual_orbit_fingerprint && !options.orbit_fingerprint_only)
            throw std::runtime_error(
                "--dual-orbit-fingerprint requires --orbit-fingerprint-only");

        std::ifstream file;
        std::istream* input = &std::cin;
        if (path != "-") {
            file.open(path);
            if (!file) throw std::runtime_error("cannot open input file: " + path);
            input = &file;
        }
        const std::vector<Action> actions = readActions(*input);
        if (actions.empty()) throw std::runtime_error("input contains no actions");
        for (const Action& action : actions) {
            try {
                processAction(action, options);
            } catch (const std::exception& error) {
                std::cout << "{\"schema\":\"AFFINE_SAXL_RESULT_V1\",\"label\":\""
                          << jsonEscape(action.label) << "\",\"status\":\"error\",\"error\":\""
                          << jsonEscape(error.what()) << "\"}\n";
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "affine_saxl_engine: " << error.what() << '\n';
        usage(argv[0]);
        return 2;
    }
}

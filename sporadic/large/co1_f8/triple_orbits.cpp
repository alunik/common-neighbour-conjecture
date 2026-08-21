// Enumerate the orbits of a pair stabilizer on the third F_2^24 coordinate.
// Vectors are packed into 24-bit words.  Three byte lookup tables make each
// matrix-vector product inexpensive, while the visited set occupies 2 MiB.

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::uint32_t dimension = 24;
constexpr std::uint32_t vector_count = 1u << dimension;
constexpr std::uint32_t vector_mask = vector_count - 1;
using Matrix = std::array<std::uint32_t, dimension>;

struct PairStabilizer {
    std::uint32_t pair_index{};
    std::uint32_t a{};
    std::uint32_t b{};
    std::uint64_t first_orbit_size{};
    std::uint64_t b_orbit_size{};
    std::uint64_t pair_stabilizer_order{};
    std::vector<Matrix> matrices;
};

struct PackedMatrix {
    Matrix rows{};
    std::array<std::array<std::uint32_t, 256>, 3> lookup{};

    explicit PackedMatrix(const Matrix& input) : rows(input) {
        for (std::uint32_t chunk = 0; chunk < 3; ++chunk) {
            for (std::uint32_t value = 0; value < 256; ++value) {
                std::uint32_t result = 0;
                for (std::uint32_t bit = 0; bit < 8; ++bit)
                    if ((value >> bit) & 1u)
                        result ^= rows[8 * chunk + bit];
                lookup[chunk][value] = result;
            }
        }
    }

    std::uint32_t image(std::uint32_t value) const noexcept {
        return lookup[0][value & 0xffu] ^
               lookup[1][(value >> 8) & 0xffu] ^
               lookup[2][(value >> 16) & 0xffu];
    }
};

class PackedBits {
  public:
    explicit PackedBits(std::uint32_t size)
        : words_((static_cast<std::size_t>(size) + 63) / 64, 0) {}

    bool test(std::uint32_t index) const noexcept {
        return (words_[index >> 6] >> (index & 63u)) & 1u;
    }

    bool set_if_new(std::uint32_t index) noexcept {
        auto& word = words_[index >> 6];
        const std::uint64_t bit = 1ULL << (index & 63u);
        if (word & bit) return false;
        word |= bit;
        return true;
    }

  private:
    std::vector<std::uint64_t> words_;
};

std::uint32_t rank(Matrix rows) {
    std::uint32_t rank = 0;
    for (std::uint32_t column = 0;
         column < dimension && rank < dimension; ++column) {
        std::uint32_t pivot = rank;
        while (pivot < dimension &&
               ((rows[pivot] >> column) & 1u) == 0)
            ++pivot;
        if (pivot == dimension) continue;
        std::swap(rows[rank], rows[pivot]);
        for (std::uint32_t row = 0; row < dimension; ++row)
            if (row != rank && ((rows[row] >> column) & 1u))
                rows[row] ^= rows[rank];
        ++rank;
    }
    return rank;
}

std::uint32_t direct_image(
    std::uint32_t value, const Matrix& matrix) noexcept {
    std::uint32_t result = 0;
    while (value != 0) {
        const auto bit =
            static_cast<std::uint32_t>(__builtin_ctz(value));
        result ^= matrix[bit];
        value &= value - 1;
    }
    return result;
}

PairStabilizer read_input(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open " + path);
    std::string token;
    input >> token;
    if (token != "CO1_F8_PAIR_STABILIZER_V1")
        throw std::runtime_error("bad pair-stabilizer header");

    PairStabilizer result;
    std::size_t generator_count = 0;
    input >> token >> result.pair_index;
    if (token != "pair_index" || result.pair_index == 0)
        throw std::runtime_error("bad pair index");
    input >> token >> result.a;
    if (token != "a" || result.a > vector_mask)
        throw std::runtime_error("bad first vector");
    input >> token >> result.b;
    if (token != "b" || result.b > vector_mask)
        throw std::runtime_error("bad second vector");
    input >> token >> result.first_orbit_size;
    if (token != "first_orbit_size" ||
        result.first_orbit_size == 0)
        throw std::runtime_error("bad first-orbit size");
    input >> token >> result.b_orbit_size;
    if (token != "b_orbit_size" || result.b_orbit_size == 0)
        throw std::runtime_error("bad second-orbit size");
    input >> token >> result.pair_stabilizer_order;
    if (token != "pair_stabilizer_order" ||
        result.pair_stabilizer_order == 0)
        throw std::runtime_error("bad pair-stabilizer order");
    input >> token >> generator_count;
    if (token != "generator_count")
        throw std::runtime_error("missing generator count");
    result.matrices.resize(generator_count);
    for (auto& matrix : result.matrices) {
        input >> token;
        if (token != "matrix")
            throw std::runtime_error("missing generator matrix");
        for (auto& row : matrix) {
            input >> row;
            if (row > vector_mask)
                throw std::runtime_error("matrix row out of range");
        }
        if (rank(matrix) != dimension)
            throw std::runtime_error("singular generator matrix");
        if (direct_image(result.a, matrix) != result.a ||
            direct_image(result.b, matrix) != result.b)
            throw std::runtime_error(
                "generator matrix does not fix pair");
    }
    if (input >> token)
        throw std::runtime_error("trailing pair-stabilizer input");
    if (result.matrices.empty() &&
        result.pair_stabilizer_order != 1)
        throw std::runtime_error(
            "nontrivial stabilizer has no generators");
    constexpr std::uint64_t co1_order = 4157776806543360000ULL;
    if (result.first_orbit_size > co1_order ||
        result.b_orbit_size >
            co1_order / result.first_orbit_size ||
        result.first_orbit_size * result.b_orbit_size >
            co1_order / result.pair_stabilizer_order ||
        result.first_orbit_size * result.b_orbit_size *
            result.pair_stabilizer_order != co1_order)
        throw std::runtime_error("pair orbit-stabilizer identity fails");
    return result;
}

struct CensusSummary {
    std::uint64_t orbit_count{};
    std::uint64_t regular_orbit_count{};
    std::uint64_t mass{};
};

CensusSummary enumerate_orbits(
    const PairStabilizer& pair,
    std::uint32_t count,
    std::ostream* output) {
    std::vector<PackedMatrix> generators;
    generators.reserve(pair.matrices.size());
    for (const auto& matrix : pair.matrices)
        generators.emplace_back(matrix);

    PackedBits seen(count);
    std::vector<std::uint32_t> queue;
    queue.reserve(1u << 20);
    CensusSummary summary;

    for (std::uint32_t seed = 0; seed < count; ++seed) {
        if (seen.test(seed)) continue;
        queue.clear();
        queue.push_back(seed);
        seen.set_if_new(seed);
        for (std::size_t head = 0; head < queue.size(); ++head) {
            for (const auto& generator : generators) {
                const auto target = generator.image(queue[head]);
                if (target >= count)
                    throw std::runtime_error(
                        "generator leaves census domain");
                if (seen.set_if_new(target)) queue.push_back(target);
            }
        }
        const auto orbit_size =
            static_cast<std::uint64_t>(queue.size());
        if (orbit_size == 0 ||
            pair.pair_stabilizer_order % orbit_size != 0)
            throw std::runtime_error(
                "nonintegral triple-stabilizer order");
        const auto triple_stabilizer_order =
            pair.pair_stabilizer_order / orbit_size;
        ++summary.orbit_count;
        if (triple_stabilizer_order == 1)
            ++summary.regular_orbit_count;
        summary.mass += orbit_size;
        if (output != nullptr) {
            *output << pair.pair_index << '\t'
                    << pair.a << '\t'
                    << pair.b << '\t'
                    << seed << '\t'
                    << pair.first_orbit_size << '\t'
                    << pair.b_orbit_size << '\t'
                    << orbit_size << '\t'
                    << triple_stabilizer_order << '\n';
        }
    }
    if (summary.mass != count)
        throw std::runtime_error("third-coordinate mass mismatch");
    return summary;
}

Matrix swap_first_two_matrix() {
    Matrix matrix{};
    for (std::uint32_t bit = 0; bit < dimension; ++bit)
        matrix[bit] = 1u << bit;
    std::swap(matrix[0], matrix[1]);
    return matrix;
}

void self_test() {
    PairStabilizer pair;
    pair.pair_index = 1;
    pair.first_orbit_size = 1;
    pair.b_orbit_size = 1;
    pair.pair_stabilizer_order = 2;
    pair.matrices.push_back(swap_first_two_matrix());
    const auto summary = enumerate_orbits(pair, 16, nullptr);
    if (summary.orbit_count != 12 ||
        summary.regular_orbit_count != 4 ||
        summary.mass != 16)
        throw std::runtime_error("small orbit self-test failed");
    const PackedMatrix packed(pair.matrices[0]);
    for (std::uint32_t value = 0; value < 16; ++value)
        if (packed.image(value) !=
            direct_image(value, pair.matrices[0]))
            throw std::runtime_error("packed image self-test failed");
    std::cout
        << "Self-test: 12 orbits, 4 regular orbits, mass 16.\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            self_test();
            return 0;
        }
        if (argc != 3) {
            std::cerr
                << "usage: triple_orbits PAIR_STABILIZER OUTPUT_TSV\n"
                << "       triple_orbits --self-test\n";
            return 1;
        }
        const auto pair = read_input(argv[1]);
        std::ofstream output(argv[2]);
        if (!output)
            throw std::runtime_error(
                "cannot open output " + std::string(argv[2]));
        output
            << "pair_index\ta\tb\tc\tfirst_orbit_size"
            << "\tb_orbit_size\tc_orbit_size"
            << "\ttriple_stabilizer_order\n";
        const auto summary =
            enumerate_orbits(pair, vector_count, &output);
        output.flush();
        if (!output)
            throw std::runtime_error("failed while writing output");
        std::cout << "Triple orbits above pair " << pair.pair_index
                  << ": " << summary.orbit_count
                  << " (regular " << summary.regular_orbit_count << ")\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "TRIPLE-ORBIT FAILURE: "
                  << error.what() << '\n';
        return 1;
    }
}

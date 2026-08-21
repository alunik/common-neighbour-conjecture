// Search for two regular vectors whose sum lies in each target orbit.
//
// The search is deterministic.  It moves representatives of regular orbits
// by short words in the standard Co1 generators and by powers of eta, then
// sends the complementary summands to the exact orbit classifier.

#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::uint32_t dimension = 24;
constexpr std::uint32_t vector_count = 1u << dimension;
constexpr std::uint32_t vector_mask = vector_count - 1;
constexpr std::uint32_t word_blocks = 64;
using Matrix = std::array<std::uint32_t, dimension>;
using Triple = std::array<std::uint32_t, 3>;

Matrix inverse(Matrix matrix) {
    std::array<std::uint64_t, dimension> augmented{};
    for (std::uint32_t row = 0; row < dimension; ++row)
        augmented[row] = matrix[row] |
            (std::uint64_t{1} << (dimension + row));
    for (std::uint32_t column = 0; column < dimension; ++column) {
        std::uint32_t pivot = column;
        while (pivot < dimension &&
               ((augmented[pivot] >> column) & 1u) == 0)
            ++pivot;
        if (pivot == dimension)
            throw std::runtime_error("singular matrix");
        std::swap(augmented[column], augmented[pivot]);
        for (std::uint32_t row = 0; row < dimension; ++row)
            if (row != column &&
                ((augmented[row] >> column) & 1u))
                augmented[row] ^= augmented[column];
    }
    Matrix result{};
    for (std::uint32_t row = 0; row < dimension; ++row) {
        if ((augmented[row] & vector_mask) !=
            (std::uint64_t{1} << row))
            throw std::runtime_error("inverse reduction failure");
        result[row] = static_cast<std::uint32_t>(
            augmented[row] >> dimension);
    }
    return result;
}

struct PackedMatrix {
    std::array<std::array<std::uint32_t, 256>, 3> lookup{};

    explicit PackedMatrix(const Matrix& matrix) {
        for (std::uint32_t chunk = 0; chunk < 3; ++chunk)
            for (std::uint32_t value = 0; value < 256; ++value) {
                std::uint32_t result = 0;
                for (std::uint32_t bit = 0; bit < 8; ++bit)
                    if ((value >> bit) & 1u)
                        result ^= matrix[8 * chunk + bit];
                lookup[chunk][value] = result;
            }
    }

    std::uint32_t image(std::uint32_t value) const noexcept {
        return lookup[0][value & 0xffu] ^
               lookup[1][(value >> 8) & 0xffu] ^
               lookup[2][(value >> 16) & 0xffu];
    }
};

std::vector<Matrix> read_primal_generators(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open Atlas input");
    std::string token;
    std::uint64_t value = 0;
    input >> token;
    if (token != "CO1_F8_ATLAS_INPUT_V1")
        throw std::runtime_error("bad Atlas header");
    input >> token >> value;
    if (token != "seed_point")
        throw std::runtime_error("missing seed point");
    input >> token >> value;
    if (token != "seed_mask")
        throw std::runtime_error("missing seed mask");
    for (int index = 0; index < 2; ++index) {
        input >> token;
        if (token != "dual_matrix")
            throw std::runtime_error("missing dual matrix");
        for (std::uint32_t row = 0; row < dimension; ++row)
            input >> value;
    }
    std::vector<Matrix> result(2);
    for (auto& matrix : result) {
        input >> token;
        if (token != "primal_matrix")
            throw std::runtime_error("missing primal matrix");
        for (auto& row : matrix) {
            input >> row;
            if (row > vector_mask)
                throw std::runtime_error("Atlas row out of range");
        }
    }
    for (int index = 0; index < 2; ++index) {
        input >> token;
        if (token != "permutation")
            throw std::runtime_error("missing permutation");
        for (std::uint32_t point = 0; point < 98280; ++point)
            input >> value;
    }
    if (input >> token)
        throw std::runtime_error("trailing Atlas input");
    return result;
}

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (true) {
        const auto tab = line.find('\t', start);
        fields.push_back(line.substr(
            start, tab == std::string::npos
                ? std::string::npos : tab - start));
        if (tab == std::string::npos) return fields;
        start = tab + 1;
    }
}

struct OrbitRow {
    std::uint32_t orbit_index{};
    Triple representative{};
    std::uint32_t eta_orbit_index{};
    bool maximal_regular{};
};

std::vector<OrbitRow> read_orbits(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot open classified orbits");
    std::string line;
    std::getline(input, line);
    if (line !=
        "orbit_index\tpair_index\ta\tb\tc\tfirst_orbit_size"
        "\tb_orbit_size\tc_orbit_size\ttriple_stabilizer_order"
        "\teta_orbit_index\tscalar_cycle_length\tmaximal_regular")
        throw std::runtime_error("bad classified-orbit header");
    std::vector<OrbitRow> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tabs(line);
        if (fields.size() != 12)
            throw std::runtime_error("bad classified-orbit row");
        OrbitRow row;
        row.orbit_index =
            static_cast<std::uint32_t>(std::stoul(fields[0]));
        for (int component = 0; component < 3; ++component) {
            const auto value = std::stoul(fields[2 + component]);
            if (value > vector_mask)
                throw std::runtime_error(
                    "orbit representative out of range");
            row.representative[component] =
                static_cast<std::uint32_t>(value);
        }
        row.eta_orbit_index =
            static_cast<std::uint32_t>(std::stoul(fields[9]));
        row.maximal_regular = std::stoul(fields[11]) == 1;
        if (row.orbit_index != result.size() + 1)
            throw std::runtime_error("noncontiguous orbit index");
        result.push_back(row);
    }
    if (result.size() != 9511)
        throw std::runtime_error("expected 9511 classified orbits");
    return result;
}

std::uint64_t splitmix64(std::uint64_t& state) noexcept {
    std::uint64_t z =
        (state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

Triple eta(Triple value) noexcept {
    return {value[2], value[0] ^ value[2], value[1]};
}

Triple group_word(
    Triple value, std::uint64_t code,
    const PackedMatrix& a,
    const PackedMatrix& b,
    const PackedMatrix& b_inverse) noexcept {
    for (std::uint32_t block = 0; block < word_blocks; ++block) {
        const auto& first =
            ((code >> block) & 1u) ? b_inverse : b;
        for (auto& component : value)
            component = first.image(component);
        for (auto& component : value)
            component = a.image(component);
    }
    return value;
}

void self_test() {
    Triple value{{1, 2, 4}};
    for (int power = 0; power < 7; ++power)
        value = eta(value);
    if (value != Triple{{1, 2, 4}})
        throw std::runtime_error("GF(8) self-test failed");
    std::cout << "Sum-search self-test passed\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            self_test();
            return 0;
        }
        if (argc != 6) {
            std::cerr
                << "usage: sum_queries ATLAS_INPUT "
                   "CLASSIFIED_ORBITS TRIALS QUERY_TSV CANDIDATE_TSV\n";
            return 1;
        }
        const auto generators = read_primal_generators(argv[1]);
        const PackedMatrix a(generators[0]);
        const PackedMatrix b(generators[1]);
        const PackedMatrix b_inverse(inverse(generators[1]));
        const auto orbits = read_orbits(argv[2]);
        const auto trials = std::stoul(argv[3]);
        if (trials == 0)
            throw std::runtime_error("trials must be positive");

        std::vector<std::size_t> regular;
        for (std::size_t index = 0; index < orbits.size(); ++index)
            if (orbits[index].maximal_regular)
                regular.push_back(index);
        if (regular.empty())
            throw std::runtime_error("no maximal regular orbits");

        std::ofstream queries(argv[4]);
        std::ofstream candidates(argv[5]);
        if (!queries || !candidates)
            throw std::runtime_error("cannot open witness output");
        queries << "query_id\tkind\tcandidate_id\ta\tb\tc\n";
        candidates
            << "candidate_id\ttarget_orbit\tsource_orbit"
            << "\tscalar_power\tword_code"
            << "\tx_a\tx_b\tx_c"
            << "\ty_a\ty_b\ty_c"
            << "\tz_a\tz_b\tz_c\n";

        std::uint64_t rng = 0xc01f8c105e5a17ULL;
        std::uint64_t candidate_id = 0;
        for (const auto& target : orbits)
            for (std::uint64_t trial = 0; trial < trials; ++trial) {
                ++candidate_id;
                const auto source_index =
                    regular[splitmix64(rng) % regular.size()];
                const auto& source = orbits[source_index];
                const auto word_code = splitmix64(rng);
                const auto scalar_power =
                    static_cast<std::uint32_t>(splitmix64(rng) % 7);
                auto y = group_word(
                    source.representative, word_code,
                    a, b, b_inverse);
                for (std::uint32_t power = 0;
                     power < scalar_power; ++power)
                    y = eta(y);
                Triple z{};
                for (int component = 0; component < 3; ++component)
                    z[component] =
                        target.representative[component] ^ y[component];

                queries << candidate_id
                        << "\tz\t" << candidate_id;
                for (const auto component : z)
                    queries << '\t' << component;
                queries << '\n';

                candidates << candidate_id
                           << '\t' << target.orbit_index
                           << '\t' << source.orbit_index
                           << '\t' << scalar_power
                           << '\t' << word_code;
                for (const auto component : target.representative)
                    candidates << '\t' << component;
                for (const auto component : y)
                    candidates << '\t' << component;
                for (const auto component : z)
                    candidates << '\t' << component;
                candidates << '\n';
            }

        std::cout << "Tried " << trials << " regular summands for each of "
                  << orbits.size() << " target orbits\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "WITNESS QUERY FAILURE: "
                  << error.what() << '\n';
        return 1;
    }
}

// Move arbitrary triples to the hierarchical Co1 orbit coordinates.
//
// Schreier forests first move coordinate a to one of four representatives,
// then move b under its stabilizer.  The remaining c coordinate can therefore
// be classified independently inside one of the 46 pair stabilizers.

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr std::uint32_t dimension = 24;
constexpr std::uint32_t vector_count = 1u << dimension;
constexpr std::uint32_t vector_mask = vector_count - 1;
constexpr std::uint32_t no_parent =
    std::numeric_limits<std::uint32_t>::max();
using Matrix = std::array<std::uint32_t, dimension>;
using Triple = std::array<std::uint32_t, 3>;

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
    for (std::uint32_t bit = 0; bit < dimension; ++bit) {
        const auto basis = std::uint32_t{1} << bit;
        if (direct_image(direct_image(basis, matrix), result) != basis ||
            direct_image(direct_image(basis, result), matrix) != basis)
            throw std::runtime_error("inverse verification failure");
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

struct OrbitDatum {
    std::uint32_t representative{};
    std::uint64_t size{};
};

class OrbitForest {
  public:
    OrbitForest(const std::vector<Matrix>& matrices,
                const std::vector<OrbitDatum>& expected)
        : parent_(vector_count, no_parent),
          parent_generator_(vector_count, 0) {
        if (matrices.empty())
            throw std::runtime_error("orbit forest has no generators");
        if (matrices.size() >
            std::numeric_limits<std::uint8_t>::max())
            throw std::runtime_error("too many orbit generators");
        generators_.reserve(matrices.size());
        inverses_.reserve(matrices.size());
        for (const auto& matrix : matrices) {
            generators_.emplace_back(matrix);
            inverses_.emplace_back(inverse(matrix));
        }

        std::vector<std::uint32_t> queue;
        queue.reserve(1u << 23);
        std::vector<OrbitDatum> actual;
        for (std::uint32_t seed = 0; seed < vector_count; ++seed) {
            if (parent_[seed] != no_parent) continue;
            queue.clear();
            queue.push_back(seed);
            parent_[seed] = seed;
            for (std::size_t head = 0; head < queue.size(); ++head) {
                const auto source = queue[head];
                for (std::size_t generator = 0;
                     generator < generators_.size(); ++generator) {
                    const auto target =
                        generators_[generator].image(source);
                    if (parent_[target] == no_parent) {
                        parent_[target] = source;
                        parent_generator_[target] =
                            static_cast<std::uint8_t>(generator);
                        queue.push_back(target);
                    }
                }
            }
            actual.push_back(
                {seed, static_cast<std::uint64_t>(queue.size())});
        }
        if (actual.size() != expected.size())
            throw std::runtime_error("orbit-count mismatch");
        for (std::size_t index = 0; index < expected.size(); ++index)
            if (actual[index].representative !=
                    expected[index].representative ||
                actual[index].size != expected[index].size)
                throw std::runtime_error(
                    "orbit representative/size mismatch");
    }

    std::uint32_t normalize(
        Triple& triple, std::size_t coordinate) const {
        while (parent_[triple[coordinate]] != triple[coordinate]) {
            const auto old = triple[coordinate];
            const auto generator = parent_generator_[old];
            for (auto& component : triple)
                component = inverses_[generator].image(component);
            if (triple[coordinate] != parent_[old])
                throw std::runtime_error(
                    "Schreier parent transport failure");
        }
        return triple[coordinate];
    }

  private:
    std::vector<PackedMatrix> generators_;
    std::vector<PackedMatrix> inverses_;
    std::vector<std::uint32_t> parent_;
    std::vector<std::uint8_t> parent_generator_;
};

struct FirstGroup {
    std::uint32_t a{};
    std::uint64_t first_orbit_size{};
    std::uint64_t stabilizer_order{};
    std::vector<Matrix> generators;
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

std::vector<FirstGroup> read_pair_level(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open pair-level input");
    std::string token;
    input >> token;
    if (token != "CO1_F8_PAIR_LEVEL_V1")
        throw std::runtime_error("bad pair-level header");
    std::vector<FirstGroup> result;
    while (input >> token) {
        if (token != "first")
            throw std::runtime_error("missing first record");
        FirstGroup group;
        std::size_t generator_count = 0;
        input >> group.a;
        input >> token >> group.first_orbit_size;
        if (token != "orbit_size")
            throw std::runtime_error("missing first-orbit size");
        input >> token >> group.stabilizer_order;
        if (token != "stabilizer_order")
            throw std::runtime_error("missing stabilizer order");
        input >> token >> generator_count;
        if (token != "generator_count")
            throw std::runtime_error("missing generator count");
        group.generators.resize(generator_count);
        for (auto& matrix : group.generators) {
            input >> token;
            if (token != "matrix")
                throw std::runtime_error("missing matrix");
            for (auto& row : matrix) {
                input >> row;
                if (row > vector_mask)
                    throw std::runtime_error("matrix row out of range");
            }
            if (direct_image(group.a, matrix) != group.a)
                throw std::runtime_error(
                    "first stabilizer generator does not fix a");
        }
        result.push_back(std::move(group));
    }
    if (result.size() != 4)
        throw std::runtime_error("expected four first-vector groups");
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

struct PairRow {
    std::uint32_t pair_index{};
    std::uint32_t a{};
    std::uint32_t b{};
    std::uint64_t first_orbit_size{};
    std::uint64_t b_orbit_size{};
    std::uint64_t pair_stabilizer_order{};
};

std::vector<PairRow> read_pairs(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open pair rows");
    std::string line;
    std::getline(input, line);
    if (line != "pair_index\ta\tb\tfirst_orbit_size\tb_orbit_size"
                "\tpair_stabilizer_order")
        throw std::runtime_error("bad pair-row header");
    std::vector<PairRow> rows;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tabs(line);
        if (fields.size() != 6)
            throw std::runtime_error("bad pair row");
        PairRow row{
            static_cast<std::uint32_t>(std::stoul(fields[0])),
            static_cast<std::uint32_t>(std::stoul(fields[1])),
            static_cast<std::uint32_t>(std::stoul(fields[2])),
            std::stoull(fields[3]),
            std::stoull(fields[4]),
            std::stoull(fields[5])};
        if (row.pair_index != rows.size() + 1)
            throw std::runtime_error("noncontiguous pair index");
        rows.push_back(row);
    }
    if (rows.size() != 46)
        throw std::runtime_error("expected 46 pair rows");
    return rows;
}

struct Query {
    std::uint64_t query_id{};
    std::string kind;
    std::uint64_t candidate_id{};
    Triple triple{};
    std::uint32_t pair_index{};
};

std::vector<Query> read_queries(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open query file");
    std::string line;
    std::getline(input, line);
    if (line != "query_id\tkind\tcandidate_id\ta\tb\tc")
        throw std::runtime_error("bad query header");
    std::vector<Query> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tabs(line);
        if (fields.size() != 6)
            throw std::runtime_error("bad query row");
        Query query;
        query.query_id = std::stoull(fields[0]);
        query.kind = fields[1];
        query.candidate_id = std::stoull(fields[2]);
        for (int component = 0; component < 3; ++component) {
            const auto value = std::stoul(fields[3 + component]);
            if (value > vector_mask)
                throw std::runtime_error("query component out of range");
            query.triple[component] =
                static_cast<std::uint32_t>(value);
        }
        if (query.query_id != result.size() + 1)
            throw std::runtime_error("noncontiguous query id");
        result.push_back(std::move(query));
    }
    return result;
}

std::uint64_t pair_key(
    std::uint32_t a, std::uint32_t b) noexcept {
    return (static_cast<std::uint64_t>(a) << dimension) | b;
}

void normalize_queries(
    const std::vector<Matrix>& co1_generators,
    const std::vector<FirstGroup>& first_groups,
    const std::vector<PairRow>& pairs,
    std::vector<Query>& queries) {
    std::vector<OrbitDatum> first_expected;
    for (const auto& group : first_groups)
        first_expected.push_back({group.a, group.first_orbit_size});
    OrbitForest first_forest(co1_generators, first_expected);

    std::unordered_map<std::uint32_t, std::size_t> group_of_a;
    for (std::size_t index = 0; index < first_groups.size(); ++index)
        group_of_a.emplace(first_groups[index].a, index);
    std::vector<std::vector<std::size_t>> query_groups(
        first_groups.size());
    for (std::size_t index = 0; index < queries.size(); ++index) {
        const auto a = first_forest.normalize(
            queries[index].triple, 0);
        const auto group = group_of_a.find(a);
        if (group == group_of_a.end())
            throw std::runtime_error("unknown normalized first vector");
        query_groups[group->second].push_back(index);
    }

    std::unordered_map<std::uint64_t, std::uint32_t> pair_of_ab;
    for (const auto& pair : pairs)
        if (!pair_of_ab.emplace(
                pair_key(pair.a, pair.b), pair.pair_index).second)
            throw std::runtime_error("duplicate pair representative");

    for (std::size_t group_index = 0;
         group_index < first_groups.size(); ++group_index) {
        const auto& group = first_groups[group_index];
        std::vector<OrbitDatum> second_expected;
        for (const auto& pair : pairs)
            if (pair.a == group.a)
                second_expected.push_back(
                    {pair.b, pair.b_orbit_size});
        OrbitForest second_forest(
            group.generators, second_expected);
        for (const auto query_index : query_groups[group_index]) {
            auto& query = queries[query_index];
            const auto b = second_forest.normalize(query.triple, 1);
            const auto pair =
                pair_of_ab.find(pair_key(group.a, b));
            if (pair == pair_of_ab.end())
                throw std::runtime_error(
                    "unknown normalized pair representative");
            query.pair_index = pair->second;
        }
        std::cout << "Normalized " << query_groups[group_index].size()
                  << " queries above first coordinate " << group.a << '\n';
    }
}

void write_queries(
    const std::filesystem::path& output_directory,
    const std::vector<Query>& queries) {
    std::filesystem::create_directories(output_directory);
    std::vector<std::ofstream> outputs(47);
    std::vector<std::uint64_t> counts(47, 0);
    for (std::uint32_t pair = 1; pair <= 46; ++pair) {
        const auto path =
            output_directory / ("pair-" + std::to_string(pair) + ".tsv");
        outputs[pair].open(path);
        if (!outputs[pair])
            throw std::runtime_error("cannot open normalized output");
        outputs[pair] << "query_id\tkind\tcandidate_id\tc\n";
    }
    for (const auto& query : queries) {
        if (query.pair_index < 1 || query.pair_index > 46)
            throw std::runtime_error("bad normalized pair index");
        outputs[query.pair_index]
            << query.query_id << '\t'
            << query.kind << '\t'
            << query.candidate_id << '\t'
            << query.triple[2] << '\n';
        ++counts[query.pair_index];
    }
    std::uint64_t total = 0;
    for (std::uint32_t pair = 1; pair <= 46; ++pair) {
        total += counts[pair];
    }
    if (total != queries.size())
        throw std::runtime_error("normalized query-count mismatch");
}

void self_test() {
    Matrix swap{};
    for (std::uint32_t bit = 0; bit < dimension; ++bit)
        swap[bit] = std::uint32_t{1} << bit;
    std::swap(swap[0], swap[1]);
    const auto inv = inverse(swap);
    if (inv != swap ||
        direct_image(5, swap) != 6 ||
        direct_image(6, inv) != 5)
        throw std::runtime_error("self-test failed");
    std::cout << "Normalizer self-test passed\n";
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
                << "usage: normalize_queries ATLAS_INPUT PAIR_LEVEL_INPUT "
                   "PAIR_ROWS QUERY_TSV OUTPUT_DIRECTORY\n";
            return 1;
        }
        const auto co1_generators = read_primal_generators(argv[1]);
        const auto first_groups = read_pair_level(argv[2]);
        const auto pairs = read_pairs(argv[3]);
        auto queries = read_queries(argv[4]);
        normalize_queries(
            co1_generators, first_groups, pairs, queries);
        write_queries(argv[5], queries);
        std::cout << "Normalized queries: " << queries.size() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "NORMALIZE FAILURE: " << error.what() << '\n';
        return 1;
    }
}

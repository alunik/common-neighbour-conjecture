// Classify normalized third coordinates under one pair stabilizer.
//
// The orbit labels are built in memory for one pair at a time.  This avoids
// writing the much larger all-pairs lookup structure to disk.

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
constexpr std::uint32_t no_label =
    std::numeric_limits<std::uint32_t>::max();
using Matrix = std::array<std::uint32_t, dimension>;

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

std::uint32_t rank(Matrix rows) {
    std::uint32_t result = 0;
    for (std::uint32_t column = 0;
         column < dimension && result < dimension; ++column) {
        std::uint32_t pivot = result;
        while (pivot < dimension &&
               ((rows[pivot] >> column) & 1u) == 0)
            ++pivot;
        if (pivot == dimension) continue;
        std::swap(rows[result], rows[pivot]);
        for (std::uint32_t row = 0; row < dimension; ++row)
            if (row != result && ((rows[row] >> column) & 1u))
                rows[row] ^= rows[result];
        ++result;
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

struct PairStabilizer {
    std::uint32_t pair_index{};
    std::uint32_t a{};
    std::uint32_t b{};
    std::uint64_t first_orbit_size{};
    std::uint64_t b_orbit_size{};
    std::uint64_t pair_stabilizer_order{};
    std::vector<Matrix> matrices;
};

PairStabilizer read_stabilizer(const std::string& path) {
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
            throw std::runtime_error("missing matrix");
        for (auto& row : matrix) {
            input >> row;
            if (row > vector_mask)
                throw std::runtime_error("matrix row out of range");
        }
        if (rank(matrix) != dimension)
            throw std::runtime_error("singular matrix");
        if (direct_image(result.a, matrix) != result.a ||
            direct_image(result.b, matrix) != result.b)
            throw std::runtime_error(
                "pair-stabilizer generator does not fix pair");
    }
    if (input >> token)
        throw std::runtime_error("trailing pair-stabilizer input");
    constexpr std::uint64_t co1_order =
        4157776806543360000ULL;
    if (result.first_orbit_size >
            co1_order / result.b_orbit_size ||
        result.first_orbit_size * result.b_orbit_size >
            co1_order / result.pair_stabilizer_order ||
        result.first_orbit_size * result.b_orbit_size *
            result.pair_stabilizer_order != co1_order)
        throw std::runtime_error(
            "pair orbit-stabilizer identity fails");
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

struct TripleRow {
    std::uint32_t pair_index{};
    std::uint32_t a{};
    std::uint32_t b{};
    std::uint32_t c{};
    std::uint64_t first_orbit_size{};
    std::uint64_t b_orbit_size{};
    std::uint64_t c_orbit_size{};
    std::uint64_t triple_stabilizer_order{};
};

std::vector<TripleRow> read_triples(
    const std::string& path, const PairStabilizer& pair) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open triple rows");
    std::string line;
    std::getline(input, line);
    if (line != "pair_index\ta\tb\tc\tfirst_orbit_size"
                "\tb_orbit_size\tc_orbit_size"
                "\ttriple_stabilizer_order")
        throw std::runtime_error("bad triple-row header");
    std::vector<TripleRow> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tabs(line);
        if (fields.size() != 8)
            throw std::runtime_error("bad triple row");
        TripleRow row{
            static_cast<std::uint32_t>(std::stoul(fields[0])),
            static_cast<std::uint32_t>(std::stoul(fields[1])),
            static_cast<std::uint32_t>(std::stoul(fields[2])),
            static_cast<std::uint32_t>(std::stoul(fields[3])),
            std::stoull(fields[4]),
            std::stoull(fields[5]),
            std::stoull(fields[6]),
            std::stoull(fields[7])};
        if (row.pair_index != pair.pair_index ||
            row.a != pair.a || row.b != pair.b ||
            row.first_orbit_size != pair.first_orbit_size ||
            row.b_orbit_size != pair.b_orbit_size)
            throw std::runtime_error("triple/pair metadata mismatch");
        if (row.c_orbit_size == 0 ||
            row.triple_stabilizer_order == 0 ||
            row.c_orbit_size >
                pair.pair_stabilizer_order /
                    row.triple_stabilizer_order ||
            row.c_orbit_size * row.triple_stabilizer_order !=
                pair.pair_stabilizer_order)
            throw std::runtime_error(
                "triple orbit-stabilizer mismatch");
        result.push_back(row);
    }
    if (result.empty())
        throw std::runtime_error("empty triple census");
    return result;
}

struct Query {
    std::uint64_t query_id{};
    std::string kind;
    std::uint64_t candidate_id{};
    std::uint32_t c{};
};

std::vector<Query> read_queries(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open normalized queries");
    std::string line;
    std::getline(input, line);
    if (line != "query_id\tkind\tcandidate_id\tc")
        throw std::runtime_error("bad normalized-query header");
    std::vector<Query> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tabs(line);
        if (fields.size() != 4)
            throw std::runtime_error("bad normalized-query row");
        const auto c = std::stoul(fields[3]);
        if (c > vector_mask)
            throw std::runtime_error("normalized c out of range");
        result.push_back(
            {std::stoull(fields[0]), fields[1],
             std::stoull(fields[2]),
             static_cast<std::uint32_t>(c)});
    }
    return result;
}

struct CensusResult {
    std::vector<std::uint32_t> label;
    std::uint64_t regular_orbits{};
};

CensusResult census(
    const PairStabilizer& pair,
    const std::vector<TripleRow>& triples) {
    std::vector<PackedMatrix> generators;
    generators.reserve(pair.matrices.size());
    for (const auto& matrix : pair.matrices)
        generators.emplace_back(matrix);
    CensusResult result;
    result.label.assign(vector_count, no_label);
    std::vector<std::uint32_t> queue;
    queue.reserve(1u << 20);
    std::size_t orbit_index = 0;
    std::uint64_t mass = 0;

    for (std::uint32_t seed = 0; seed < vector_count; ++seed) {
        if (result.label[seed] != no_label) continue;
        if (orbit_index >= triples.size())
            throw std::runtime_error("more c-orbits than census rows");
        const auto& expected = triples[orbit_index];
        if (expected.c != seed)
            throw std::runtime_error(
                "c-orbit representative mismatch");
        queue.clear();
        queue.push_back(seed);
        result.label[seed] =
            static_cast<std::uint32_t>(orbit_index);
        for (std::size_t head = 0; head < queue.size(); ++head)
            for (const auto& generator : generators) {
                const auto target =
                    generator.image(queue[head]);
                if (result.label[target] == no_label) {
                    result.label[target] =
                        static_cast<std::uint32_t>(orbit_index);
                    queue.push_back(target);
                } else if (result.label[target] != orbit_index) {
                    throw std::runtime_error("c-orbit overlap");
                }
            }
        if (queue.size() != expected.c_orbit_size)
            throw std::runtime_error("c-orbit size mismatch");
        mass += queue.size();
        if (expected.triple_stabilizer_order == 1)
            ++result.regular_orbits;
        ++orbit_index;
    }
    if (orbit_index != triples.size())
        throw std::runtime_error("fewer c-orbits than census rows");
    if (mass != vector_count)
        throw std::runtime_error("c-orbit mass mismatch");
    return result;
}

void classify(
    const PairStabilizer& pair,
    const std::vector<TripleRow>& triples,
    const std::vector<Query>& queries,
    const CensusResult& census_result,
    const std::string& output_path) {
    std::ofstream output(output_path);
    if (!output)
        throw std::runtime_error("cannot open classification output");
    output << "query_id\tkind\tcandidate_id\tpair_index\tc_rep"
              "\ttriple_stabilizer_order\n";
    for (const auto& query : queries) {
        const auto label = census_result.label[query.c];
        if (label == no_label || label >= triples.size())
            throw std::runtime_error("unclassified query");
        const auto& row = triples[label];
        output << query.query_id << '\t'
               << query.kind << '\t'
               << query.candidate_id << '\t'
               << pair.pair_index << '\t'
               << row.c << '\t'
               << row.triple_stabilizer_order << '\n';
    }
}

void self_test() {
    Matrix swap{};
    for (std::uint32_t bit = 0; bit < dimension; ++bit)
        swap[bit] = std::uint32_t{1} << bit;
    std::swap(swap[0], swap[1]);
    PackedMatrix packed(swap);
    if (rank(swap) != dimension ||
        packed.image(5) != 6 ||
        direct_image(5, swap) != 6)
        throw std::runtime_error("self-test failed");
    std::cout << "Classifier self-test passed\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            self_test();
            return 0;
        }
        if (argc != 5) {
            std::cerr
                << "usage: classify_queries PAIR_STABILIZER "
                   "TRIPLE_ROWS NORMALIZED_QUERIES OUTPUT_TSV\n";
            return 1;
        }
        const auto pair = read_stabilizer(argv[1]);
        const auto triples = read_triples(argv[2], pair);
        const auto queries = read_queries(argv[3]);
        const auto census_result = census(pair, triples);
        classify(
            pair, triples, queries, census_result, argv[4]);
        std::cout << "Classified " << queries.size()
                  << " queries above pair " << pair.pair_index << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "PAIR CLASSIFICATION FAILURE: "
                  << error.what() << '\n';
        return 1;
    }
}

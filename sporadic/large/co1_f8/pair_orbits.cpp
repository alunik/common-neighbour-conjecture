// Enumerate Co1-orbits on ordered pairs in the 24-dimensional F_2 module.
// For each first-coordinate orbit, GAP supplies generators for its stabilizer;
// a packed breadth-first search then enumerates the second-coordinate orbits.

#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::uint32_t vector_count = 1u << 24;
using Matrix = std::array<std::uint32_t, 24>;

struct FirstGroup {
    std::uint32_t a{};
    std::uint64_t first_orbit_size{};
    std::uint64_t stabilizer_order{};
    std::vector<Matrix> generators;
};

std::uint32_t image(std::uint32_t vector, const Matrix& matrix) {
    std::uint32_t result = 0;
    while (vector != 0) {
        const unsigned bit = static_cast<unsigned>(__builtin_ctz(vector));
        result ^= matrix[bit];
        vector &= vector - 1;
    }
    return result;
}

std::vector<FirstGroup> read_input(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open pair-level input");
    std::string token;
    input >> token;
    if (token != "CO1_F8_PAIR_LEVEL_V1")
        throw std::runtime_error("bad pair-level header");
    std::vector<FirstGroup> result;
    while (input >> token) {
        if (token != "first")
            throw std::runtime_error("missing first-group record");
        FirstGroup group;
        std::size_t generator_count = 0;
        input >> group.a;
        input >> token >> group.first_orbit_size;
        if (token != "orbit_size")
            throw std::runtime_error("missing first orbit size");
        input >> token >> group.stabilizer_order;
        if (token != "stabilizer_order")
            throw std::runtime_error("missing first stabilizer order");
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
                if (row >= vector_count)
                    throw std::runtime_error("matrix row out of range");
            }
            if (image(group.a, matrix) != group.a)
                throw std::runtime_error("generator does not fix first mask");
        }
        result.push_back(std::move(group));
    }
    if (result.size() != 4)
        throw std::runtime_error("expected four first-vector groups");
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            std::cerr
                << "usage: pair_orbits PAIR_LEVEL_INPUT OUTPUT_TSV\n";
            return 1;
        }
        const auto groups = read_input(argv[1]);
        std::ofstream output(argv[2]);
        if (!output) throw std::runtime_error("cannot open pair output");
        output
            << "pair_index\ta\tb\tfirst_orbit_size"
            << "\tb_orbit_size\tpair_stabilizer_order\n";

        std::vector<std::uint32_t> queue;
        std::vector<std::uint8_t> seen(vector_count, 0);
        std::uint64_t pair_index = 0;
        for (const auto& group : groups) {
            std::fill(seen.begin(), seen.end(), 0);
            std::uint64_t orbit_count = 0;
            std::uint64_t mass = 0;
            for (std::uint32_t seed = 0; seed < vector_count; ++seed) {
                if (seen[seed]) continue;
                queue.clear();
                queue.push_back(seed);
                seen[seed] = 1;
                for (std::size_t head = 0; head < queue.size(); ++head) {
                    for (const auto& generator : group.generators) {
                        const auto target =
                            image(queue[head], generator);
                        if (!seen[target]) {
                            seen[target] = 1;
                            queue.push_back(target);
                        }
                    }
                }
                const auto orbit_size =
                    static_cast<std::uint64_t>(queue.size());
                if (orbit_size == 0 ||
                    group.stabilizer_order % orbit_size != 0)
                    throw std::runtime_error(
                        "nonintegral pair stabilizer order");
                ++pair_index;
                ++orbit_count;
                mass += orbit_size;
                output << pair_index << '\t'
                    << group.a << '\t' << seed << '\t'
                    << group.first_orbit_size << '\t'
                    << orbit_size << '\t'
                    << group.stabilizer_order / orbit_size << '\n';
            }
            if (mass != vector_count)
                throw std::runtime_error("incomplete second-vector census");
            std::cout << "Pair orbits above " << group.a << ": "
                      << orbit_count << '\n';
        }
        std::cout << "Co1 pair orbits on (F_2^24)^2: "
                  << pair_index << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "PAIR-ORBIT FAILURE: " << error.what() << '\n';
        return 1;
    }
}

// Enumerate the four Co1-orbits on the 24-dimensional F_2 module.
//
// The small orbit is also followed in the aligned 98280-point permutation
// representation.  GAP uses this point-to-vector table when computing exact
// stabilizers of pairs of vectors.

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

constexpr std::uint32_t mask24 = (1u << 24) - 1;
constexpr std::size_t degree = 98280;
using Matrix = std::array<std::uint32_t, 24>;
using Permutation = std::vector<std::uint32_t>;

struct AtlasInput {
    std::uint32_t seed_point{};
    std::uint32_t seed_mask{};
    std::array<Matrix, 2> dual_matrices{};
    std::array<Matrix, 2> primal_matrices{};
    std::array<Permutation, 2> permutations{};
};

AtlasInput read_atlas_input(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open Atlas input");
    std::string token;
    in >> token;
    if (token != "CO1_F8_ATLAS_INPUT_V1")
        throw std::runtime_error("bad Atlas input header");

    AtlasInput out;
    in >> token >> out.seed_point;
    if (token != "seed_point" || out.seed_point < 1 || out.seed_point > degree)
        throw std::runtime_error("bad seed point");
    --out.seed_point;
    in >> token >> out.seed_mask;
    if (token != "seed_mask" || out.seed_mask == 0 || out.seed_mask > mask24)
        throw std::runtime_error("bad seed mask");
    for (auto& matrix : out.dual_matrices) {
        in >> token;
        if (token != "dual_matrix") throw std::runtime_error("missing matrix");
        for (auto& row : matrix) {
            in >> row;
            if (row > mask24) throw std::runtime_error("bad matrix row");
        }
    }
    for (auto& matrix : out.primal_matrices) {
        in >> token;
        if (token != "primal_matrix") throw std::runtime_error("missing matrix");
        for (auto& row : matrix) {
            in >> row;
            if (row > mask24) throw std::runtime_error("bad matrix row");
        }
    }
    for (auto& permutation : out.permutations) {
        in >> token;
        if (token != "permutation")
            throw std::runtime_error("missing permutation");
        permutation.resize(degree);
        std::vector<bool> seen(degree, false);
        for (auto& image : permutation) {
            in >> image;
            if (image < 1 || image > degree)
                throw std::runtime_error("bad permutation image");
            --image;
            if (seen[image]) throw std::runtime_error("duplicate image");
            seen[image] = true;
        }
    }
    if (in >> token) throw std::runtime_error("trailing Atlas input");
    return out;
}

std::uint32_t image(std::uint32_t vector, const Matrix& matrix) {
    std::uint32_t result = 0;
    while (vector != 0) {
        const unsigned bit = static_cast<unsigned>(__builtin_ctz(vector));
        result ^= matrix[bit];
        vector &= vector - 1;
    }
    return result;
}

struct VectorOrbitCensus {
    std::vector<std::uint16_t> label;
    std::vector<std::uint32_t> representative;
    std::vector<std::uint32_t> size;
};

VectorOrbitCensus classify_vector_orbits(
    const std::array<Matrix, 2>& generators) {
    constexpr std::uint16_t unassigned = 0xffff;
    constexpr std::uint32_t vector_count = 1u << 24;
    VectorOrbitCensus out;
    out.label.assign(vector_count, unassigned);
    out.label[0] = 0;
    out.representative.push_back(0);
    out.size.push_back(1);

    std::vector<std::uint32_t> queue;
    for (std::uint32_t seed = 1; seed < vector_count; ++seed) {
        if (out.label[seed] != unassigned) continue;
        if (out.representative.size() == unassigned)
            throw std::runtime_error("too many vector orbits");
        const auto orbit_label =
            static_cast<std::uint16_t>(out.representative.size());
        out.representative.push_back(seed);
        queue.clear();
        queue.push_back(seed);
        out.label[seed] = orbit_label;
        for (std::size_t head = 0; head < queue.size(); ++head) {
            for (const auto& generator : generators) {
                const auto target = image(queue[head], generator);
                if (out.label[target] == unassigned) {
                    out.label[target] = orbit_label;
                    queue.push_back(target);
                } else if (out.label[target] != orbit_label) {
                    throw std::runtime_error("vector orbit overlap");
                }
            }
        }
        out.size.push_back(static_cast<std::uint32_t>(queue.size()));
    }
    std::uint64_t total = 0;
    for (const auto size : out.size) total += size;
    if (total != vector_count)
        throw std::runtime_error("vector orbit census is incomplete");
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: vector_orbits GROUP_DATA OUTPUT_PREFIX\n";
        return 1;
    }
    const std::string input_path = argv[1];
    const std::string prefix = argv[2];
    const auto input = read_atlas_input(input_path);

    std::vector<std::uint32_t> orbit(degree, 0);
    std::vector<bool> assigned(degree, false);
    std::unordered_map<std::uint32_t, std::uint32_t> point_of_vector;
    std::queue<std::uint32_t> pending;
    orbit[input.seed_point] = input.seed_mask;
    assigned[input.seed_point] = true;
    point_of_vector.emplace(input.seed_mask, input.seed_point);
    pending.push(input.seed_point);

    std::size_t assigned_count = 1;
    while (!pending.empty()) {
        const auto point = pending.front();
        pending.pop();
        for (int generator = 0; generator < 2; ++generator) {
            const auto target_point = input.permutations[generator][point];
            const auto target_vector =
                image(orbit[point], input.dual_matrices[generator]);
            if (target_vector == 0)
                throw std::runtime_error("singular generator");
            if (!assigned[target_point]) {
                if (!point_of_vector.emplace(target_vector, target_point).second)
                    throw std::runtime_error("two points have the same line");
                orbit[target_point] = target_vector;
                assigned[target_point] = true;
                ++assigned_count;
                pending.push(target_point);
            } else if (orbit[target_point] != target_vector) {
                throw std::runtime_error("matrix/permutation alignment failure");
            }
        }
    }
    if (assigned_count != degree || point_of_vector.size() != degree)
        throw std::runtime_error("dual orbit is incomplete");

    const auto dual_census =
        classify_vector_orbits(input.dual_matrices);
    const auto vector_census =
        classify_vector_orbits(input.primal_matrices);
    if (dual_census.size.size() != 4)
        throw std::runtime_error("unexpected number of dual-vector orbits");
    if (vector_census.size.size() != 4)
        throw std::runtime_error("unexpected number of vector orbits");
    std::array<std::uint16_t, 3> dual_orbit_labels{{1,2,3}};
    std::sort(dual_orbit_labels.begin(), dual_orbit_labels.end(),
        [&](std::uint16_t a, std::uint16_t b) {
            return dual_census.size[a] < dual_census.size[b];
        });
    std::array<std::uint32_t, 3> dual_orbit_sizes{};
    for (int i = 0; i < 3; ++i)
        dual_orbit_sizes[i] = dual_census.size[dual_orbit_labels[i]];
    if (dual_orbit_sizes !=
        std::array<std::uint32_t, 3>{{98280,8292375,8386560}})
        throw std::runtime_error("unexpected dual-vector orbit sizes");
    const auto small_orbit_label = dual_orbit_labels[0];
    for (const auto functional : orbit)
        if (dual_census.label[functional] != small_orbit_label)
            throw std::runtime_error(
                "degree-98280 orbit does not match dual census");

    {
        std::ofstream out(prefix + "_orbit.tsv");
        out << "point\tfunctional_mask\n";
        for (std::size_t point = 0; point < degree; ++point)
            out << (point + 1) << '\t' << orbit[point] << '\n';
    }
    {
        std::ofstream out(prefix + "_dual_vector_orbits.tsv");
        out << "orbit\trepresentative_mask\tsize\n";
        for (std::size_t i = 0; i < dual_census.size.size(); ++i)
            out << i << '\t' << dual_census.representative[i]
                << '\t' << dual_census.size[i] << '\n';
    }
    {
        std::ofstream out(prefix + "_vector_orbits.tsv");
        out << "orbit\trepresentative_mask\tsize\n";
        for (std::size_t i = 0; i < vector_census.size.size(); ++i)
            out << i << '\t' << vector_census.representative[i]
                << '\t' << vector_census.size[i] << '\n';
    }

    std::cout << "Co1 vector orbits on F_2^24: "
              << vector_census.size.size() << " (sizes";
    for (const auto size : vector_census.size) std::cout << ' ' << size;
    std::cout << ")\n";
    return 0;
}

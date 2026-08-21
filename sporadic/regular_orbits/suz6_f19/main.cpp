#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Construct the 32760-point orbit of projective functionals used to test
// regular vectors for the two dual 12-dimensional F_19 modules of 6.Suz.
// The output is deliberately transient: run.sh puts it in a temporary
// directory for GAP to consume and removes it at the end of the run.

namespace {

constexpr int field_order = 19;
constexpr int dimension = 12;
using Vector = std::array<std::uint8_t, dimension>;
using Matrix = std::array<Vector, dimension>;

// Standard generators of 6.Suz on its 12-dimensional F_19 module.
constexpr Matrix generator_a{{
    {{0,1,0,0,0,0,0,0,0,0,0,0}},
    {{18,0,0,0,0,0,0,0,0,0,0,0}},
    {{0,0,0,1,0,0,0,0,0,0,0,0}},
    {{0,0,18,0,0,0,0,0,0,0,0,0}},
    {{0,0,0,0,0,0,1,0,0,0,0,0}},
    {{0,0,0,0,0,0,0,1,0,0,0,0}},
    {{0,0,0,0,18,0,0,0,0,0,0,0}},
    {{0,0,0,0,0,18,0,0,0,0,0,0}},
    {{18,16,0,3,17,10,15,7,0,0,12,3}},
    {{11,14,0,5,3,4,6,18,0,0,1,13}},
    {{14,8,5,0,6,18,16,15,6,3,0,0}},
    {{0,0,0,0,0,0,0,0,1,7,0,0}}
}};

constexpr Matrix generator_b{{
    {{11,0,0,0,0,0,0,0,0,0,0,0}},
    {{0,0,1,0,0,0,0,0,0,0,0,0}},
    {{0,0,0,0,1,0,0,0,0,0,0,0}},
    {{0,0,0,0,0,1,0,0,0,0,0,0}},
    {{0,1,0,0,0,0,0,0,0,0,0,0}},
    {{5,11,14,0,13,1,3,4,13,16,0,0}},
    {{0,0,0,0,0,0,0,0,0,1,0,0}},
    {{0,0,0,0,0,0,0,0,11,0,0,0}},
    {{12,8,13,9,17,10,12,5,3,7,0,0}},
    {{15,11,14,7,13,12,15,4,13,5,0,0}},
    {{7,0,0,0,0,0,0,0,0,0,1,0}},
    {{9,14,15,10,9,9,11,8,7,8,0,1}}
}};

int mod(int value) {
    value %= field_order;
    return value < 0 ? value + field_order : value;
}

int inverse_scalar(int value) {
    if (value == 0) {
        throw std::runtime_error("division by zero");
    }
    for (int candidate = 1; candidate < field_order; ++candidate) {
        if ((value * candidate) % field_order == 1) {
            return candidate;
        }
    }
    throw std::runtime_error("nonzero field element has no inverse");
}

Matrix identity_matrix() {
    Matrix result{};
    for (int i = 0; i < dimension; ++i) {
        result[i][i] = 1;
    }
    return result;
}

Matrix inverse(Matrix matrix) {
    Matrix result = identity_matrix();
    for (int column = 0; column < dimension; ++column) {
        int pivot = column;
        while (pivot < dimension && matrix[pivot][column] == 0) {
            ++pivot;
        }
        if (pivot == dimension) {
            throw std::runtime_error("singular generator matrix");
        }
        if (pivot != column) {
            std::swap(matrix[pivot], matrix[column]);
            std::swap(result[pivot], result[column]);
        }

        const int scale = inverse_scalar(matrix[column][column]);
        for (int j = 0; j < dimension; ++j) {
            matrix[column][j] = mod(matrix[column][j] * scale);
            result[column][j] = mod(result[column][j] * scale);
        }
        for (int row = 0; row < dimension; ++row) {
            if (row == column || matrix[row][column] == 0) {
                continue;
            }
            const int scale_row = matrix[row][column];
            for (int j = 0; j < dimension; ++j) {
                matrix[row][j] =
                    mod(matrix[row][j] - scale_row * matrix[column][j]);
                result[row][j] =
                    mod(result[row][j] - scale_row * result[column][j]);
            }
        }
    }
    return result;
}

Matrix transpose(const Matrix& matrix) {
    Matrix result{};
    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            result[i][j] = matrix[j][i];
        }
    }
    return result;
}

Vector multiply(const Vector& vector, const Matrix& matrix) {
    Vector result{};
    for (int j = 0; j < dimension; ++j) {
        int value = 0;
        for (int i = 0; i < dimension; ++i) {
            value += vector[i] * matrix[i][j];
        }
        result[j] = static_cast<std::uint8_t>(value % field_order);
    }
    return result;
}

Vector normalize(Vector vector) {
    int first = 0;
    while (first < dimension && vector[first] == 0) {
        ++first;
    }
    if (first == dimension) {
        throw std::runtime_error("cannot normalize the zero vector");
    }
    const int scale = inverse_scalar(vector[first]);
    for (auto& entry : vector) {
        entry = static_cast<std::uint8_t>((entry * scale) % field_order);
    }
    return vector;
}

std::uint64_t encode(const Vector& vector) {
    std::uint64_t code = 0;
    for (int i = dimension - 1; i >= 0; --i) {
        code = code * field_order + vector[i];
    }
    return code;
}

void write_vector(std::ostream& output, const Vector& vector) {
    for (int i = 0; i < dimension; ++i) {
        if (i != 0) {
            output << ',';
        }
        output << static_cast<int>(vector[i]);
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: " << argv[0] << " OUTPUT_PREFIX A|B\n";
        return 1;
    }

    const std::string output_prefix = argv[1];
    const std::string module = argv[2];

    // Functionals on module A transform by g^(-T).  The second module is
    // contragredient, so its functionals transform by g.
    std::array<Matrix, 2> functional_generators{};
    Vector seed{};
    if (module == "A") {
        functional_generators = {
            transpose(inverse(generator_a)),
            transpose(inverse(generator_b))};
        seed = normalize(Vector{{1,18,0,1,8,0,7,11,18,12,0,0}});
    } else if (module == "B") {
        functional_generators = {generator_a, generator_b};
        seed = normalize(Vector{{1,4,14,5,8,12,0,16,2,5,14,14}});
    } else {
        std::cerr << "module must be A or B\n";
        return 1;
    }

    std::vector<Vector> orbit{seed};
    std::unordered_map<std::uint64_t, std::uint32_t> index;
    std::queue<std::uint32_t> pending;
    index.emplace(encode(seed), 0);
    pending.push(0);

    while (!pending.empty()) {
        const std::uint32_t point = pending.front();
        pending.pop();
        for (const Matrix& generator : functional_generators) {
            const Vector image = normalize(multiply(orbit[point], generator));
            const std::uint64_t key = encode(image);
            if (index.emplace(key, static_cast<std::uint32_t>(orbit.size())).second) {
                orbit.push_back(image);
                pending.push(static_cast<std::uint32_t>(orbit.size() - 1));
            }
        }
    }
    if (orbit.size() != 32760) {
        std::cerr << "unexpected projective orbit size: " << orbit.size() << '\n';
        return 2;
    }

    std::array<std::vector<std::uint32_t>, 2> permutations;
    for (auto& permutation : permutations) {
        permutation.resize(orbit.size());
    }
    for (std::uint32_t point = 0; point < orbit.size(); ++point) {
        for (int k = 0; k < 2; ++k) {
            const Vector image =
                normalize(multiply(orbit[point], functional_generators[k]));
            const auto position = index.find(encode(image));
            if (position == index.end()) {
                throw std::runtime_error("the projective orbit is not closed");
            }
            permutations[k][point] = position->second;
        }
    }

    std::ofstream orbit_file(output_prefix + "_orbit.tsv");
    if (!orbit_file) {
        throw std::runtime_error("cannot open orbit output");
    }
    orbit_file << "index\tfunctional\n";
    for (std::uint32_t point = 0; point < orbit.size(); ++point) {
        orbit_file << (point + 1) << '\t';
        write_vector(orbit_file, orbit[point]);
        orbit_file << '\n';
    }

    std::ofstream permutation_file(output_prefix + "_perms.g");
    if (!permutation_file) {
        throw std::runtime_error("cannot open permutation output");
    }
    for (int k = 0; k < 2; ++k) {
        permutation_file << "SUZ6_F19_PERM_" << (k == 0 ? 'A' : 'B')
                         << " := PermList([";
        for (std::uint32_t point = 0; point < orbit.size(); ++point) {
            if (point != 0) {
                permutation_file << ',';
            }
            permutation_file << (permutations[k][point] + 1);
        }
        permutation_file << "]);;\n";
    }

    std::cout << "Module " << module
              << ": constructed a projective orbit of " << orbit.size()
              << " functionals.\n";
    return 0;
}

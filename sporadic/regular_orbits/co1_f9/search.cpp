#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iterator>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Candidate search for 2.Co1 on GF(9)^24.
//
// We work over the GF(3)-form of the module.  If v is the vector fixed by
// an M24 complement and omega is outside GF(3), a candidate has the form
//
//                         x = v + omega*w.
//
// The orbit of a dual vector has length 196560.  Counting the members of
// this orbit that vanish on each of the four lines in <v,w> gives a quick
// 2.Co1-invariant test.  We retain w when the line <v> is distinguished
// from the other three lines.  GAP then does the exact stabiliser and orbit
// tests; this program only supplies a deterministic stream of candidates.

namespace {

constexpr int field_size = 3;
constexpr int dimension = 24;
constexpr int block_size = 6;
constexpr int block_count = dimension / block_size;
constexpr int block_states = 729;  // 3^6
constexpr std::size_t dual_orbit_size = 196560;

using Vector = std::array<std::uint8_t, dimension>;
using Matrix = std::array<Vector, dimension>;

struct Representation {
    Vector fixed_vector{};
    Vector dual_seed{};
    std::array<Matrix, 2> generators{};
};

std::vector<std::string> split(const std::string& text, char separator) {
    std::vector<std::string> fields;
    std::stringstream input(text);
    std::string field;
    while (std::getline(input, field, separator)) fields.push_back(field);
    return fields;
}

Vector parse_vector(const std::string& text) {
    const auto fields = split(text, ',');
    if (fields.size() != dimension)
        throw std::runtime_error("a vector has the wrong length");

    Vector vector{};
    for (int i = 0; i < dimension; ++i) {
        const int value = std::stoi(fields[i]);
        if (value < 0 || value >= field_size)
            throw std::runtime_error("an entry is not in GF(3)");
        vector[i] = static_cast<std::uint8_t>(value);
    }
    return vector;
}

Matrix parse_matrix(const std::string& text) {
    const auto fields = split(text, ',');
    if (fields.size() != dimension * dimension)
        throw std::runtime_error("a matrix has the wrong size");

    Matrix matrix{};
    for (int row = 0; row < dimension; ++row) {
        for (int column = 0; column < dimension; ++column) {
            const int value = std::stoi(fields[dimension * row + column]);
            if (value < 0 || value >= field_size)
                throw std::runtime_error("an entry is not in GF(3)");
            matrix[row][column] = static_cast<std::uint8_t>(value);
        }
    }
    return matrix;
}

Representation read_representation(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open representation.txt");

    const std::string raw((std::istreambuf_iterator<char>(input)),
                          std::istreambuf_iterator<char>());
    std::string text;
    text.reserve(raw.size());
    for (std::size_t i = 0; i < raw.size(); ++i) {
        if (raw[i] == '\\' && i + 1 < raw.size() && raw[i + 1] == '\n') {
            ++i;
        } else if (raw[i] != '\r') {
            text.push_back(raw[i]);
        }
    }

    std::stringstream lines(text);
    std::string line;
    std::getline(lines, line);
    if (line != "CO1_F9_REPRESENTATION")
        throw std::runtime_error("unexpected representation header");

    Representation representation;
    int entries_seen = 0;
    while (std::getline(lines, line)) {
        const auto space = line.find(' ');
        if (space == std::string::npos) continue;
        const std::string name = line.substr(0, space);
        const std::string value = line.substr(space + 1);
        if (name == "fixed_vector") {
            representation.fixed_vector = parse_vector(value);
            entries_seen |= 1;
        } else if (name == "dual_seed") {
            representation.dual_seed = parse_vector(value);
            entries_seen |= 2;
        } else if (name == "generator_1") {
            representation.generators[0] = parse_matrix(value);
            entries_seen |= 4;
        } else if (name == "generator_2") {
            representation.generators[1] = parse_matrix(value);
            entries_seen |= 8;
        }
    }
    if (entries_seen != 15)
        throw std::runtime_error("representation.txt is incomplete");
    return representation;
}

Matrix identity_matrix() {
    Matrix matrix{};
    for (int i = 0; i < dimension; ++i) matrix[i][i] = 1;
    return matrix;
}

int inverse_scalar(int value) {
    if (value == 1) return 1;
    if (value == 2) return 2;
    throw std::runtime_error("division by zero");
}

Matrix inverse(Matrix matrix) {
    Matrix result = identity_matrix();
    for (int column = 0; column < dimension; ++column) {
        int pivot = column;
        while (pivot < dimension && matrix[pivot][column] == 0) ++pivot;
        if (pivot == dimension) throw std::runtime_error("singular generator");
        if (pivot != column) {
            std::swap(matrix[pivot], matrix[column]);
            std::swap(result[pivot], result[column]);
        }

        const int scale = inverse_scalar(matrix[column][column]);
        for (int j = 0; j < dimension; ++j) {
            matrix[column][j] =
                static_cast<std::uint8_t>(matrix[column][j] * scale % 3);
            result[column][j] =
                static_cast<std::uint8_t>(result[column][j] * scale % 3);
        }
        for (int row = 0; row < dimension; ++row) {
            if (row == column || matrix[row][column] == 0) continue;
            const int coefficient = matrix[row][column];
            for (int j = 0; j < dimension; ++j) {
                matrix[row][j] = static_cast<std::uint8_t>(
                    (matrix[row][j] + 3 -
                     coefficient * matrix[column][j] % 3) % 3);
                result[row][j] = static_cast<std::uint8_t>(
                    (result[row][j] + 3 -
                     coefficient * result[column][j] % 3) % 3);
            }
        }
    }
    return result;
}

Matrix transpose(const Matrix& matrix) {
    Matrix result{};
    for (int row = 0; row < dimension; ++row)
        for (int column = 0; column < dimension; ++column)
            result[row][column] = matrix[column][row];
    return result;
}

Vector multiply(const Vector& vector, const Matrix& matrix) {
    Vector result{};
    for (int column = 0; column < dimension; ++column) {
        unsigned value = 0;
        for (int row = 0; row < dimension; ++row)
            value += vector[row] * matrix[row][column];
        result[column] = static_cast<std::uint8_t>(value % 3);
    }
    return result;
}

int dot(const Vector& left, const Vector& right) {
    unsigned value = 0;
    for (int i = 0; i < dimension; ++i) value += left[i] * right[i];
    return static_cast<int>(value % 3);
}

std::uint64_t encode(const Vector& vector) {
    std::uint64_t value = 0;
    for (int i = dimension - 1; i >= 0; --i)
        value = 3 * value + vector[i];
    return value;
}

std::array<std::uint16_t, block_count> blocks(const Vector& vector) {
    std::array<std::uint16_t, block_count> result{};
    for (int block = 0; block < block_count; ++block) {
        int value = 0;
        for (int i = block_size - 1; i >= 0; --i)
            value = 3 * value + vector[block * block_size + i];
        result[block] = static_cast<std::uint16_t>(value);
    }
    return result;
}

// A fixed seed makes the candidate stream identical on every machine.
std::uint64_t next_random(std::uint64_t& state) {
    std::uint64_t value = (state += 0x9e3779b97f4a7c15ULL);
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

void print_vector(const Vector& vector) {
    for (int i = 0; i < dimension; ++i) {
        if (i != 0) std::cout << ',';
        std::cout << static_cast<int>(vector[i]);
    }
    std::cout << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            std::cerr << "usage: search REPRESENTATION NUMBER_OF_CANDIDATES\n";
            return 1;
        }
        const Representation representation = read_representation(argv[1]);
        const int requested = std::stoi(argv[2]);
        if (requested < 1) throw std::runtime_error("empty candidate search");

        // Row vectors in the dual module act by inverse transposes.
        const std::array<Matrix, 2> dual_generators{
            transpose(inverse(representation.generators[0])),
            transpose(inverse(representation.generators[1]))};

        std::vector<Vector> orbit;
        orbit.reserve(dual_orbit_size);
        std::unordered_map<std::uint64_t, std::uint32_t> position;
        position.reserve(2 * dual_orbit_size);
        std::queue<std::uint32_t> pending;
        orbit.push_back(representation.dual_seed);
        position.emplace(encode(representation.dual_seed), 0);
        pending.push(0);

        while (!pending.empty()) {
            const auto current = pending.front();
            pending.pop();
            for (const Matrix& generator : dual_generators) {
                const Vector image = multiply(orbit[current], generator);
                const auto inserted = position.emplace(
                    encode(image), static_cast<std::uint32_t>(orbit.size()));
                if (inserted.second) {
                    orbit.push_back(image);
                    if (orbit.size() > dual_orbit_size)
                        throw std::runtime_error("dual orbit is too large");
                    pending.push(static_cast<std::uint32_t>(orbit.size() - 1));
                }
            }
        }
        if (orbit.size() != dual_orbit_size)
            throw std::runtime_error("dual orbit has the wrong size");

        std::vector<std::uint8_t> fixed_values(orbit.size());
        std::vector<std::array<std::uint16_t, block_count>> orbit_blocks(
            orbit.size());
        for (std::size_t i = 0; i < orbit.size(); ++i) {
            fixed_values[i] = static_cast<std::uint8_t>(
                dot(orbit[i], representation.fixed_vector));
            orbit_blocks[i] = blocks(orbit[i]);
        }
        const std::uint32_t fixed_zero_count =
            static_cast<std::uint32_t>(std::count(
                fixed_values.begin(), fixed_values.end(), 0));
        if (fixed_zero_count != 65826)
            throw std::runtime_error("unexpected zero count on the fixed line");

        std::cout << "vector\n";
        std::unordered_set<std::uint64_t> seen_candidates;
        seen_candidates.reserve(2 * static_cast<std::size_t>(requested));
        std::uint64_t state = 0x5f3759df6b8b4567ULL;
        int accepted = 0;

        while (accepted < requested) {
            Vector candidate{};
            do {
                for (auto& entry : candidate)
                    entry = static_cast<std::uint8_t>(next_random(state) % 3);
            } while (encode(candidate) == 0);
            if (!seen_candidates.insert(encode(candidate)).second) continue;

            // Six-coordinate lookup tables make the 196560 dot products
            // inexpensive.  The table entry is ell(w) on one block.
            std::array<std::array<std::uint8_t, block_states>, block_count>
                dot_table{};
            for (int block = 0; block < block_count; ++block) {
                for (int code = 0; code < block_states; ++code) {
                    int remainder = code;
                    int value = 0;
                    for (int i = 0; i < block_size; ++i) {
                        value += (remainder % 3) *
                                 candidate[block * block_size + i];
                        remainder /= 3;
                    }
                    dot_table[block][code] =
                        static_cast<std::uint8_t>(value % 3);
                }
            }

            std::array<std::uint32_t, 9> counts{};
            for (std::size_t i = 0; i < orbit.size(); ++i) {
                int value = 0;
                for (int block = 0; block < block_count; ++block)
                    value += dot_table[block][orbit_blocks[i][block]];
                const int fixed_value = fixed_values[i];
                const int candidate_value = value % 3;
                ++counts[fixed_value + 3 * candidate_value];
            }

            // The four quantities count zeros on <w>, <v>, <v-w>, <v+w>.
            const std::array<std::uint32_t, 4> non_joint{
                counts[1] + counts[2],
                counts[3] + counts[6],
                counts[4] + counts[8],
                counts[5] + counts[7]};
            if (counts[0] + non_joint[1] != fixed_zero_count)
                throw std::runtime_error("the fixed-line count changed");
            if (non_joint[1] == non_joint[0] ||
                non_joint[1] == non_joint[2] ||
                non_joint[1] == non_joint[3])
                continue;

            ++accepted;
            print_vector(candidate);
            if (accepted % 1000 == 0)
                std::cerr << "candidates " << accepted << '/' << requested
                          << '\n';
        }

        std::cerr << "dual_orbit_size " << orbit.size() << '\n'
                  << "fixed_line_zero_count " << fixed_zero_count << '\n'
                  << "candidates_written " << accepted << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "co1_f9: " << error.what() << '\n';
        return 2;
    }
}

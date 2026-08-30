// Exact product-action stabilizer-orbit classifier.
//
// The point stabilizer contains K^k, where K is the component point
// stabilizer.  Its orbits on Delta^k are Cartesian products of K-orbits.
// The quotient generators act on this compact type space.  Thus H-orbits
// can be enumerated exactly without materialising degree-sized permutations.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;

namespace {

using Permutation = std::vector<u32>;

struct QuotientGenerator {
    Permutation top;
    std::vector<Permutation> components;
};

struct Action {
    std::string label;
    u64 stabilizer_order = 0;
    u64 quotient_order = 0;
    u64 check_index = 0;
    std::vector<QuotientGenerator> quotient_generators;
};

struct Layer {
    u32 component_degree = 0;
    u32 exponent = 0;
    u64 socle_order = 0;
    u64 outer_order = 0;
    std::vector<Permutation> socle_generators;
    std::vector<Permutation> stabilizer_generators;
    std::vector<Action> actions;
};

u64 checkedInteger(const std::string& text, u64 maximum,
                   const std::string& description, bool positive = true) {
    if (text.empty()) throw std::runtime_error("missing " + description);
    u64 value = 0;
    for (unsigned char c : text) {
        if (c < '0' || c > '9') throw std::runtime_error("bad " + description);
        const u64 digit = static_cast<u64>(c - '0');
        if (value > (maximum - digit) / 10)
            throw std::runtime_error(description + " is too large");
        value = 10 * value + digit;
    }
    if (positive && value == 0) throw std::runtime_error("bad " + description);
    return value;
}

void expect(std::istream& input, const std::string& expected) {
    std::string token;
    if (!(input >> token) || token != expected)
        throw std::runtime_error("expected '" + expected + "'");
}

void validatePermutation(const Permutation& permutation, u32 degree,
                         const std::string& description) {
    if (permutation.size() != degree)
        throw std::runtime_error(description + " has wrong degree");
    std::vector<unsigned char> seen(degree, 0);
    for (u32 image : permutation) {
        if (image >= degree || seen[image])
            throw std::runtime_error(description + " is not a permutation");
        seen[image] = 1;
    }
}

Permutation readPackedPermutation(std::istream& input, u32 degree,
                                  const std::string& prefix) {
    expect(input, prefix);
    u32 width = 0;
    std::string packed;
    if (!(input >> width >> packed) || width < 1 || width > 5 ||
        packed.size() != static_cast<std::size_t>(degree) * width)
        throw std::runtime_error("bad packed " + prefix);
    Permutation permutation(degree);
    for (u32 point = 0; point < degree; ++point) {
        u64 image = 0, place = 1;
        for (u32 digit = 0; digit < width; ++digit) {
            const unsigned char code = static_cast<unsigned char>(
                packed[static_cast<std::size_t>(point) * width + digit]);
            if (code < 33 || code > 126)
                throw std::runtime_error("bad packed permutation digit");
            image += static_cast<u64>(code - 33) * place;
            place *= 94;
        }
        if (image >= degree)
            throw std::runtime_error("packed permutation image is out of range");
        permutation[point] = static_cast<u32>(image);
    }
    validatePermutation(permutation, degree, prefix);
    return permutation;
}

Permutation readTopPermutation(std::istream& input, u32 degree) {
    expect(input, "top_gen");
    Permutation permutation(degree);
    for (u32 point = 0; point < degree; ++point) {
        std::string image;
        if (!(input >> image)) throw std::runtime_error("missing top image");
        permutation[point] = static_cast<u32>(
            checkedInteger(image, degree, "top image") - 1);
    }
    validatePermutation(permutation, degree, "top generator");
    return permutation;
}

Layer readLayer(std::istream& input) {
    expect(input, "PRIMITIVE_SAXL_PA_STRUCTURED_V1");
    expect(input, "layer");
    Layer layer;
    std::string value;
    expect(input, "component_degree"); input >> value;
    layer.component_degree = static_cast<u32>(checkedInteger(
        value, std::numeric_limits<u32>::max(), "component degree"));
    expect(input, "exponent"); input >> value;
    layer.exponent = static_cast<u32>(checkedInteger(value, 32, "exponent"));
    expect(input, "socle_order"); input >> value;
    layer.socle_order = checkedInteger(value, std::numeric_limits<u64>::max(), "socle order");
    expect(input, "outer_order"); input >> value;
    layer.outer_order = checkedInteger(value, std::numeric_limits<u64>::max(), "outer order");

    u64 count = 0;
    expect(input, "socle_gens");
    if (!(input >> count) || count == 0 || count > 100)
        throw std::runtime_error("bad socle generator count");
    for (u64 index = 0; index < count; ++index)
        layer.socle_generators.push_back(readPackedPermutation(
            input, layer.component_degree, "component_gen"));
    expect(input, "socle_stabilizer_gens");
    if (!(input >> count) || count > 100)
        throw std::runtime_error("bad stabilizer generator count");
    for (u64 index = 0; index < count; ++index)
        layer.stabilizer_generators.push_back(readPackedPermutation(
            input, layer.component_degree, "component_gen"));

    u64 action_count = 0;
    expect(input, "actions");
    if (!(input >> action_count) || action_count == 0 || action_count > 1000000)
        throw std::runtime_error("bad action count");
    layer.actions.reserve(static_cast<std::size_t>(action_count));
    for (u64 action_index = 0; action_index < action_count; ++action_index) {
        expect(input, "action");
        Action action;
        expect(input, "label"); input >> action.label;
        expect(input, "stabilizer_order"); input >> value;
        action.stabilizer_order = checkedInteger(
            value, std::numeric_limits<u64>::max(), "stabilizer order");
        expect(input, "quotient_order"); input >> value;
        action.quotient_order = checkedInteger(
            value, std::numeric_limits<u64>::max(), "quotient order");
        expect(input, "smallq_check_index");
        if (!(input >> action.check_index) || action.check_index == 0)
            throw std::runtime_error("bad check index");
        u64 generator_count = 0;
        expect(input, "quotient_gens");
        if (!(input >> generator_count) || generator_count == 0 || generator_count > 100)
            throw std::runtime_error("bad quotient generator count");
        for (u64 generator_index = 0; generator_index < generator_count; ++generator_index) {
            expect(input, "qgen");
            QuotientGenerator generator;
            generator.top = readTopPermutation(input, layer.exponent);
            u64 component_count = 0;
            expect(input, "component_gens");
            if (!(input >> component_count) || component_count != layer.exponent)
                throw std::runtime_error("bad quotient component count");
            for (u64 coordinate = 0; coordinate < component_count; ++coordinate)
                generator.components.push_back(readPackedPermutation(
                    input, layer.component_degree, "component_gen"));
            expect(input, "end_qgen");
            action.quotient_generators.push_back(std::move(generator));
        }
        expect(input, "end");
        layer.actions.push_back(std::move(action));
    }
    expect(input, "end_layer");
    std::string trailing;
    if (input >> trailing) throw std::runtime_error("trailing structured input");
    return layer;
}

Permutation inverse(const Permutation& permutation) {
    Permutation result(permutation.size());
    for (u32 point = 0; point < permutation.size(); ++point)
        result[permutation[point]] = point;
    return result;
}

Permutation compose(const Permutation& left, const Permutation& right) {
    if (left.size() != right.size()) throw std::runtime_error("composition degree mismatch");
    Permutation result(left.size());
    for (u32 point = 0; point < left.size(); ++point)
        result[point] = right[left[point]];
    return result;
}

std::vector<Permutation> componentTransportersToZero(const Layer& layer) {
    const u32 n = layer.component_degree;
    std::vector<Permutation> moves = layer.socle_generators;
    for (const auto& generator : layer.socle_generators)
        moves.push_back(inverse(generator));
    std::vector<u32> parent(n, std::numeric_limits<u32>::max());
    std::vector<u32> parent_move(n, std::numeric_limits<u32>::max());
    std::queue<u32> queue;
    parent[0] = 0; queue.push(0);
    while (!queue.empty()) {
        const u32 point = queue.front(); queue.pop();
        for (u32 move = 0; move < moves.size(); ++move) {
            const u32 image = moves[move][point];
            if (parent[image] != std::numeric_limits<u32>::max()) continue;
            parent[image] = point; parent_move[image] = move; queue.push(image);
        }
    }
    if (std::find(parent.begin(), parent.end(), std::numeric_limits<u32>::max()) != parent.end())
        throw std::runtime_error("component socle is not transitive");
    std::vector<Permutation> inverse_moves;
    for (const auto& move : moves) inverse_moves.push_back(inverse(move));
    Permutation identity(n);
    for (u32 point = 0; point < n; ++point) identity[point] = point;
    std::vector<Permutation> transporters(n, Permutation(n));
    for (u32 target = 0; target < n; ++target) {
        Permutation transporter = identity;
        u32 point = target;
        while (point != 0) {
            const u32 move = parent_move[point];
            transporter = compose(transporter, inverse_moves[move]);
            point = parent[point];
        }
        if (transporter[target] != 0)
            throw std::runtime_error("component transporter certification failed");
        transporters[target] = std::move(transporter);
    }
    return transporters;
}

struct ComponentOrbits {
    std::vector<u32> index;
    std::vector<u32> sizes;
    std::vector<u32> representatives;
};

ComponentOrbits componentOrbits(const Layer& layer) {
    const u32 n = layer.component_degree;
    ComponentOrbits result;
    result.index.assign(n, std::numeric_limits<u32>::max());
    std::vector<u32> queue;
    for (u32 seed = 0; seed < n; ++seed) {
        if (result.index[seed] != std::numeric_limits<u32>::max()) continue;
        const u32 orbit = static_cast<u32>(result.sizes.size());
        queue.clear(); queue.push_back(seed); result.index[seed] = orbit;
        for (std::size_t head = 0; head < queue.size(); ++head) {
            for (const auto& generator : layer.stabilizer_generators) {
                const u32 image = generator[queue[head]];
                if (result.index[image] != std::numeric_limits<u32>::max()) continue;
                result.index[image] = orbit; queue.push_back(image);
            }
        }
        result.representatives.push_back(seed);
        result.sizes.push_back(static_cast<u32>(queue.size()));
    }
    if (result.index[0] != 0 || result.sizes[0] != 1)
        throw std::runtime_error("component base orbit is not singleton");
    return result;
}

u64 checkedPower(u64 base, u32 exponent, const std::string& description) {
    u64 value = 1;
    for (u32 index = 0; index < exponent; ++index) {
        if (value > std::numeric_limits<u64>::max() / base)
            throw std::runtime_error(description + " overflow");
        value *= base;
    }
    return value;
}

struct TypeGenerator {
    Permutation top_inverse;
    std::vector<Permutation> orbit_maps;
};

std::vector<TypeGenerator> inducedGenerators(
    const Layer& layer, const Action& action, const ComponentOrbits& orbits,
    const std::vector<Permutation>& transporters) {
    const u32 n = layer.component_degree;
    const u32 rank = static_cast<u32>(orbits.sizes.size());
    std::vector<TypeGenerator> result;
    for (const auto& generator : action.quotient_generators) {
        TypeGenerator induced;
        induced.top_inverse = inverse(generator.top);
        induced.orbit_maps.resize(layer.exponent, Permutation(rank));
        for (u32 source = 0; source < layer.exponent; ++source) {
            const u32 base_image = generator.components[source][0];
            const Permutation corrected = compose(
                generator.components[source], transporters[base_image]);
            if (corrected[0] != 0)
                throw std::runtime_error("corrected quotient component does not fix base");
            for (u32 orbit = 0; orbit < rank; ++orbit) {
                const u32 target = orbits.index[corrected[orbits.representatives[orbit]]];
                induced.orbit_maps[source][orbit] = target;
                if (orbits.sizes[target] != orbits.sizes[orbit])
                    throw std::runtime_error("quotient component changes orbit size");
            }
            for (u32 point = 0; point < n; ++point) {
                const u32 source_orbit = orbits.index[point];
                if (orbits.index[corrected[point]] != induced.orbit_maps[source][source_orbit])
                    throw std::runtime_error("quotient component does not permute K-orbits");
            }
            validatePermutation(induced.orbit_maps[source], rank, "induced orbit map");
        }
        result.push_back(std::move(induced));
    }
    return result;
}

u64 imageType(u64 code, const TypeGenerator& generator, u32 rank, u32 exponent) {
    std::vector<u32> digits(exponent);
    u64 remaining = code;
    for (u32 coordinate = 0; coordinate < exponent; ++coordinate) {
        digits[coordinate] = static_cast<u32>(remaining % rank);
        remaining /= rank;
    }
    u64 image = 0, place = 1;
    for (u32 output = 0; output < exponent; ++output) {
        const u32 source = generator.top_inverse[output];
        image += static_cast<u64>(generator.orbit_maps[source][digits[source]]) * place;
        place *= rank;
    }
    return image;
}

u64 typeWeight(u64 code, const ComponentOrbits& orbits, u32 exponent) {
    const u32 rank = static_cast<u32>(orbits.sizes.size());
    u64 weight = 1;
    for (u32 coordinate = 0; coordinate < exponent; ++coordinate) {
        weight *= orbits.sizes[static_cast<u32>(code % rank)];
        code /= rank;
    }
    return weight;
}

std::string jsonEscape(const std::string& text) {
    std::ostringstream out;
    for (unsigned char c : text) {
        if (c == '\\' || c == '"') out << '\\' << c;
        else if (c < 0x20) throw std::runtime_error("control character in label");
        else out << static_cast<char>(c);
    }
    return out.str();
}

void classify(const Layer& layer, const Action& action,
              const ComponentOrbits& orbits,
              const std::vector<Permutation>& transporters) {
    const auto started = std::chrono::steady_clock::now();
    const u32 rank = static_cast<u32>(orbits.sizes.size());
    const u64 degree = checkedPower(layer.component_degree, layer.exponent, "degree");
    const u64 type_degree = checkedPower(rank, layer.exponent, "type degree");
    const u64 component_stabilizer_order = layer.socle_order / layer.component_degree;
    if (component_stabilizer_order * layer.component_degree != layer.socle_order)
        throw std::runtime_error("socle order is not divisible by component degree");
    const u64 expected_h_order = checkedPower(
        component_stabilizer_order, layer.exponent, "base group stabilizer order") *
        action.quotient_order;
    if (expected_h_order != action.stabilizer_order)
        throw std::runtime_error("stabilizer order does not match product-action contract");

    const auto generators = inducedGenerators(layer, action, orbits, transporters);
    std::vector<unsigned char> seen(static_cast<std::size_t>(type_degree), 0);
    std::vector<unsigned char> regular_type(static_cast<std::size_t>(type_degree), 0);
    std::vector<u64> queue;
    std::vector<u64> orbit_representatives, orbit_masses;
    u64 orbit_count = 0, regular_orbits = 0, regular_points = 0;
    u64 mass = 0, maximum_orbit = 0;
    for (u64 seed = 0; seed < type_degree; ++seed) {
        if (seen[static_cast<std::size_t>(seed)]) continue;
        ++orbit_count;
        orbit_representatives.push_back(seed);
        queue.clear(); queue.push_back(seed); seen[static_cast<std::size_t>(seed)] = 1;
        const u64 weight = typeWeight(seed, orbits, layer.exponent);
        for (std::size_t head = 0; head < queue.size(); ++head) {
            for (const auto& generator : generators) {
                const u64 image = imageType(queue[head], generator, rank, layer.exponent);
                if (image >= type_degree) throw std::runtime_error("type image out of range");
                if (typeWeight(image, orbits, layer.exponent) != weight)
                    throw std::runtime_error("type generator changes product orbit size");
                if (seen[static_cast<std::size_t>(image)]) continue;
                seen[static_cast<std::size_t>(image)] = 1; queue.push_back(image);
            }
        }
        if (action.quotient_order % queue.size() != 0)
            throw std::runtime_error("type orbit length does not divide quotient order");
        const u64 h_orbit_size = weight * queue.size();
        orbit_masses.push_back(h_orbit_size);
        mass += h_orbit_size;
        maximum_orbit = std::max(maximum_orbit, h_orbit_size);
        if (h_orbit_size == action.stabilizer_order) {
            ++regular_orbits;
            regular_points += h_orbit_size;
            for (const u64 code : queue)
                regular_type[static_cast<std::size_t>(code)] = 1;
        }
    }
    if (mass != degree) throw std::runtime_error("H-orbit mass does not equal degree");
    std::string disposition;
    if (regular_orbits == 0) disposition = "base_gt_2";
    else if (regular_points == degree - 1) disposition = "complete";
    else if (2 * regular_points > degree) disposition = "density";
    else disposition = "graph_residual";

    bool common_neighbour_all_pairs = disposition != "graph_residual";
    u64 common_neighbour_failures = 0;
    u64 witness_type_code = std::numeric_limits<u64>::max();
    if (disposition == "graph_residual") {
        // For a target component K-orbit d, choose its stored representative
        // beta and the certified transporter t with beta^t=0.  The Boolean
        // relation relation[d][a][b] records whether some x in K-orbit a has
        // x^t in K-orbit b.  Coordinate independence then gives an exact
        // product relation between K^k-orbit types, without visiting Delta^k.
        std::vector<std::vector<std::vector<u32>>> relation(
            rank, std::vector<std::vector<u32>>(rank));
        for (u32 target = 0; target < rank; ++target) {
            const Permutation& transporter =
                transporters[orbits.representatives[target]];
            std::vector<std::vector<unsigned char>> present(
                rank, std::vector<unsigned char>(rank, 0));
            for (u32 point = 0; point < layer.component_degree; ++point) {
                const u32 source_orbit = orbits.index[point];
                const u32 image_orbit = orbits.index[transporter[point]];
                present[source_orbit][image_orbit] = 1;
            }
            for (u32 source = 0; source < rank; ++source)
                for (u32 image = 0; image < rank; ++image)
                    if (present[source][image]) relation[target][source].push_back(image);
        }
        std::vector<u64> regular_codes;
        for (u64 code = 0; code < type_degree; ++code)
            if (regular_type[static_cast<std::size_t>(code)]) regular_codes.push_back(code);
        if (regular_codes.empty()) throw std::runtime_error("graph residual has empty R");

        std::vector<u32> target_digits(layer.exponent), source_digits(layer.exponent);
        auto decode = [rank, &layer](u64 code, std::vector<u32>& digits) {
            for (u32 coordinate = 0; coordinate < layer.exponent; ++coordinate) {
                digits[coordinate] = static_cast<u32>(code % rank); code /= rank;
            }
        };
        for (std::size_t target_index = 0;
             target_index < orbit_representatives.size(); ++target_index) {
            const u64 target_code = orbit_representatives[target_index];
            decode(target_code, target_digits);
            bool common = false;
            for (const u64 source_code : regular_codes) {
                decode(source_code, source_digits);
                std::function<bool(u32, u64, u64)> search = [&](u32 coordinate,
                                                                 u64 image_code,
                                                                 u64 place) {
                    if (coordinate == layer.exponent)
                        return regular_type[static_cast<std::size_t>(image_code)] != 0;
                    const auto& images = relation[target_digits[coordinate]]
                                                 [source_digits[coordinate]];
                    for (const u32 image : images)
                        if (search(coordinate + 1, image_code + place * image,
                                   place * rank)) return true;
                    return false;
                };
                if (search(0, 0, 1)) { common = true; break; }
            }
            if (!common) {
                common_neighbour_failures += orbit_masses[target_index];
                if (witness_type_code == std::numeric_limits<u64>::max())
                    witness_type_code = target_code;
            }
        }
        common_neighbour_all_pairs = common_neighbour_failures == 0;
        disposition = common_neighbour_all_pairs ? "graph_certificate"
                                                 : "graph_counterexample";
    }
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout << "{\"schema\":\"PA8_STRUCTURED_ORBIT_SUMMARY_V1\""
              << ",\"label\":\"" << jsonEscape(action.label) << "\""
              << ",\"degree\":" << degree
              << ",\"stabilizer_order_decimal\":\"" << action.stabilizer_order << "\""
              << ",\"quotient_order\":" << action.quotient_order
              << ",\"component_rank\":" << rank
              << ",\"type_degree\":" << type_degree
              << ",\"stabilizer_orbits\":" << orbit_count
              << ",\"regular_orbits\":" << regular_orbits
              << ",\"regular_points\":" << regular_points
              << ",\"maximum_stabilizer_orbit\":" << maximum_orbit
              << ",\"disposition\":\"" << disposition << "\""
              << ",\"common_neighbour_all_pairs\":";
    if (regular_orbits == 0) std::cout << "null";
    else std::cout << (common_neighbour_all_pairs ? "true" : "false");
    std::cout << ",\"common_neighbour_failures_from_basepoint\":";
    if (regular_orbits == 0) std::cout << "null";
    else std::cout << common_neighbour_failures;
    std::cout << ",\"witness_type_code\":";
    if (witness_type_code == std::numeric_limits<u64>::max()) std::cout << "null";
    else std::cout << witness_type_code;
    std::cout
              << ",\"seconds\":" << seconds << "}\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    try {
        u64 selected_action = 0;
        if (argc == 3 && std::string(argv[1]) == "--action-index")
            selected_action = checkedInteger(
                argv[2], std::numeric_limits<u64>::max(), "selected action index");
        else if (argc != 1)
            throw std::runtime_error(
                "usage: pa_structured_orbit_engine [--action-index ONE_BASED_INDEX]");
        const Layer layer = readLayer(std::cin);
        if (selected_action > layer.actions.size())
            throw std::runtime_error("selected action index exceeds pack size");
        const ComponentOrbits orbits = componentOrbits(layer);
        const auto transporters = componentTransportersToZero(layer);
        for (std::size_t index = 0; index < layer.actions.size(); ++index) {
            if (selected_action != 0 && index + 1 != selected_action) continue;
            classify(layer, layer.actions[index], orbits, transporters);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "pa_structured_orbit_engine: " << error.what() << '\n';
        return 1;
    }
}

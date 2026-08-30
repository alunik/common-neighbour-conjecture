// Expand a compact exact product-action description into the existing
// PRIMITIVE_SAXL_V1 permutation-stream contract.  The downstream scientific
// engine is intentionally unchanged.

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <queue>
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

struct StructuredAction {
    std::string label;
    std::string stabilizer_order;
    std::string quotient_order;
    u64 smallq_check_index = 0;
    std::vector<QuotientGenerator> quotient_generators;
};

struct Layer {
    u32 component_degree = 0;
    u32 exponent = 0;
    std::string socle_order;
    std::string outer_order;
    std::vector<Permutation> socle_generators;
    std::vector<Permutation> socle_stabilizer_generators;
    std::vector<StructuredAction> actions;
};

u64 checkedInteger(const std::string& text, u64 maximum,
                   const std::string& description, bool positive = true) {
    if (text.empty()) throw std::runtime_error("missing " + description);
    u64 value = 0;
    for (const unsigned char c : text) {
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
        throw std::runtime_error(description + " has the wrong degree");
    std::vector<unsigned char> seen(degree, 0);
    for (const u32 image : permutation) {
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
        if (!(input >> image)) throw std::runtime_error("missing top generator image");
        const u64 value = checkedInteger(image, degree, "top generator image");
        permutation[point] = static_cast<u32>(value - 1);
    }
    validatePermutation(permutation, degree, "top generator");
    return permutation;
}

Layer readLayer(std::istream& input) {
    expect(input, "PRIMITIVE_SAXL_PA_STRUCTURED_V1");
    expect(input, "layer");
    Layer layer;
    std::string value;
    expect(input, "component_degree");
    if (!(input >> value)) throw std::runtime_error("missing component degree");
    layer.component_degree = static_cast<u32>(checkedInteger(
        value, std::numeric_limits<u32>::max(), "component degree"));
    expect(input, "exponent");
    if (!(input >> value)) throw std::runtime_error("missing exponent");
    layer.exponent = static_cast<u32>(checkedInteger(value, 32, "exponent"));
    expect(input, "socle_order");
    if (!(input >> layer.socle_order)) throw std::runtime_error("missing socle order");
    checkedInteger(layer.socle_order, std::numeric_limits<u64>::max(), "socle order");
    expect(input, "outer_order");
    if (!(input >> layer.outer_order)) throw std::runtime_error("missing outer order");
    checkedInteger(layer.outer_order, std::numeric_limits<u64>::max(), "outer order");

    u64 count = 0;
    expect(input, "socle_gens");
    if (!(input >> count) || count == 0 || count > 100)
        throw std::runtime_error("bad socle generator count");
    for (u64 index = 0; index < count; ++index)
        layer.socle_generators.push_back(readPackedPermutation(
            input, layer.component_degree, "component_gen"));

    expect(input, "socle_stabilizer_gens");
    if (!(input >> count) || count > 100)
        throw std::runtime_error("bad socle stabilizer generator count");
    for (u64 index = 0; index < count; ++index)
        layer.socle_stabilizer_generators.push_back(readPackedPermutation(
            input, layer.component_degree, "component_gen"));

    u64 action_count = 0;
    expect(input, "actions");
    if (!(input >> action_count) || action_count == 0 || action_count > 1000000)
        throw std::runtime_error("bad action count");
    layer.actions.reserve(static_cast<std::size_t>(action_count));
    for (u64 action_index = 0; action_index < action_count; ++action_index) {
        expect(input, "action");
        StructuredAction action;
        expect(input, "label");
        if (!(input >> action.label) || action.label.empty())
            throw std::runtime_error("missing action label");
        expect(input, "stabilizer_order");
        if (!(input >> action.stabilizer_order))
            throw std::runtime_error("missing stabilizer order");
        checkedInteger(action.stabilizer_order, std::numeric_limits<u64>::max(),
                       "stabilizer order");
        expect(input, "quotient_order");
        if (!(input >> action.quotient_order))
            throw std::runtime_error("missing quotient order");
        checkedInteger(action.quotient_order, std::numeric_limits<u64>::max(),
                       "quotient order");
        expect(input, "smallq_check_index");
        if (!(input >> action.smallq_check_index) || action.smallq_check_index == 0)
            throw std::runtime_error("bad small quotient check index");
        u64 quotient_generator_count = 0;
        expect(input, "quotient_gens");
        if (!(input >> quotient_generator_count) ||
            quotient_generator_count == 0 || quotient_generator_count > 100)
            throw std::runtime_error("bad quotient generator count");
        for (u64 generator_index = 0;
             generator_index < quotient_generator_count; ++generator_index) {
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

// Right-action composition: result[x] = (x^left)^right.
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
    parent[0] = 0;
    queue.push(0);
    while (!queue.empty()) {
        const u32 point = queue.front();
        queue.pop();
        for (u32 move = 0; move < moves.size(); ++move) {
            const u32 image = moves[move][point];
            if (parent[image] != std::numeric_limits<u32>::max()) continue;
            parent[image] = point;
            parent_move[image] = move;
            queue.push(image);
        }
    }
    if (std::find(parent.begin(), parent.end(), std::numeric_limits<u32>::max()) !=
        parent.end())
        throw std::runtime_error("component socle generators are not transitive");

    std::vector<Permutation> inverse_moves;
    inverse_moves.reserve(moves.size());
    for (const auto& move : moves) inverse_moves.push_back(inverse(move));
    std::vector<Permutation> transporters(n, Permutation(n));
    const Permutation identity = [&] {
        Permutation value(n);
        for (u32 point = 0; point < n; ++point) value[point] = point;
        return value;
    }();
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

u32 packedWidth(u32 degree) {
    u32 width = 1;
    u64 capacity = 94;
    while (capacity < degree) {
        ++width;
        capacity *= 94;
    }
    return width;
}

template <class ImageFunction>
void writePackedGenerator(u32 degree, ImageFunction imageFunction) {
    const u32 width = packedWidth(degree);
    std::cout << "packed_gen " << width << ' ';
    std::string buffer;
    buffer.reserve(1U << 20);
    for (u32 point = 0; point < degree; ++point) {
        u64 image = imageFunction(point);
        if (image >= degree) throw std::runtime_error("expanded image is out of range");
        for (u32 digit = 0; digit < width; ++digit) {
            buffer.push_back(static_cast<char>(33 + image % 94));
            image /= 94;
        }
        if (image != 0) throw std::runtime_error("expanded image exceeds packed width");
        if (buffer.size() >= (1U << 20)) {
            std::cout.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
            buffer.clear();
        }
    }
    if (!buffer.empty())
        std::cout.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    std::cout << '\n';
}

u32 productDegree(u32 n, u32 k) {
    u64 degree = 1;
    for (u32 index = 0; index < k; ++index) {
        degree *= n;
        if (degree > std::numeric_limits<u32>::max())
            throw std::runtime_error("product degree exceeds uint32");
    }
    return static_cast<u32>(degree);
}

void writeCoordinateLift(const Layer& layer, u32 coordinate,
                         const Permutation& component) {
    const u32 degree = productDegree(layer.component_degree, layer.exponent);
    u64 stride = 1;
    for (u32 index = 0; index < coordinate; ++index)
        stride *= layer.component_degree;
    writePackedGenerator(degree, [&](u32 point) -> u32 {
        const u32 digit = static_cast<u32>((point / stride) % layer.component_degree);
        const std::int64_t delta = static_cast<std::int64_t>(component[digit]) - digit;
        return static_cast<u32>(static_cast<std::int64_t>(point) +
                                delta * static_cast<std::int64_t>(stride));
    });
}

void writeCorrectedQuotientLift(
    const Layer& layer, const QuotientGenerator& generator,
    const std::vector<Permutation>& transporters) {
    const u32 n = layer.component_degree;
    const u32 k = layer.exponent;
    const u32 degree = productDegree(n, k);
    const Permutation top_inverse = inverse(generator.top);
    std::vector<u64> strides(k, 1);
    for (u32 coordinate = 1; coordinate < k; ++coordinate)
        strides[coordinate] = strides[coordinate - 1] * n;
    std::vector<const Permutation*> corrections(k);
    for (u32 output_coordinate = 0; output_coordinate < k; ++output_coordinate) {
        const u32 source_coordinate = top_inverse[output_coordinate];
        const u32 base_image = generator.components[source_coordinate][0];
        corrections[output_coordinate] = &transporters[base_image];
    }
    writePackedGenerator(degree, [&](u32 point) -> u32 {
        u64 image = 0;
        for (u32 output_coordinate = 0; output_coordinate < k; ++output_coordinate) {
            const u32 source_coordinate = top_inverse[output_coordinate];
            const u32 digit = static_cast<u32>(
                (point / strides[source_coordinate]) % n);
            const u32 lifted = generator.components[source_coordinate][digit];
            const u32 corrected = (*corrections[output_coordinate])[lifted];
            image += static_cast<u64>(corrected) * strides[output_coordinate];
        }
        return static_cast<u32>(image);
    });
}

void expand(const Layer& layer, u64 selected_action) {
    const u32 degree = productDegree(layer.component_degree, layer.exponent);
    const auto transporters = componentTransportersToZero(layer);
    std::cout << "PRIMITIVE_SAXL_V1\n";
    if (selected_action > layer.actions.size())
        throw std::runtime_error("selected action index exceeds pack size");
    for (std::size_t action_index = 0;
         action_index < layer.actions.size(); ++action_index) {
        if (selected_action != 0 && action_index + 1 != selected_action) continue;
        const auto& action = layer.actions[action_index];
        const u64 hgen_count =
            static_cast<u64>(layer.exponent) *
                layer.socle_stabilizer_generators.size() +
            action.quotient_generators.size();
        const u64 gen_count =
            static_cast<u64>(layer.exponent) * layer.socle_generators.size();
        if (hgen_count == 0 || gen_count == 0)
            throw std::runtime_error("empty expanded generator family");
        std::cout << "action\n"
                  << "label " << action.label << '\n'
                  << "degree " << degree << '\n'
                  << "stabilizer_order " << action.stabilizer_order << '\n'
                  << "classification compute\n"
                  << "regular_orbits 0\n"
                  << "regular_count 0\n"
                  << "hgens " << hgen_count << '\n';
        for (u32 coordinate = 0; coordinate < layer.exponent; ++coordinate)
            for (const auto& generator : layer.socle_stabilizer_generators)
                writeCoordinateLift(layer, coordinate, generator);
        for (const auto& generator : action.quotient_generators)
            writeCorrectedQuotientLift(layer, generator, transporters);
        std::cout << "gens " << gen_count << '\n';
        for (u32 coordinate = 0; coordinate < layer.exponent; ++coordinate)
            for (const auto& generator : layer.socle_generators)
                writeCoordinateLift(layer, coordinate, generator);
        std::cout << "end\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    try {
        u64 selected_action = 0;
        bool validate_only = false;
        if (argc == 3 && std::string(argv[1]) == "--action-index")
            selected_action = checkedInteger(
                argv[2], std::numeric_limits<u64>::max(), "selected action index");
        else if (argc == 2 && std::string(argv[1]) == "--validate-only")
            validate_only = true;
        else if (argc != 1)
            throw std::runtime_error(
                "usage: pa_structured_expander [--validate-only | --action-index ONE_BASED_INDEX]");
        Layer layer = readLayer(std::cin);
        if (!validate_only) expand(layer, selected_action);
        if (!std::cout) throw std::runtime_error("failed to write expanded stream");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "pa_structured_expander: " << error.what() << '\n';
        return 1;
    }
}

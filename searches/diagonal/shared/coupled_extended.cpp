// Exact closure of the four missing sign-coupled quotient classes:
//
//   PSL(2,p), p=11,13,17, k=3, graph S3 -> Out(T)=C2;
//   A7, k=3, graph S3 -> Out(A7)=C2.

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace psl_coupled_source {
#include "diagonal_starstar_psl2_prime.cpp"
#undef BG_DIAGONAL_COUPLED_ENGINE_HPP
#include "coupled_engine.hpp"
#undef BG_DIAGONAL_COUPLED_SIGN_RUNNER_HPP
#include "coupled_sign_runner.hpp"
}  // namespace psl_coupled_source

namespace a7_coupled_source {
#include "diagonal_starstar_a7.cpp"
#undef BG_DIAGONAL_COUPLED_ENGINE_HPP
#include "coupled_engine.hpp"
#undef BG_DIAGONAL_COUPLED_SIGN_RUNNER_HPP
#include "coupled_sign_runner.hpp"
}  // namespace a7_coupled_source

namespace {

struct ExtendedArguments {
  int selected_case = -1;
  bool census_only = false;
  bool progress = false;
  bool verbose_targets = false;
};

ExtendedArguments parse_extended_arguments(int argc, char** argv) {
  ExtendedArguments answer;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--case") {
      if (++index >= argc) {
        throw std::runtime_error("--case requires an integer 0..3");
      }
      answer.selected_case = std::stoi(argv[index]);
    } else if (argument == "--census") {
      answer.census_only = true;
    } else if (argument == "--progress") {
      answer.progress = true;
    } else if (argument == "--verbose-targets") {
      answer.verbose_targets = true;
    } else if (argument == "--help") {
      std::cout
          << "usage: coupled_extended [--case 0..3] [--census] "
             "[--progress] [--verbose-targets]\n"
          << "cases 0,1,2 are PSL(2,11/13/17); case 3 is A7\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown argument: " + argument);
    }
  }
  if (answer.selected_case < -1 || answer.selected_case > 3) {
    throw std::runtime_error("--case must lie in 0..3");
  }
  return answer;
}

struct ExtendedTotals {
  int cases = 0;
  int applicable = 0;
  int starstar_fail = 0;
  int hall_fail = 0;
};

void accumulate(
    const psl_coupled_source::CoupledCaseResult& result,
    const ExtendedArguments& arguments,
    ExtendedTotals& totals) {
  ++totals.cases;
  totals.applicable += result.applicable;
  if (!arguments.census_only && result.applicable) {
    totals.starstar_fail += !result.starstar;
    totals.hall_fail += !result.hall;
  }
}

void accumulate(
    const a7_coupled_source::CoupledCaseResult& result,
    const ExtendedArguments& arguments,
    ExtendedTotals& totals) {
  ++totals.cases;
  totals.applicable += result.applicable;
  if (!arguments.census_only && result.applicable) {
    totals.starstar_fail += !result.starstar;
    totals.hall_fail += !result.hall;
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const ExtendedArguments arguments =
        parse_extended_arguments(argc, argv);
    std::cout << "ENGINE|name=coupled_extended|schema=1"
              << "|selected_case=" << arguments.selected_case
              << "|mode="
              << (arguments.census_only ? "census" : "exact")
              << "|exactness=combined_generators_integer_H_orbits\n";
    psl_coupled_source::coupled_matching_self_test();
    a7_coupled_source::coupled_matching_self_test();
    std::cout << "MATCHING_SELF_TEST|implementations=2|status=PASS\n";

    ExtendedTotals totals;
    if (arguments.selected_case < 0 ||
        arguments.selected_case == 0) {
      accumulate(
          psl_coupled_source::run_supported_sign_graph(
              psl_coupled_source::psl2_prime(11), 3,
              "PSL211_k3_sign_graph",
              arguments.census_only, arguments.progress,
              arguments.verbose_targets),
          arguments, totals);
    }
    if (arguments.selected_case < 0 ||
        arguments.selected_case == 1) {
      accumulate(
          psl_coupled_source::run_supported_sign_graph(
              psl_coupled_source::psl2_prime(13), 3,
              "PSL213_k3_sign_graph",
              arguments.census_only, arguments.progress,
              arguments.verbose_targets),
          arguments, totals);
    }
    if (arguments.selected_case < 0 ||
        arguments.selected_case == 2) {
      accumulate(
          psl_coupled_source::run_supported_sign_graph(
              psl_coupled_source::psl2_prime(17), 3,
              "PSL217_k3_sign_graph",
              arguments.census_only, arguments.progress,
              arguments.verbose_targets),
          arguments, totals);
    }
    if (arguments.selected_case < 0 ||
        arguments.selected_case == 3) {
      accumulate(
          a7_coupled_source::run_supported_sign_graph(
              a7_coupled_source::alternating_group_dynamic(7), 3,
              "A7_k3_sign_graph",
              arguments.census_only, arguments.progress,
              arguments.verbose_targets),
          arguments, totals);
    }
    std::cout << "COUPLED_RUN_RESULT|cases=" << totals.cases
              << "|applicable=" << totals.applicable
              << "|starstar_fail=" << totals.starstar_fail
              << "|hall_fail=" << totals.hall_fail
              << "|status="
              << (totals.hall_fail == 0
                      ? (arguments.census_only
                             ? "CENSUS_PASS"
                             : "ALL_PASS")
                      : "COUNTEREXAMPLE_FOUND")
              << '\n';
    return totals.hall_fail == 0 ? 0 : 2;
  } catch (const std::exception& error) {
    std::cerr << "ERROR|" << error.what() << '\n';
    return 1;
  }
}

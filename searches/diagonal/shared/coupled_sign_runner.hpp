#ifndef BG_DIAGONAL_COUPLED_SIGN_RUNNER_HPP
#define BG_DIAGONAL_COUPLED_SIGN_RUNNER_HPP

// Included after coupled_engine.hpp inside a model namespace.

CoupledCaseResult run_supported_sign_graph(
    FiniteGroup group,
    int k,
    const std::string& label,
    bool census_only,
    bool progress,
    bool verbose_targets) {
  verify_group(group);
  if (group.outer_order != 2) {
    throw std::runtime_error(
        "sign graph requires outer automorphism of order two");
  }
  const StateSpace space = make_state_space(group.order, k);
  const ActionMaps maps = build_action_maps(group, space, k);
  if (!coupled_maps_commute(maps.symmetric, maps.outer)) {
    throw std::runtime_error(
        "outer and odd-top actions do not commute");
  }
  const std::vector<Code> combined =
      coupled_compose_maps(maps.symmetric, maps.outer);
  Digits origin_digits{};
  origin_digits.fill(group.identity);
  verify_action_map(
      combined, encode_state(space, origin_digits),
      label + "_combined_sign_outer");
  if (coupled_map_order(combined, 2) != 2) {
    throw std::runtime_error(
        "combined sign/outer map is not an involution");
  }
  if (maps.base.size() !=
      group.generators.size() +
          static_cast<std::size_t>(k - 2)) {
    throw std::runtime_error(
        "unexpected inner/alternating generator count");
  }

  std::vector<const std::vector<Code>*> actions;
  for (const auto& action : maps.base) actions.push_back(&action);
  actions.push_back(&combined);
  std::cout << "COUPLING_CHECK|label=" << label
            << "|class=graph_of_sign"
            << "|outer_order=2"
            << "|top=S" << k
            << "|quotient_order=" << factorial(k)
            << "|outer_top_commute=YES"
            << "|combined_order=2"
            << "|standalone_outer=NO"
            << "|standalone_odd_top=NO"
            << "|status=PASS\n";
  return run_coupled_case(
      group, space, k, label, factorial(k), actions,
      census_only, progress, verbose_targets);
}

#endif  // BG_DIAGONAL_COUPLED_SIGN_RUNNER_HPP

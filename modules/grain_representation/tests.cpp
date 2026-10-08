// SPDX-License-Identifier: LGPL-2.1-or-later
#include "identity.h"
#include "remapping.h"

#include <iostream>
#include <limits>
#include <stdexcept>

using namespace bravais::grains;

namespace
{
void require(const bool condition, const char *message)
{
  if (!condition) throw std::runtime_error(message);
}

InitialIdentity stripes(const unsigned int dimension, const unsigned int count)
{
  InitialIdentity input{{dimension, {12, 1, 1}}, count, {}, {}};
  for (unsigned int grain = 0; grain < 12; ++grain)
    {
      input.channels.push_back(grain % 3);
      input.cell_grains.push_back(grain + 1);
    }
  return input;
}

Fields initial_fields(const InitialIdentity &input)
{
  Fields fields(input.count, std::vector<double>(input.grid.node_count()));
  std::vector<unsigned int> neighbors(input.grid.node_count());
  for (std::size_t cell = 0; cell < input.grid.cell_count(); ++cell)
    for (const auto node : input.grid.corners(cell))
      {
        fields[input.channels[input.cell_grains[cell] - 1]][node] += 1;
        ++neighbors[node];
      }
  for (auto &field : fields)
    for (std::size_t node = 0; node < field.size(); ++node) field[node] /= neighbors[node];
  return fields;
}

void test_remapping(const unsigned int dimension)
{
  const auto input = stripes(dimension, 5);
  Identity identity(input);
  auto fields = initial_fields(input);
  // Nonzero negative tails are retained, not silently clipped to a threshold.
  fields[input.channels.front()][0] = -1e-200;
  identity.observe(fields);
  const auto labels = identity.cell_grains(fields);
  const auto remap = prepare_remapping(identity, fields);
  require(remap.changed_grains > 0, "Expected an actual channel remapping");
  for (unsigned int channel = 0; channel < fields.size(); ++channel)
    for (std::size_t node = 0; node < input.grid.node_count(); ++node)
      {
        const auto grain = identity.owners()[channel][node];
        if (grain != 0)
          require(fields[channel][node] == remap.fields[remap.channels[grain - 1]][node],
                  "Remapping must preserve each grain's values exactly");
      }
  identity.replace(remap.channels, remap.predicted_owners);
  identity.observe(remap.fields);
  require(identity.cell_grains(remap.fields) == labels, "Remapping changed physical labels");
  const auto before = identity.owners();
  auto invalid = remap.fields;
  invalid[0][0] = std::numeric_limits<double>::quiet_NaN();
  bool rejected = false;
  try { identity.observe(invalid); }
  catch (const std::runtime_error &) { rejected = true; }
  require(rejected && identity.owners() == before, "Invalid observations must not mutate identity");
  invalid = remap.fields;
  bool injected = false;
  for (std::size_t channel = 0; channel < before.size() && !injected; ++channel)
    for (std::size_t node = 0; node < before[channel].size() && !injected; ++node)
      if (before[channel][node] == 0)
        {
          invalid[channel][node] = 1e-200;
          injected = true;
        }
  rejected = false;
  try { identity.observe(invalid); }
  catch (const std::runtime_error &) { rejected = true; }
  require(injected && rejected && identity.owners() == before,
          "Even a tiny field outside owned support must be rejected");
}

void test_local_identity(const unsigned int dimension)
{
  const auto input = stripes(dimension, 5);
  Identity identity(input);
  const auto fields = initial_fields(input);
  identity.observe(fields);
  const auto before = identity.owners();
  require(identity.grain_at(0, {.5, .5, dimension == 3 ? .5 : 0}, .5) == 1 &&
          identity.grain_at(0, {3.5, .5, dimension == 3 ? .5 : 0}, .5) == 4,
          "A reused channel lost its distinct local grains");
  require(identity.grain_at(0, {1, 0, 0}, .5) == 1 &&
          identity.grain_at(0, {1, 0, 0}, 0) == 0,
          "Q1 boundary/zero-field identity is incorrect");
  require(identity.owners() == before, "Read-only local identity query mutated tracking");
  auto owners = before;
  for (const auto node : input.grid.corners(0)) owners[0][node] = 0;
  identity.replace(input.channels, owners);
  bool rejected = false;
  try { identity.grain_at(0, {.5, .5, 0}, 1e-200); }
  catch (const std::runtime_error &) { rejected = true; }
  require(rejected, "Nonzero field without local identity accepted");
  owners = before;
  owners[0][1] = 4;
  identity.replace(input.channels, owners);
  rejected = false;
  try { identity.grain_at(0, {.5, .5, 0}, .5); }
  catch (const std::runtime_error &) { rejected = true; }
  require(rejected, "Conflicting local identities accepted");
  require(identity.grain_at(0, {0, 0, 0}, .5) == 1,
          "Zero Q1 weights leaked an adjacent identity into a boundary query");
}

void test_insufficient_count()
{
  const auto input = stripes(2, 3);
  Identity identity(input);
  const auto fields = initial_fields(input);
  identity.observe(fields);
  const auto before = identity.owners();
  bool rejected = false;
  try { prepare_remapping(identity, fields); }
  catch (const std::runtime_error &) { rejected = true; }
  require(rejected, "Insufficient exact-support count was accepted");
  require(identity.channels() == input.channels && identity.owners() == before,
          "Failed allocation changed physical identity");
}

void test_tail_policy()
{
  const auto input = stripes(2, 5);
  Identity identity(input);
  auto fields = initial_fields(input);
  fields[0][0] = -1e-13;
  identity.observe(fields);
  const auto before = identity.owners();
  const auto remap = prepare_remapping(identity, fields, {1e-12, 100000});
  require(remap.removed.count == 1 && remap.removed.maximum == 1e-13 &&
          remap.removed.absolute_sum == 1e-13, "Tail removal must be explicitly measured");
  require(remap.fields[remap.channels[0]][0] == 0, "Declared tail cutoff was not applied");
  require(fields[0][0] == -1e-13 && identity.owners() == before,
          "Planning must not mutate source fields or identity");
}

void test_backtracking()
{
  // Three-colourable graph for which this ordering's greedy DSATUR fails.
  const std::vector<std::set<unsigned int>> graph = {
    {8,1,4,5}, {0,8,3}, {4,7}, {8,1,7}, {0,2,5,6},
    {0,4,7}, {8,4,7}, {2,3,5,6}, {0,1,3,6}};
  const auto channels = allocate_channels(graph, std::vector<unsigned int>(9, 0), 3, 100000);
  for (unsigned int grain = 0; grain < graph.size(); ++grain)
    for (const auto other : graph[grain])
      require(channels[grain] != channels[other], "Backtracking returned an invalid colouring");
  bool exhausted = false;
  try { allocate_channels(graph, std::vector<unsigned int>(9, 0), 3, 1); }
  catch (const std::runtime_error &error)
    { exhausted = std::string(error.what()).find("feasibility is unknown") != std::string::npos; }
  require(exhausted, "Search exhaustion must not claim impossibility");
}

void test_small_graph_feasibility()
{
  // Exhaust all 1024 simple five-vertex graphs against direct enumeration.
  // This checks both backtracking and unused-colour symmetry pruning.
  for (unsigned int mask = 0; mask < 1024; ++mask)
    {
      std::vector<std::set<unsigned int>> graph(5);
      unsigned int edge = 0;
      for (unsigned int a = 0; a < 5; ++a)
        for (unsigned int b = a + 1; b < 5; ++b, ++edge)
          if (mask & (1U << edge)) { graph[a].insert(b); graph[b].insert(a); }
      for (const unsigned int count : {2U, 3U})
        {
          bool feasible = false;
          unsigned int combinations = 1;
          for (unsigned int i = 0; i < 5; ++i) combinations *= count;
          for (unsigned int encoded = 0; encoded < combinations && !feasible; ++encoded)
            {
              auto remaining = encoded;
              std::vector<unsigned int> channels(5);
              for (auto &channel : channels) { channel = remaining % count; remaining /= count; }
              bool valid = true;
              for (unsigned int a = 0; a < 5; ++a)
                for (const auto b : graph[a])
                  if (channels[a] == channels[b]) valid = false;
              feasible = valid;
            }
          bool allocated = false;
          try
            {
              const auto channels = allocate_channels(graph, {0, count - 1, 0, 1, 0}, count, 100000);
              for (unsigned int a = 0; a < 5; ++a)
                for (const auto b : graph[a])
                  require(channels[a] != channels[b], "Invalid exhaustive-test colouring");
              allocated = true;
            }
          catch (const std::runtime_error &error)
            {
              require(std::string(error.what()).find("no coloring") != std::string::npos,
                      "Small graph unexpectedly exhausted search budget");
            }
          require(allocated == feasible, "Allocator disagrees with exhaustive feasibility");
        }
    }
}
}

int main()
{
  InitialIdentity independent{{2, {2, 1, 1}}, 2, {0, 1}, {1, 2}};
  Identity independent_identity(independent);
  Fields smooth{{0.9, 0.5, 0.1, 0.9, 0.5, 0.1}, {0.1, 0.5, 0.9, 0.1, 0.5, 0.9}};
  independent_identity.observe(smooth);
  require(independent_identity.cell_grains(smooth) == std::vector<unsigned int>({1, 2}),
          "Independent smooth initial tails lost grain identity");
  require(snap_grid_coordinate(3.0 - 1e-14, 12) == 3.0, "Grid roundoff was not snapped");
  require(snap_grid_coordinate(3.25, 12) == 3.25, "Real interpolation location was changed");
  test_remapping(2);
  test_remapping(3);
  test_local_identity(2);
  test_local_identity(3);
  test_insufficient_count();
  test_tail_policy();
  test_backtracking();
  test_small_graph_feasibility();
  std::cout << "grain representation: 2D/3D remap, exact values, identity and rejection passed\n";
}

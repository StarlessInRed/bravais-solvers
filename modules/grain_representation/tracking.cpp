// SPDX-License-Identifier: LGPL-2.1-or-later
#include "tracking.h"

namespace bravais::grains
{
SupportPrediction predict_support(const Identity &identity)
{
  const auto &grid = identity.grid();
  SupportPrediction prediction;
  prediction.nodes_by_grain.resize(identity.channels().size());
  prediction.conflicts.resize(identity.channels().size());
  std::vector<std::set<unsigned int>> current_at_node(grid.node_count());
  for (const auto &owners : identity.owners())
    for (std::size_t node = 0; node < grid.node_count(); ++node)
      if (owners[node] != 0) current_at_node[node].insert(owners[node] - 1);
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
    {
      const auto corners = grid.corners(cell);
      std::set<unsigned int> grains;
      for (const auto node : corners)
        grains.insert(current_at_node[node].begin(), current_at_node[node].end());
      for (const auto grain : grains)
        prediction.nodes_by_grain[grain].insert(corners.begin(), corners.end());
    }
  // Keep predicted supports element-disjoint, not only node-disjoint. This
  // protects both nonlinear quadrature and cell-centre identity after the step.
  std::vector<std::set<unsigned int>> at_node(grid.node_count());
  for (unsigned int grain = 0; grain < prediction.nodes_by_grain.size(); ++grain)
    for (const auto node : prediction.nodes_by_grain[grain])
      at_node[node].insert(grain);
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
    {
      std::set<unsigned int> grains;
      for (const auto node : grid.corners(cell))
        grains.insert(at_node[node].begin(), at_node[node].end());
      for (const auto grain : grains)
        for (const auto other : grains)
          if (grain != other) prediction.conflicts[grain].insert(other);
    }
  return prediction;
}
}

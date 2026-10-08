// SPDX-License-Identifier: LGPL-2.1-or-later
#include "identity.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>

namespace bravais::grains
{
double snap_grid_coordinate(const double coordinate, const unsigned int cells)
{
  const double rounded = std::round(coordinate);
  const double tolerance = 64 * std::numeric_limits<double>::epsilon() * std::max(1.0, double(cells));
  return std::abs(coordinate - rounded) <= tolerance ? rounded : coordinate;
}

void Grid::validate() const
{
  if ((dimension != 2 && dimension != 3) || cells[0] == 0 || cells[1] == 0 ||
      cells[2] == 0 || (dimension == 2 && cells[2] != 1))
    throw std::invalid_argument("Grain identity requires a 2D/3D nonempty grid");
  std::size_t nodes = 1;
  for (unsigned int d = 0; d < dimension; ++d)
    {
      const auto extent = std::size_t(cells[d]) + 1;
      if (extent > std::numeric_limits<std::size_t>::max() / nodes)
        throw std::invalid_argument("Grain identity grid overflows node indexing");
      nodes *= extent;
    }
}

std::size_t Grid::node_count() const
{
  return (std::size_t(cells[0]) + 1) * (std::size_t(cells[1]) + 1) *
         (dimension == 3 ? std::size_t(cells[2]) + 1 : 1);
}

std::size_t Grid::cell_count() const
{
  return std::size_t(cells[0]) * cells[1] * cells[2];
}

std::vector<std::size_t> Grid::corners(const std::size_t cell) const
{
  const std::array<std::size_t, 3> origin = {
    cell % cells[0], (cell / cells[0]) % cells[1], cell / cells[0] / cells[1]};
  std::vector<std::size_t> result;
  for (unsigned int corner = 0; corner < (1u << dimension); ++corner)
    {
      const auto x = origin[0] + (corner & 1);
      const auto y = origin[1] + ((corner >> 1) & 1);
      const auto z = origin[2] + (dimension == 3 ? (corner >> 2) & 1 : 0);
      result.push_back(x + (std::size_t(cells[0]) + 1) *
                            (y + (std::size_t(cells[1]) + 1) * z));
    }
  return result;
}

void InitialIdentity::validate() const
{
  grid.validate();
  if (count == 0 || channels.empty() || cell_grains.size() != grid.cell_count())
    throw std::invalid_argument("Invalid grain identity dimensions");
  for (const auto channel : channels)
    if (channel >= count)
      throw std::invalid_argument("Grain channel is outside fixed count");
  std::vector<bool> present(channels.size(), false);
  for (const auto grain : cell_grains)
    {
      if (grain == 0 || grain > channels.size())
        throw std::invalid_argument("Initial physical grain ID is out of range");
      present[grain - 1] = true;
    }
  for (const bool found : present)
    if (!found)
      throw std::invalid_argument("Initial grain identity has an unused grain");
}

Identity::Identity(const InitialIdentity &input)
  : mesh(input.grid), channel_by_grain(input.channels)
{
  input.validate();
  owner_by_channel_node.assign(input.count, std::vector<unsigned int>(mesh.node_count()));
  if (std::set<unsigned int>(input.channels.begin(), input.channels.end()).size() == input.channels.size())
    {
      // An independent channel has one unambiguous grain even for a smooth
      // analytical initial interface with nonzero tails outside labelled cells.
      for (unsigned int grain = 0; grain < input.channels.size(); ++grain)
        std::fill(owner_by_channel_node[input.channels[grain]].begin(),
                  owner_by_channel_node[input.channels[grain]].end(), grain + 1);
      return;
    }
  for (std::size_t cell = 0; cell < mesh.cell_count(); ++cell)
    {
      const auto grain = input.cell_grains[cell];
      auto &owners = owner_by_channel_node[channel_by_grain[grain - 1]];
      for (const auto node : mesh.corners(cell))
        {
          if (owners[node] != 0 && owners[node] != grain)
            throw std::invalid_argument("Initial grains overlap in one channel");
          owners[node] = grain;
        }
    }
}

const Grid &Identity::grid() const { return mesh; }
const std::vector<unsigned int> &Identity::channels() const { return channel_by_grain; }
const Owners &Identity::owners() const { return owner_by_channel_node; }

unsigned int Identity::grain_at(const unsigned int channel,
                                const std::array<double, 3> &point,
                                const double value) const
{
  if (channel >= owner_by_channel_node.size() || !std::isfinite(value))
    throw std::invalid_argument("Invalid local grain channel or field value");
  std::array<std::size_t, 3> lower{};
  std::array<double, 3> fraction{};
  for (unsigned int d = 0; d < mesh.dimension; ++d)
    {
      const auto coordinate = snap_grid_coordinate(point[d], mesh.cells[d]);
      if (!std::isfinite(coordinate) || coordinate < 0 || coordinate > mesh.cells[d])
        throw std::invalid_argument("Local grain query outside Q1 grid");
      lower[d] = std::min(std::size_t(std::floor(coordinate)), std::size_t(mesh.cells[d] - 1));
      fraction[d] = coordinate - lower[d];
    }
  if (value == 0) return 0;
  const auto location = [&] {
    return " channel=" + std::to_string(channel) + " grid_point=(" +
      std::to_string(point[0]) + "," + std::to_string(point[1]) + "," + std::to_string(point[2]) + ")";
  };
  unsigned int grain = 0;
  for (unsigned int corner = 0; corner < (1u << mesh.dimension); ++corner)
    {
      auto node = lower;
      double weight = 1;
      for (unsigned int d = 0; d < mesh.dimension; ++d)
        {
          const auto upper = (corner >> d) & 1;
          node[d] += upper;
          weight *= upper ? fraction[d] : 1 - fraction[d];
        }
      // On a Q1 face/node only contributing nodes establish identity. Adjacent
      // cells cannot change a boundary value's identity through zero weights.
      if (weight == 0) continue;
      const auto owner = owner_by_channel_node[channel].at(
        node[0] + (std::size_t(mesh.cells[0]) + 1) *
                    (node[1] + (std::size_t(mesh.cells[1]) + 1) * node[2]));
      if (owner == 0) continue;
      if (owner > channel_by_grain.size() || channel_by_grain[owner - 1] != channel)
        throw std::runtime_error("Local grain identity disagrees with channel mapping" + location());
      if (grain != 0 && grain != owner)
        throw std::runtime_error("Ambiguous local grain identity in Q1 support" + location());
      grain = owner;
    }
  if (grain == 0)
    throw std::runtime_error("Nonzero field has no local grain identity" + location());
  return grain;
}

void Identity::observe(const Fields &fields)
{
  if (fields.size() != owner_by_channel_node.size())
    throw std::invalid_argument("Observed field count differs from fixed count");
  for (std::size_t channel = 0; channel < fields.size(); ++channel)
    {
      if (fields[channel].size() != mesh.node_count())
        throw std::invalid_argument("Observed field grid differs from identity grid");
      for (std::size_t node = 0; node < mesh.node_count(); ++node)
        {
          const auto value = fields[channel][node];
          if (!std::isfinite(value))
            throw std::runtime_error("Nonfinite field in grain tracking");
          if (value != 0 && owner_by_channel_node[channel][node] == 0)
            throw std::runtime_error("Nonzero field outside predicted Q1 grain support");
        }
    }
  for (std::size_t channel = 0; channel < fields.size(); ++channel)
    for (std::size_t node = 0; node < mesh.node_count(); ++node)
      if (fields[channel][node] == 0)
        owner_by_channel_node[channel][node] = 0;
}

void Identity::replace(std::vector<unsigned int> channels, Owners owners)
{
  channel_by_grain = std::move(channels);
  owner_by_channel_node = std::move(owners);
}

std::vector<unsigned int> Identity::cell_grains(const Fields &fields) const
{
  std::vector<unsigned int> result;
  for (std::size_t cell = 0; cell < mesh.cell_count(); ++cell)
    {
      double best = -std::numeric_limits<double>::infinity();
      unsigned int winner = 0;
      for (std::size_t channel = 0; channel < fields.size(); ++channel)
        {
          unsigned int grain = 0;
          double sum = 0;
          for (const auto node : mesh.corners(cell))
            {
              const auto owner = owner_by_channel_node[channel][node];
              if (owner != 0 && grain != 0 && owner != grain)
                throw std::runtime_error("Two physical grains share one Q1 element/channel");
              if (owner != 0) grain = owner;
              sum += fields[channel][node];
            }
          // Q1 cell-centre weights are equal. Tie-break by grain, never channel.
          if (grain != 0 && (sum > best || (sum == best && grain < winner)))
            {
              best = sum;
              winner = grain;
            }
        }
      if (winner == 0)
        throw std::runtime_error("No physical grain at a Q1 cell centre");
      result.push_back(winner);
    }
  return result;
}
}

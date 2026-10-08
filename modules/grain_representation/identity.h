// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace bravais::grains
{
using Fields = std::vector<std::vector<double>>;
using Owners = std::vector<std::vector<unsigned int>>;

// Remove coordinate-arithmetic roundoff at integer grid nodes, not field tails.
double snap_grid_coordinate(double coordinate, unsigned int cells);

// Fixed, nonperiodic Q1 grid. Node/cell indices are x-fastest; grain 0 is absent.
struct Grid
{
  unsigned int dimension;
  std::array<unsigned int, 3> cells;
  std::size_t node_count() const;
  std::size_t cell_count() const;
  std::vector<std::size_t> corners(std::size_t cell) const;
  void validate() const;
};

struct InitialIdentity
{
  Grid grid;
  unsigned int count;
  std::vector<unsigned int> channels;
  std::vector<unsigned int> cell_grains;
  void validate() const;
};

class Identity
{
public:
  explicit Identity(const InitialIdentity &input);
  const Grid &grid() const;
  const std::vector<unsigned int> &channels() const;
  const Owners &owners() const;
  // Grid coordinates, not a nearest-grain lookup. Zero field has no identity.
  unsigned int grain_at(unsigned int channel, const std::array<double, 3> &point,
                        double value) const;
  void observe(const Fields &fields);
  void replace(std::vector<unsigned int> channels, Owners owners);
  std::vector<unsigned int> cell_grains(const Fields &fields) const;

private:
  Grid mesh;
  std::vector<unsigned int> channel_by_grain;
  Owners owner_by_channel_node;
};
}

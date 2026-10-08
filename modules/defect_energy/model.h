// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include "recovery.h"

#include <array>
#include <filesystem>
#include <span>
#include <vector>

namespace bravais::defects
{
struct Input
{
  unsigned int dimension;
  std::array<unsigned int, 3> cells;
  std::array<double, 3> extent;
  double energy_scale; // J/m^3 per dimensionless energy density
  double energy_per_density; // alpha_rho G b^2, J/m
  std::vector<unsigned int> parent_by_grain; // grain_id 1..G, not channel 0..Q-1
  std::vector<double> density; // immutable initial parent branch, m^-2, x fastest
  double time_scale_s;
  std::shared_ptr<const Recovery> recovery;
};

struct LocalState
{
  double h_parent;
  double h_rx;
  double rho_effective;
  double energy; // dimensionless
  std::vector<double> derivative;
};

Input read_input(const std::filesystem::path &path, unsigned int dimension,
                 unsigned int grain_count, const std::array<double, 3> &extent);
double density_at(const Input &input, const std::array<double, 3> &point, double time_s);
LocalState evaluate(const Input &input, std::span<const double> eta,
                    std::span<const unsigned int> local_grains, double density);
}

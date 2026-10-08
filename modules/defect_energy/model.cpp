// SPDX-License-Identifier: LGPL-2.1-or-later
#include "model.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

namespace bravais::defects
{
Input read_input(const std::filesystem::path &path, const unsigned int dimension,
                 const unsigned int grain_count, const std::array<double, 3> &extent)
{
  std::ifstream stream(path);
  Input input{};
  std::string format;
  unsigned int grains = 0;
  stream >> format >> input.dimension >> input.cells[0] >> input.cells[1] >> input.cells[2]
         >> grains >> input.energy_scale >> input.energy_per_density >> input.time_scale_s;
  if (!stream || format != "bravais_parent_rx_local_v1" || input.dimension != dimension ||
      grains != grain_count || grains == 0 || (dimension != 2 && dimension != 3) ||
      !std::isfinite(input.energy_scale) || input.energy_scale <= 0 ||
      !std::isfinite(input.energy_per_density) || input.energy_per_density <= 0 ||
      !std::isfinite(input.time_scale_s) || input.time_scale_s <= 0 ||
      (dimension == 2 && input.cells[2] != 1))
    throw std::invalid_argument("Invalid parent/RX stored-energy header");
  input.recovery = read_recovery(stream);
  input.extent = extent;
  std::size_t size = 1;
  for (unsigned int d = 0; d < 3; ++d)
    {
      if (input.cells[d] == 0 || !std::isfinite(extent[d]) || extent[d] <= 0)
        throw std::invalid_argument("Invalid stored-energy grid");
      size *= input.cells[d];
    }
  input.parent_by_grain.resize(grains);
  for (auto &parent : input.parent_by_grain)
    if (!(stream >> parent) || parent > 1)
      throw std::invalid_argument("Stored-energy category must be parent=1 or rx=0");
  input.density.resize(size);
  for (auto &rho : input.density)
    if (!(stream >> rho) || !std::isfinite(rho) || rho < 0 ||
        !std::isfinite(rho * input.energy_per_density / input.energy_scale))
      throw std::invalid_argument("Initial parent density must be finite, nonnegative and representable");
  std::string extra;
  if (stream >> extra) throw std::invalid_argument("Unexpected stored-energy input data");
  return input;
}

double density_at(const Input &input, const std::array<double, 3> &point, const double time_s)
{
  if (!std::isfinite(time_s) || time_s < 0 || !input.recovery)
    throw std::invalid_argument("Stored-energy query requires nonnegative physical time and a recovery model");
  std::array<std::size_t, 3> low{}, high{};
  std::array<double, 3> fraction{};
  for (unsigned int d = 0; d < input.dimension; ++d)
    {
      if (!std::isfinite(point[d]) || point[d] < 0 || point[d] > input.extent[d])
        throw std::invalid_argument("Stored-energy query outside supplied domain");
      // Cell-centred values linearly reconstructed between centres; the outer
      // half-cells extend their nearest centre. No values outside the domain.
      const auto coordinate = std::clamp(point[d] / input.extent[d] * input.cells[d] - 0.5,
                                         0.0, double(input.cells[d] - 1));
      low[d] = static_cast<std::size_t>(std::floor(coordinate));
      high[d] = std::min(low[d] + 1, std::size_t(input.cells[d] - 1));
      fraction[d] = coordinate - low[d];
    }
  double value = 0;
  for (unsigned int corner = 0; corner < (1u << input.dimension); ++corner)
    {
      std::array<std::size_t, 3> cell{};
      double weight = 1;
      for (unsigned int d = 0; d < input.dimension; ++d)
        {
          const bool upper = (corner >> d) & 1;
          cell[d] = upper ? high[d] : low[d];
          weight *= upper ? fraction[d] : 1 - fraction[d];
        }
      const auto initial = input.density.at(cell[0] + input.cells[0] * (cell[1] + input.cells[1] * cell[2]));
      // Evolve each supplied cell first; interpolation does not commute with a
      // nonzero recovery target when initial densities straddle that target.
      value += weight * input.recovery->density(initial, time_s);
    }
  return value;
}

LocalState evaluate(const Input &input, const std::span<const double> eta,
                    const std::span<const unsigned int> local_grains, const double density)
{
  if (eta.size() != local_grains.size() || !std::isfinite(density) || density < 0)
    throw std::invalid_argument("Stored-energy evaluation shape or density is invalid");
  double sum = 0, parent = 0;
  for (std::size_t i = 0; i < eta.size(); ++i)
    {
      if (!std::isfinite(eta[i])) throw std::runtime_error("Nonfinite order parameter");
      if ((eta[i] != 0 && local_grains[i] == 0) || local_grains[i] > input.parent_by_grain.size())
        throw std::runtime_error("Nonzero stored-energy field requires a valid local grain identity");
      sum += eta[i] * eta[i];
      if (eta[i] != 0) parent += eta[i] * eta[i] * input.parent_by_grain[local_grains[i] - 1];
    }
  if (!std::isfinite(sum) || sum <= 0)
    throw std::runtime_error("Stored-energy fraction undefined: sum(eta^2) must be positive");
  const double h_parent = parent / sum;
  const double e_parent = density * input.energy_per_density / input.energy_scale;
  LocalState state{h_parent, 1 - h_parent, h_parent * density, h_parent * e_parent, {}};
  state.derivative.resize(eta.size());
  for (std::size_t i = 0; i < eta.size(); ++i)
    {
      const auto branch = eta[i] == 0 ? 0 : input.parent_by_grain[local_grains[i] - 1];
      state.derivative[i] = 2 * eta[i] / sum * (branch * e_parent - state.energy);
    }
  return state;
}
}

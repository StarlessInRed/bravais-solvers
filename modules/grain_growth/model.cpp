// SPDX-License-Identifier: LGPL-2.1-or-later
#include "model.h"

namespace phasefield::growth
{
double potential(const std::span<const double> eta, const double interaction)
{
  double energy = 0.25;
  for (std::size_t i = 0; i < eta.size(); ++i)
    {
      const auto square = eta[i] * eta[i];
      energy += square * square / 4 - square / 2;
      for (std::size_t j = 0; j < i; ++j)
        energy += interaction * square * eta[j] * eta[j];
    }
  return energy;
}

double potential_derivative(const double eta, const double sum_squares, const double interaction)
{
  auto derivative = eta * eta * eta - eta;
  derivative += 2.0 * interaction * eta * (sum_squares - eta * eta);
  return derivative;
}
}

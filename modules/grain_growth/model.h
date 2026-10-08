// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <span>

namespace phasefield::growth
{
// Dimensionless multiphase potential; the constant offset does not affect its derivative.
double potential(std::span<const double> eta, double interaction);
double potential_derivative(double eta, double sum_squares, double interaction);
}

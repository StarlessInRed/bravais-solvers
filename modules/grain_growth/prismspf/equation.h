// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <prismspf/core/pde_operator_base.h>
#include <array>
#include <span>
#include <vector>

namespace phasefield::growth
{
// Ordinary Allen-Cahn update and Q1 nodal input. Callers own tracking and any extra driving force.
template <unsigned int dim>
class Equation
{
public:
  using Value = dealii::VectorizedArray<double>;
  Equation(const prismspf::UserInputParameters<dim> &input, unsigned int count);
  const std::array<unsigned int, 3> &grid_cells() const;
  double initial_value(unsigned int field, const dealii::Point<dim> &point) const;
  void compute_rhs(prismspf::FieldContainer<dim, 1, double> &fields,
                   const prismspf::SimulationTimer &timer,
                   std::span<const Value> extra_derivatives) const;

  const unsigned int count;
  const double mobility, gradient, interaction;

private:
  std::array<unsigned int, 3> cells;
  std::array<double, dim> extent;
  std::vector<double> initial;
};

extern template class Equation<2>;
extern template class Equation<3>;
}

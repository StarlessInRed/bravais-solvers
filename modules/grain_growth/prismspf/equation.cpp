// SPDX-License-Identifier: LGPL-2.1-or-later
#include "equation.h"
#include "../model.h"
#include <grain_representation/identity.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>

namespace phasefield::growth
{
template <unsigned int dim>
Equation<dim>::Equation(const prismspf::UserInputParameters<dim> &input, const unsigned int count)
  : count(count), mobility(input.user_constants.get_double("mobility")),
    gradient(input.user_constants.get_double("gradient_energy")),
    interaction(input.user_constants.get_double("interaction"))
{
  AssertThrow(count > 0 && std::isfinite(mobility) && mobility > 0 &&
                std::isfinite(gradient) && gradient > 0 &&
                std::isfinite(interaction) && interaction > 0,
              dealii::ExcMessage("Field count and model coefficients must be positive and finite"));
  for (unsigned int d = 0; d < dim; ++d)
    extent[d] = input.spatial_discretization.rectangular_mesh.size[d];
  std::ifstream file(input.user_constants.get_string("initial_fields"));
  std::string format;
  unsigned int fields = 0;
  file >> format >> cells[0] >> cells[1] >> cells[2] >> fields;
  AssertThrow(file && format == "bravais_nodal_fields_v1" && fields == count,
              dealii::ExcMessage("Invalid initial fields header"));
  AssertThrow(dim == 3 || cells[2] == 1, dealii::ExcMessage("2D input requires one storage layer"));
  std::size_t nodes = 1;
  for (unsigned int d = 0; d < dim; ++d)
    {
      AssertThrow(cells[d] > 0, dealii::ExcMessage("Initial grid axis must be positive"));
      nodes *= cells[d] + 1;
    }
  initial.resize(nodes * count);
  for (auto &value : initial)
    {
      file >> value;
      AssertThrow(file && std::isfinite(value) && value >= 0 && value <= 1,
                  dealii::ExcMessage("Invalid initial field value"));
    }
  std::string extra;
  AssertThrow(!(file >> extra), dealii::ExcMessage("Unexpected initial field data"));
}

template <unsigned int dim>
const std::array<unsigned int, 3> &Equation<dim>::grid_cells() const { return cells; }

template <unsigned int dim>
double Equation<dim>::initial_value(const unsigned int field, const dealii::Point<dim> &point) const
{
  std::array<unsigned int, 3> lower = {};
  std::array<double, 3> fraction = {};
  for (unsigned int d = 0; d < dim; ++d)
    {
      const double coordinate = bravais::grains::snap_grid_coordinate(
        std::clamp(point[d] / extent[d] * cells[d], 0.0, double(cells[d])), cells[d]);
      lower[d] = std::min(unsigned(std::floor(coordinate)), cells[d] - 1);
      fraction[d] = coordinate - lower[d];
    }
  double value = 0;
  for (unsigned int corner = 0; corner < (1u << dim); ++corner)
    {
      std::array<unsigned int, 3> node = {};
      double weight = 1;
      for (unsigned int d = 0; d < dim; ++d)
        {
          const auto offset = (corner >> d) & 1;
          node[d] = lower[d] + offset;
          weight *= offset ? fraction[d] : 1 - fraction[d];
        }
      const auto index = node[0] + (cells[0] + 1) * (node[1] + (cells[1] + 1) * node[2]);
      value += weight * initial[index * count + field];
    }
  return value;
}

template <unsigned int dim>
void Equation<dim>::compute_rhs(prismspf::FieldContainer<dim, 1, double> &fields,
                                const prismspf::SimulationTimer &timer,
                                const std::span<const Value> extra_derivatives) const
{
  AssertThrow(extra_derivatives.size() == count, dealii::ExcMessage("Driving force count differs from Q"));
  std::vector<Value> terms(count);
  std::vector<dealii::Tensor<1, dim, Value>> gradients(count);
  Value sum_squares = {};
  for (unsigned int i = 0; i < count; ++i)
    {
      const auto ni = fields.template get_value<prismspf::Scalar, prismspf::OldOne>(i);
      sum_squares += ni * ni;
    }
  for (unsigned int i = 0; i < count; ++i)
    {
      const auto ni = fields.template get_value<prismspf::Scalar, prismspf::OldOne>(i);
      Value derivative;
      for (unsigned int lane = 0; lane < Value::size(); ++lane)
        derivative[lane] = potential_derivative(ni[lane], sum_squares[lane], interaction);
      derivative += extra_derivatives[i];
      terms[i] = ni - timer.get_timestep() * mobility * derivative;
      gradients[i] = -timer.get_timestep() * mobility * gradient *
                     fields.template get_gradient<prismspf::Scalar, prismspf::OldOne>(i);
    }
  // Read all OldOne fields before publishing any residual term.
  for (unsigned int i = 0; i < count; ++i)
    {
      fields.set_value_term(i, terms[i]);
      fields.set_gradient_term(i, gradients[i]);
    }
}

template class Equation<2>;
template class Equation<3>;
}

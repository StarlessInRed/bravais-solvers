// SPDX-License-Identifier: LGPL-2.1-or-later
// Shared stored-energy Allen-Cahn operator. Apps choose the grain evolution policy.
#pragma once
#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/pde_operator_base.h>
#include <prismspf/core/problem.h>
#include <prismspf/core/timer.h>
#include <grain_growth/prismspf/equation.h>
#include <grain_representation/io.h>
#include <grain_representation/prismspf/field_access.h>
#include <defect_energy/model.h>
#include <defect_energy/prismspf/diagnostics.h>
#include <array>
#include <algorithm>
#include <iostream>
#include <optional>
#include <set>
#include <string>

using namespace prismspf;

template <unsigned int dim>
class StoredEnergyGrowth : public PDEOperatorBase<dim, 1, double>
{
public:
  StoredEnergyGrowth(const UserInputParameters<dim> &input, PhaseFieldTools<dim> &tools,
              const unsigned int count)
    : PDEOperatorBase<dim, 1, double>(input, tools), count(count), equation(input, count),
      tracking(bravais::grains::read_identity(input.user_constants.get_string("initial_identity")))
  {
    const auto names = UserConstants<dim>::get_names(input.user_constants.file_name);
    if (names.contains("stored_energy"))
      {
        const auto path = input.user_constants.get_string("stored_energy");
        const auto identity = bravais::grains::read_identity(input.user_constants.get_string("initial_identity"));
        std::array<double, 3> extent{1, 1, 1};
        for (unsigned int d = 0; d < dim; ++d)
          extent[d] = input.spatial_discretization.rectangular_mesh.size[d];
        stored_energy = bravais::defects::read_input(path, dim, identity.channels.size(), extent);
        AssertThrow(stored_energy->cells == equation.grid_cells(), dealii::ExcMessage("Stored energy grid differs from initial fields"));
      }

  }

protected:
  virtual void update_grains(SolveContext<dim, 1, double> &context) = 0;

private:
  void pre_solve_block(SolveContext<dim, 1, double> &context, unsigned int) override
  {
    tracking.before_step(context);
  }

  void post_solve_block(SolveContext<dim, 1, double> &context, unsigned int) override
  {
    tracking.observe_step(context);
    update_grains(context);
    tracking.write_step(context);
    if (stored_energy)
      bravais::defects::write_diagnostics(context, *stored_energy, tracking, equation.interaction, equation.gradient);
  }

  void compute_rhs(FieldContainer<dim, 1, double> &fields,
                   const SimulationTimer &timer, unsigned int) const override
  {
    using Value = dealii::VectorizedArray<double>;
    Value sum_squares = {};
    std::vector<Value> defect_derivatives(count);
    for (auto &derivative : defect_derivatives) derivative = 0.0;
    for (unsigned int i = 0; i < count; ++i)
      {
        const auto ni = fields.template get_value<prismspf::Scalar, OldOne>(i);
        sum_squares += ni * ni;
      }
    if (stored_energy)
      {
        // v4.0.1 advances the timer before this explicit update. OldOne is at
        // t_n, while timer is t_(n+1); cycle zero only sets initial conditions.
        const auto time_s = std::max(0.0, timer.get_time() - timer.get_timestep()) * stored_energy->time_scale_s;
        const auto location = fields.get_q_point_location();
        std::vector<double> eta(count);
        for (unsigned int lane = 0; lane < Value::size(); ++lane)
          {
            // Padding lanes may be zero or duplicate a real lane. Zero lanes
            // have no identity. Every active Q1 Lobatto state is checked by
            // diagnostics; duplicate padding only performs the same const query.
            if (sum_squares[lane] == 0) continue;
            std::array<double, 3> point{};
            for (unsigned int d = 0; d < dim; ++d) point[d] = location[d][lane];
            for (unsigned int i = 0; i < count; ++i)
              eta[i] = fields.template get_value<prismspf::Scalar, OldOne>(i)[lane];
            const auto state = bravais::defects::evaluate(*stored_energy, eta, tracking.local_grains(point, eta),
              bravais::defects::density_at(*stored_energy, point, time_s));
            for (unsigned int i = 0; i < count; ++i) defect_derivatives[i][lane] = state.derivative[i];
          }
      }
    equation.compute_rhs(fields, timer, defect_derivatives);
  }

  void set_initial_condition(const unsigned int &field, const unsigned int &,
                             const dealii::Point<dim> &point, double &value,
                             double &) const override
  {
    value = equation.initial_value(field, point);
  }

protected:
  const unsigned int count;
  phasefield::growth::Equation<dim> equation;
  bravais::grains::NativeTracking<dim> tracking;
  std::optional<bravais::defects::Input> stored_energy;
};


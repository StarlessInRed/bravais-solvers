// SPDX-License-Identifier: LGPL-2.1-or-later
#include "diagnostics.h"
#include <grain_growth/model.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <prismspf/core/system_wide.h>
#include <fstream>
#include <iomanip>

namespace bravais::defects
{
template <unsigned int dim>
void write_diagnostics(prismspf::SolveContext<dim, 1, double> &context,
                       const Input &input, const grains::NativeTracking<dim> &tracking,
                       const double interaction, const double gradient)
{
  const auto &timer = context.get_simulation_timer();
  const auto &output = context.get_user_inputs().output_parameters;
  const bool save = output.should_output(timer.get_increment());
  const auto count = context.get_field_attributes().size();
  const auto &first = context.get_dof_manager().get_field_dof_handler(0);
  const auto &quadrature = prismspf::SystemWide<dim, 1>::quadrature;
  dealii::FEValues<dim> fe(first.get_fe(), quadrature,
    dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points | dealii::update_JxW_values);
  std::vector<typename dealii::DoFHandler<dim>::active_cell_iterator> cells;
  for (unsigned int i = 0; i < count; ++i)
    cells.push_back(context.get_dof_manager().get_field_dof_handler(i).begin_active());
  std::vector<std::vector<double>> values(count, std::vector<double>(quadrature.size()));
  std::vector<std::vector<dealii::Tensor<1, dim>>> gradients(count, std::vector<dealii::Tensor<1, dim>>(quadrature.size()));
  const auto path = std::filesystem::path(output.folder) /
    (output.file_name + "-energy-" + std::to_string(timer.get_increment()) + ".csv");
  std::ofstream stream;
  if (save)
    {
      stream.open(path);
      if (!stream) throw std::runtime_error("Cannot write stored-energy diagnostics");
      stream << std::setprecision(17);
      stream << "cycle,time_star,time_s,x_star,y_star,z_star,weight_star,rho_p_m2,rho_eff_m2,h_p,h_rx,e_loc_j_m3,f_well_j_m3,f_grad_j_m3\n";
    }
  for (; cells[0] != first.end(); )
    {
      for (unsigned int i = 0; i < count; ++i)
        {
          fe.reinit(cells[i]);
          const auto &solution = context.get_solution_indexer().get_solution_vector(i);
          fe.get_function_values(solution, values[i]);
          fe.get_function_gradients(solution, gradients[i]);
        }
      for (unsigned int q = 0; q < quadrature.size(); ++q)
        {
          std::array<double, 3> point{};
          for (unsigned int d = 0; d < dim; ++d) point[d] = fe.quadrature_point(q)[d];
          const auto time_s = timer.get_time() * input.time_scale_s;
          const auto rho = density_at(input, point, time_s);
          std::vector<double> eta(count);
          double grad = 0;
          for (unsigned int i = 0; i < count; ++i)
            {
              eta[i] = values[i][q];
              grad += gradient / 2 * gradients[i][q].norm_square();
            }
          const auto well = phasefield::growth::potential(eta, interaction);
          const auto local = evaluate(input, eta, tracking.local_grains(point, eta), rho);
          if (!save) continue;
          stream << timer.get_increment() << ',' << timer.get_time() << ',' << time_s;
          for (const auto x : point) stream << ',' << x;
          stream << ',' << fe.JxW(q) << ',' << rho << ',' << local.rho_effective << ',' << local.h_parent << ',' << local.h_rx
                 << ',' << local.energy * input.energy_scale << ',' << well * input.energy_scale
                 << ',' << grad * input.energy_scale << '\n';
        }
      for (auto &cell : cells) ++cell;
    }
  if (save && !stream) throw std::runtime_error("Stored-energy diagnostics write failed");
}
template void write_diagnostics<2>(prismspf::SolveContext<2, 1, double> &, const Input &, const grains::NativeTracking<2> &, double, double);
template void write_diagnostics<3>(prismspf::SolveContext<3, 1, double> &, const Input &, const grains::NativeTracking<3> &, double, double);
}

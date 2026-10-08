// SPDX-License-Identifier: LGPL-2.1-or-later
#include "field_access.h"
#include "../io.h"
#include "../remapping.h"

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace bravais::grains
{
template <unsigned int dim>
NativeTracking<dim>::NativeTracking(const InitialIdentity &input) : identity(input)
{
  if (input.grid.dimension != dim)
    throw std::invalid_argument("Grain identity dimension differs from native app");
  // Candidate absolute field-tail cutoff. Independent channels remain the
  // untruncated reference; reduced channels report every removal in output.
  if (input.count < input.channels.size()) support.cutoff = reused_field_tail_cutoff;
}

template <unsigned int dim>
void NativeTracking<dim>::bind_grid(prismspf::SolveContext<dim, 1, double> &context)
{
  const auto &input = context.get_user_inputs();
  if (context.get_field_attributes().size() != identity.owners().size())
    throw std::invalid_argument("Native field count differs from initial identity count");
  if (input.spatial_discretization.has_adaptivity)
    throw std::invalid_argument("Grain tracking requires the fixed Q1 input grid");
  const auto &extent = input.spatial_discretization.rectangular_mesh.size;
  for (unsigned int d = 0; d < dim; ++d) domain_extent[d] = extent[d];
  const auto nodes = identity.grid().node_count();
  const auto missing = dealii::numbers::invalid_dof_index;
  dofs.assign(identity.owners().size(), std::vector<dealii::types::global_dof_index>(nodes, missing));
  for (unsigned int channel = 0; channel < dofs.size(); ++channel)
    {
      if (context.get_constraint_manager().get_field_constraint(channel).n_constraints() != 0)
        throw std::invalid_argument("Grain tracking currently requires natural nonperiodic boundaries");
      const auto &handler = context.get_dof_manager().get_field_dof_handler(channel);
      if (handler.n_dofs() != nodes)
        throw std::invalid_argument("Native grid differs from initial grain grid");
      const auto points = dealii::DoFTools::map_dofs_to_support_points(dealii::MappingQ1<dim>(), handler);
      for (const auto &[dof, point] : points)
        {
          std::array<std::size_t, 3> coordinate{};
          for (unsigned int d = 0; d < dim; ++d)
            {
              const double scaled = snap_grid_coordinate(point[d] / extent[d] * identity.grid().cells[d],
                                                         identity.grid().cells[d]);
              const double rounded = std::round(scaled);
              if (!std::isfinite(scaled) || rounded < 0 || rounded > identity.grid().cells[d] ||
                  scaled != rounded)
                throw std::invalid_argument("Native Q1 node is not on the initial organization grid");
              coordinate[d] = static_cast<std::size_t>(rounded);
            }
          const auto node = coordinate[0] + (std::size_t(identity.grid().cells[0]) + 1) *
                            (coordinate[1] + (std::size_t(identity.grid().cells[1]) + 1) * coordinate[2]);
          if (dofs[channel][node] != missing)
            throw std::invalid_argument("Duplicate native Q1 grid node");
          dofs[channel][node] = dof;
        }
      for (const auto dof : dofs[channel])
        if (dof == missing) throw std::invalid_argument("Missing native Q1 grid node");
    }
}

template <unsigned int dim>
std::vector<unsigned int> NativeTracking<dim>::local_grains(
  const std::array<double, 3> &point, const std::span<const double> values) const
{
  if (dofs.empty() || values.size() != identity.owners().size())
    throw std::invalid_argument("Local grain query requires the bound native grid and field count");
  std::array<double, 3> coordinate{};
  for (unsigned int d = 0; d < dim; ++d)
    coordinate[d] = point[d] / domain_extent[d] * identity.grid().cells[d];
  std::vector<unsigned int> grains(values.size());
  for (unsigned int channel = 0; channel < values.size(); ++channel)
    grains[channel] = identity.grain_at(channel, coordinate, values[channel]);
  return grains;
}

template <unsigned int dim>
Fields NativeTracking<dim>::read_fields(prismspf::SolveContext<dim, 1, double> &context,
                                        const bool old) const
{
  Fields fields(dofs.size(), std::vector<double>(identity.grid().node_count()));
  for (unsigned int channel = 0; channel < dofs.size(); ++channel)
    {
      const auto &vector = old ? context.get_solution_indexer().get_old_solution_vector(0, channel)
                               : context.get_solution_indexer().get_solution_vector(channel);
      for (std::size_t node = 0; node < identity.grid().node_count(); ++node)
        fields[channel][node] = vector[dofs[channel][node]];
    }
  return fields;
}

template <unsigned int dim>
void NativeTracking<dim>::before_step(prismspf::SolveContext<dim, 1, double> &context)
{
  const auto cycle = context.get_simulation_timer().get_increment();
  if (cycle == 0) return;
  const auto fields = read_fields(context, true);
  identity.observe(fields);
  auto remapping = prepare_remapping(identity, fields, support);
  if (remapping.changed_grains != 0 || remapping.removed.count != 0)
    {
      for (unsigned int channel = 0; channel < dofs.size(); ++channel)
        {
          // The explicit RHS reads OldOne. The current buffer is scratch after
          // PRISMS-PF's end-of-increment swap; it must not be used as history.
          auto &vector = context.get_solution_indexer().get_old_solution_vector(0, channel);
          vector.zero_out_ghost_values();
          for (std::size_t node = 0; node < identity.grid().node_count(); ++node)
            vector[dofs[channel][node]] = remapping.fields[channel][node];
          vector.compress(dealii::VectorOperation::insert);
          vector.update_ghost_values();
        }
      std::cout << "Bravais grain remap cycle=" << cycle
                << " changed_grains=" << remapping.changed_grains
                << " removed_values=" << remapping.removed.count
                << " removed_max=" << remapping.removed.maximum << '\n';
    }
  removed.count += remapping.removed.count;
  removed.maximum = std::max(removed.maximum, remapping.removed.maximum);
  removed.absolute_sum += remapping.removed.absolute_sum;
  identity.replace(std::move(remapping.channels), std::move(remapping.predicted_owners));
}

template <unsigned int dim>
void NativeTracking<dim>::after_step(prismspf::SolveContext<dim, 1, double> &context)
{
  observe_step(context);
  write_step(context);
}

template <unsigned int dim>
void NativeTracking<dim>::observe_step(prismspf::SolveContext<dim, 1, double> &context)
{
  const auto &timer = context.get_simulation_timer();
  if (timer.get_increment() == 0) bind_grid(context);
  const auto fields = read_fields(context, false);
  identity.observe(fields);
}

template <unsigned int dim>
void NativeTracking<dim>::write_step(prismspf::SolveContext<dim, 1, double> &context) const
{
  const auto &timer = context.get_simulation_timer();
  const auto fields = read_fields(context, false);
  const auto &output = context.get_user_inputs().output_parameters;
  if (output.should_output(timer.get_increment()))
    write_identity(std::filesystem::path(output.folder) /
                     (output.file_name + "-grains-" + std::to_string(timer.get_increment()) + ".json"),
                   identity, fields, timer.get_increment(), timer.get_time(), support, removed);
}

template <unsigned int dim>
const Identity &NativeTracking<dim>::current_identity() const { return identity; }

template <unsigned int dim>
void NativeTracking<dim>::enable_channel_reuse() { support.cutoff = reused_field_tail_cutoff; }

template <unsigned int dim>
Fields NativeTracking<dim>::current_fields(prismspf::SolveContext<dim, 1, double> &context) const
{ return read_fields(context, false); }

template <unsigned int dim>
void NativeTracking<dim>::commit_current(prismspf::SolveContext<dim, 1, double> &context,
                                        Identity next, const Fields &fields)
{
  next.observe(fields);
  for (unsigned int channel = 0; channel < dofs.size(); ++channel)
    {
      auto &vector = context.get_solution_indexer().get_solution_vector(channel);
      vector.zero_out_ghost_values();
      for (std::size_t node = 0; node < identity.grid().node_count(); ++node)
        vector[dofs[channel][node]] = fields[channel][node];
      vector.compress(dealii::VectorOperation::insert);
      vector.update_ghost_values();
      // PRISMS-PF fills initial history before the cycle-zero post hook.
      if (context.get_simulation_timer().get_increment() == 0)
        context.get_solution_indexer().get_old_solution_vector(0, channel) = vector;
    }
  identity = std::move(next);
}

template class NativeTracking<2>;
template class NativeTracking<3>;
}

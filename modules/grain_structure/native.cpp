// SPDX-License-Identifier: LGPL-2.1-or-later
// v4 lifecycle/IO adaptation; grain algorithms in upstream/ retain their source.
#include "native.h"
#include <grains/FloodFiller.h>
#include <grains/OrderParameterRemapper.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/grid_tools.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>

namespace phasefield::grains
{
template <unsigned int dim>
NativeStructure<dim>::NativeStructure(const prismspf::UserInputParameters<dim> &input,
                                      unsigned int field_count) : count(field_count)
{
  const auto &c = input.user_constants;
  const auto cycles = c.get_int("grain_smoothing_cycles");
  const auto interval = c.get_int("grain_remap_interval");
  threshold = c.get_double("grain_threshold");
  buffer = c.get_double("grain_buffer");
  AssertThrow(count > 0 && cycles >= 0 && interval > 0 &&
              std::isfinite(threshold) && threshold > 0 && threshold < 1 &&
              std::isfinite(buffer) && buffer >= 0,
              dealii::ExcMessage("Invalid native grain controls"));
  smoothing_cycles = cycles;
  remap_interval = interval;
  AssertThrow(!input.spatial_discretization.has_adaptivity,
              dealii::ExcMessage("This v4 grain adapter requires a fixed Q1 grid"));
  std::ifstream file(c.get_string("grain_structure"));
  std::string header;
  unsigned int dimension = 0;
  file >> header >> dimension >> grid[0] >> grid[1] >> grid[2] >> grain_count;
  AssertThrow(file && header == "prisms_feature_grid_v1" && dimension == dim &&
              grain_count > 0 && grid[0] > 0 && grid[1] > 0 && grid[2] > 0 &&
              (dim == 3 || grid[2] == 1), dealii::ExcMessage("Invalid feature grid"));
  const std::size_t cells = std::size_t(grid[0]) * grid[1] * grid[2];
  AssertThrow(cells <= std::filesystem::file_size(c.get_string("grain_structure")) / 2 &&
              grain_count <= cells, dealii::ExcMessage("Feature grid exceeds payload"));
  features.resize(cells);
  std::set<unsigned int> ids;
  for (auto &id : features)
    {
      file >> id;
      AssertThrow(file && id > 0 && id <= grain_count, dealii::ExcMessage("Invalid FeatureId"));
      ids.insert(id);
    }
  std::string extra;
  AssertThrow(!(file >> extra) && ids.size() == grain_count,
              dealii::ExcMessage("Feature grid IDs must cover the orientation table"));
  for (unsigned int d = 0; d < dim; ++d)
    size[d] = input.spatial_discretization.rectangular_mesh.size[d];
  channels.assign(grain_count, 0);
}

template <unsigned int dim>
unsigned int NativeStructure<dim>::feature_at(const dealii::Point<dim> &point) const
{
  std::array<unsigned int, 3> p{};
  for (unsigned int d = 0; d < dim; ++d)
    {
      double coordinate = point[d] * grid[d] / size[d];
      const double integer = std::round(coordinate);
      // Four floating operations form the normalized coordinate; snap only
      // their machine-roundoff neighborhood, never a geometric tolerance.
      const double roundoff = 4 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(coordinate));
      if (std::abs(coordinate - integer) <= roundoff) coordinate = integer;
      p[d] = std::min(grid[d] - 1, static_cast<unsigned int>(std::max(0.0, std::floor(coordinate))));
    }
  return features[p[0] + grid[0] * (p[1] + grid[1] * p[2])];
}

template <unsigned int dim>
std::vector<typename NativeStructure<dim>::Vector *>
NativeStructure<dim>::fields(Context &context, bool old) const
{
  std::vector<Vector *> result;
  for (unsigned int i = 0; i < count; ++i)
    result.push_back(old ? &context.get_solution_indexer().get_old_solution_vector(0, i)
                         : &context.get_solution_indexer().get_solution_vector(i));
  return result;
}

template <unsigned int dim>
void NativeStructure<dim>::check_assignment(const std::vector<Grain> &items) const
{
  // Upstream's greedy routine may leave conflicts. Never silently accept a
  // grain merger or increase the requested count in the product adapter.
  for (std::size_t i = 0; i < items.size(); ++i)
    for (std::size_t j = 0; j < i; ++j)
      if (items[i].getOrderParameterId() == items[j].getOrderParameterId())
        AssertThrow(items[i].getCenter().distance(items[j].getCenter()) >=
                     items[i].getRadius() + items[j].getRadius() + 2 * buffer,
                     dealii::ExcMessage("Upstream grain reassignment left overlapping spheres for grain " +
                       std::to_string(items[j].getGrainId()) + " and " + std::to_string(items[i].getGrainId()) +
                       " in channel " + std::to_string(items[i].getOrderParameterId()) +
                       " with count=" + std::to_string(count) + "; this assignment did not separate the grains"));
}

template <unsigned int dim>
void NativeStructure<dim>::initialize(Context &context)
{
  // v4 exposes a const DoFHandler; the legacy algorithms only mutate cell
  // traversal flags, not topology/DoFs. Preserve and restore those flags.
  auto &handler = const_cast<dealii::DoFHandler<dim> &>(context.get_dof_manager().get_field_dof_handler(0));
  auto &fe = const_cast<dealii::FESystem<dim> &>(dynamic_cast<const dealii::FESystem<dim> &>(handler.get_fe()));
  AssertThrow(handler.get_triangulation().n_active_cells() == features.size(),
              dealii::ExcMessage("Feature and native mesh cell counts differ"));
  for (unsigned int i = 0; i < count; ++i)
    AssertThrow(context.get_constraint_manager().get_field_constraint(i).n_constraints() == 0,
                dealii::ExcMessage("Native grain input currently requires natural boundaries"));
  auto vectors = fields(context, false);
  Vector labels(*vectors[0]);
  const auto points = dealii::DoFTools::map_dofs_to_support_points(dealii::MappingQ1<dim>(), handler);
  for (const auto &[dof, point] : points) labels[dof] = feature_at(point);
  labels.update_ghost_values();
  for (auto *v : vectors) *v = 0;
  // FeatureIds are cell data. Build the initial grain geometry directly from
  // those cells: taking a modal nodal label first can erase one-cell grains.
  std::vector<GrainSet<dim>> sets(grain_count);
  for (unsigned int id = 1; id <= grain_count; ++id)
    {
      sets[id - 1].setGrainIndex(id);
      sets[id - 1].setOrderParameterIndex(0);
    }
  for (const auto &cell : handler.active_cell_iterators())
    {
      std::vector<dealii::Point<dim>> vertices;
      for (unsigned int v = 0; v < dealii::GeometryInfo<dim>::vertices_per_cell; ++v)
        vertices.push_back(cell->vertex(v));
      sets[feature_at(cell->center()) - 1].addVertexList(vertices);
    }
  for (const auto &set : sets) grains.emplace_back(set);
  std::vector<unsigned int> ops(count);
  std::iota(ops.begin(), ops.end(), 0);
  if (count >= grain_count)
    for (auto &grain : grains) grain.setOrderParameterId(grain.getGrainId() - 1);
  else
    SimplifiedGrainManipulator<dim>().reassignGrains(grains, buffer, ops);
  check_assignment(grains);
  OrderParameterRemapper<dim>().remap_from_index_field(grains, &labels, vectors, handler, fe.dofs_per_cell);
  for (const auto &g : grains) channels[g.getGrainId() - 1] = g.getOrderParameterId();
  smooth(context, vectors);
  for (unsigned int i = 0; i < count; ++i)
    context.get_solution_indexer().get_old_solution_vector(0, i) = *vectors[i];
  std::cout << "Native grain initialization: grains=" << grains.size() << " fields=" << count
            << " smoothing=" << smoothing_cycles << '\n';
}

template <unsigned int dim>
void NativeStructure<dim>::smooth(Context &context, const std::vector<Vector *> &vectors)
{
  const auto &handler = context.get_dof_manager().get_field_dof_handler(0);
  const auto &fe = handler.get_fe();
  dealii::QGaussLobatto<dim> q(2);
  dealii::FEValues<dim> values(fe, q, dealii::update_values | dealii::update_gradients | dealii::update_JxW_values);
  const auto dx = dealii::GridTools::minimal_cell_diameter(handler.get_triangulation()) / std::sqrt(double(dim));
  const double dt = 0.25 * dx * dx / dim;
  std::vector<dealii::types::global_dof_index> dofs(fe.dofs_per_cell);
  for (unsigned int cycle = 0; cycle < smoothing_cycles; ++cycle)
    for (auto *field : vectors)
      {
        Vector rhs(*field), mass(*field);
        rhs = 0; mass = 0;
        for (const auto &cell : handler.active_cell_iterators())
          {
            values.reinit(cell); cell->get_dof_indices(dofs);
            std::vector<dealii::Tensor<1, dim>> gradients(q.size());
            values.get_function_gradients(*field, gradients);
            for (unsigned int k = 0; k < q.size(); ++k)
              for (unsigned int i = 0; i < dofs.size(); ++i)
                {
                  rhs[dofs[i]] += values.shape_grad(i, k) * gradients[k] * values.JxW(k);
                  mass[dofs[i]] += values.shape_value(i, k) * values.JxW(k);
                }
          }
        for (dealii::types::global_dof_index i = 0; i < field->size(); ++i)
          (*field)[i] -= dt * rhs[i] / mass[i];
        field->update_ghost_values();
      }
}

template <unsigned int dim>
void NativeStructure<dim>::before_step(Context &context)
{
  // Independent fields have a permanent physical identity, including split
  // portions of one input grain. They never need channel reassignment.
  if (count >= grain_count) return;
  const auto cycle = context.get_simulation_timer().get_increment();
  if (cycle < 2 || (cycle - 1) % remap_interval != 0) return;
  identify(context, true);
  auto vectors = fields(context, true);
  auto &handler = const_cast<dealii::DoFHandler<dim> &>(context.get_dof_manager().get_field_dof_handler(0));
  std::vector<unsigned int> ops(count);
  std::iota(ops.begin(), ops.end(), 0);
  SimplifiedGrainManipulator<dim>().reassignGrains(grains, buffer, ops);
  check_assignment(grains);
  OrderParameterRemapper<dim>().remap(grains, vectors, handler, handler.get_fe().dofs_per_cell);
  for (auto *v : vectors) v->update_ghost_values();
  for (const auto &g : grains) channels[g.getGrainId() - 1] = g.getOrderParameterId();
  std::cout << "Native grain reassignment before cycle=" << cycle << " grains=" << grains.size() << '\n';
}

template <unsigned int dim>
void NativeStructure<dim>::identify(Context &context, bool old)
{
  auto vectors = fields(context, old);
  auto &handler = const_cast<dealii::DoFHandler<dim> &>(context.get_dof_manager().get_field_dof_handler(0));
  auto &fe = const_cast<dealii::FESystem<dim> &>(dynamic_cast<const dealii::FESystem<dim> &>(handler.get_fe()));
  FloodFiller<dim, 1> flood(fe, dealii::QGaussLobatto<dim>(2));
  std::vector<Grain> next;
  std::vector<GrainSet<dim>> components;
  std::vector<bool> flags;
  handler.get_triangulation().save_user_flags(flags);
  for (unsigned int op = 0; op < count; ++op)
    {
      std::vector<GrainSet<dim>> sets;
      flood.calcGrainSets(fe, handler, vectors[op], threshold, 1 + threshold, 0, op, sets);
      for (const auto &s : sets) { next.emplace_back(s); components.push_back(s); }
    }
  const_cast<dealii::Triangulation<dim> &>(handler.get_triangulation()).load_user_flags(flags);
  AssertThrow(!next.empty(), dealii::ExcMessage("No grains identified above threshold"));
  SimplifiedGrainManipulator<dim> manipulator;
  manipulator.transferGrainIds(grains, next);
  std::set<unsigned int> ids;
  for (const auto &g : next)
    {
      AssertThrow(ids.insert(g.getGrainId()).second && channels[g.getGrainId() - 1] == g.getOrderParameterId(),
                  dealii::ExcMessage("Upstream nearest-center identity is ambiguous or crosses channels; refusing guessed physical grain IDs"));
    }
  cell_owners.assign(count, std::vector<unsigned int>(features.size(), 0));
  for (std::size_t i = 0; i < next.size(); ++i)
    for (const auto &vertices : components[i].getVertexList())
      {
        const auto center = (vertices.front() + vertices.back()) / 2;
        std::array<unsigned int, 3> p{};
        for (unsigned int d = 0; d < dim; ++d)
          p[d] = std::min(grid[d] - 1, unsigned(center[d] / size[d] * grid[d]));
        const auto index = p[0] + grid[0] * (p[1] + grid[1] * p[2]);
        auto &id = cell_owners[next[i].getOrderParameterId()][index];
        AssertThrow(id == 0, dealii::ExcMessage("Multiple connected components own one channel cell"));
        id = next[i].getGrainId();
      }
  for (const auto &g : next) channels[g.getGrainId() - 1] = g.getOrderParameterId();
  grains = std::move(next);
}

template <unsigned int dim>
void NativeStructure<dim>::write_frame(Context &context)
{
  if (count < grain_count) identify(context, false);
  const auto &timer = context.get_simulation_timer();
  const auto &output = context.get_user_inputs().output_parameters;
  if (!output.should_output(timer.get_increment())) return;
  const auto &handler = context.get_dof_manager().get_field_dof_handler(0);
  dealii::FEValues<dim> sample(handler.get_fe(), dealii::QGauss<dim>(1), dealii::update_values);
  const auto vectors = fields(context, false);
  std::vector<unsigned int> labels(features.size(), 0);
  std::vector<GrainSet<dim>> independent(count >= grain_count ? grain_count : 0);
  for (unsigned int id = 1; id <= independent.size(); ++id)
    {
      independent[id - 1].setGrainIndex(id);
      independent[id - 1].setOrderParameterIndex(id - 1);
    }
  for (const auto &cell : handler.active_cell_iterators())
    {
      sample.reinit(cell);
      unsigned int best_op = 0;
      double best = -std::numeric_limits<double>::infinity();
      for (unsigned int op = 0; op < count; ++op)
        {
          std::vector<double> value(1);
          sample.get_function_values(*vectors[op], value);
          AssertThrow(std::isfinite(value[0]), dealii::ExcMessage("Nonfinite order parameter"));
          if (value[0] > best) { best = value[0]; best_op = op; }
        }
      std::array<unsigned int, 3> p{};
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = std::min(grid[d] - 1, unsigned(cell->center()[d] / size[d] * grid[d]));
      const auto index = p[0] + grid[0] * (p[1] + grid[1] * p[2]);
      // With one field per grain, disconnected portions retain that field's
      // identity. Flood-fill centroids and modal nodal values cannot override it.
      const auto owner = count >= grain_count ? best_op + 1 : cell_owners[best_op][index];
      AssertThrow(owner > 0 && owner <= grain_count && best > threshold,
                  dealii::ExcMessage("No current connected grain owns the dominant cell field"));
      labels[index] = owner;
      if (!independent.empty())
        {
          std::vector<dealii::Point<dim>> vertices;
          for (unsigned int v = 0; v < dealii::GeometryInfo<dim>::vertices_per_cell; ++v)
            vertices.push_back(cell->vertex(v));
          independent[owner - 1].addVertexList(std::move(vertices));
        }
    }
  if (!independent.empty())
    {
      grains.clear();
      for (const auto &grain : independent)
        if (!grain.getVertexList().empty()) grains.emplace_back(grain);
    }
  const auto path = std::filesystem::path(output.folder) /
    (output.file_name + "-grains-" + std::to_string(timer.get_increment()) + ".json");
  std::filesystem::create_directories(path.parent_path());
  std::ofstream file(path);
  file << std::setprecision(17) << "{\"schema_version\":1,\"cycle\":" << timer.get_increment()
       << ",\"time\":" << timer.get_time() << ",\"dimension\":" << dim << ",\"grid\":["
       << grid[0] << ',' << grid[1] << ',' << grid[2]
       << "],\"support\":\"prisms_grain_spheres\",\"channel_by_grain\":[";
  for (std::size_t i = 0; i < channels.size(); ++i) file << (i ? "," : "") << channels[i];
  file << "],\"grain_id\":[";
  for (std::size_t i = 0; i < labels.size(); ++i) file << (i ? "," : "") << labels[i];
  file << "],\"identification_threshold\":" << threshold << ",\"grains\":[";
  for (std::size_t i = 0; i < grains.size(); ++i)
    {
      const auto &g = grains[i];
      file << (i ? "," : "") << "{\"grain_id\":" << g.getGrainId() << ",\"channel\":"
           << g.getOrderParameterId() << ",\"radius\":" << g.getRadius() << ",\"center\":[";
      for (unsigned int d = 0; d < dim; ++d) file << (d ? "," : "") << g.getCenter()[d];
      file << "]}";
    }
  file << "]}\n";
  file.close();
  AssertThrow(file, dealii::ExcMessage("Cannot write native grain identity"));
}
template class NativeStructure<2>;
template class NativeStructure<3>;
}

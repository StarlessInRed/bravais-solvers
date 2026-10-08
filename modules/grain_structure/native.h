// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <prismspf/solvers/solve_context.h>
#include <grains/SimplifiedGrainRepresentation.h>
#include <array>
#include <filesystem>
#include <vector>

namespace phasefield::grains
{
// v4 integration of the v2.4.2 grainGrowth_dream3d grain algorithms.
// Feature IDs enter here; no channel assignment is accepted from the client.
template <unsigned int dim>
class NativeStructure
{
public:
  using Context = prismspf::SolveContext<dim, 1, double>;
  NativeStructure(const prismspf::UserInputParameters<dim> &, unsigned int count);
  void initialize(Context &);
  void before_step(Context &);
  void write_frame(Context &);

private:
  using Vector = dealii::LinearAlgebra::distributed::Vector<double>;
  using Grain = SimplifiedGrainRepresentation<dim>;
  std::vector<Vector *> fields(Context &, bool old) const;
  void smooth(Context &, const std::vector<Vector *> &);
  void identify(Context &, bool old);
  void check_assignment(const std::vector<Grain> &) const;
  unsigned int feature_at(const dealii::Point<dim> &) const;
  unsigned int count, smoothing_cycles, remap_interval, grain_count;
  double threshold, buffer;
  std::array<unsigned int, 3> grid{};
  std::array<double, dim> size{};
  std::vector<unsigned int> features, channels;
  std::vector<std::vector<unsigned int>> cell_owners;
  std::vector<Grain> grains;
};
}

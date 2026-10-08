// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include "../identity.h"
#include "../remapping.h"
#include <prismspf/solvers/solve_context.h>
#include <span>

namespace bravais::grains
{
template <unsigned int dim>
class NativeTracking
{
public:
  explicit NativeTracking(const InitialIdentity &input);
  void before_step(prismspf::SolveContext<dim, 1, double> &context);
  void after_step(prismspf::SolveContext<dim, 1, double> &context);
  void observe_step(prismspf::SolveContext<dim, 1, double> &context);
  void write_step(prismspf::SolveContext<dim, 1, double> &context) const;
  const Identity &current_identity() const;
  void enable_channel_reuse();
  Fields current_fields(prismspf::SolveContext<dim, 1, double> &context) const;
  void commit_current(prismspf::SolveContext<dim, 1, double> &context,
                      Identity next, const Fields &fields);
  std::vector<unsigned int> local_grains(const std::array<double, 3> &point,
                                       std::span<const double> values) const;

private:
  void bind_grid(prismspf::SolveContext<dim, 1, double> &context);
  Fields read_fields(prismspf::SolveContext<dim, 1, double> &context, bool old) const;
  Identity identity;
  SupportPolicy support;
  TailRemoval removed;
  std::vector<std::vector<dealii::types::global_dof_index>> dofs;
  std::array<double, 3> domain_extent{};
};

extern template class NativeTracking<2>;
extern template class NativeTracking<3>;
}

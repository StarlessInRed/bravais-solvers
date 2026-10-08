// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include "../model.h"
#include <grain_representation/prismspf/field_access.h>
#include <prismspf/solvers/solve_context.h>

namespace bravais::defects
{
// Native Q1 quadrature diagnostics, using the same energy evaluator as the RHS.
template <unsigned int dim>
void write_diagnostics(prismspf::SolveContext<dim, 1, double> &context,
                       const Input &input, const grains::NativeTracking<dim> &tracking,
                       double interaction, double gradient);
extern template void write_diagnostics<2>(prismspf::SolveContext<2, 1, double> &, const Input &, const grains::NativeTracking<2> &, double, double);
extern template void write_diagnostics<3>(prismspf::SolveContext<3, 1, double> &, const Input &, const grains::NativeTracking<3> &, double, double);
}

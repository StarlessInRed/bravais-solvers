// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include "identity.h"
#include <set>

namespace bravais::grains
{
struct SupportPrediction
{
  std::vector<std::set<std::size_t>> nodes_by_grain;
  std::vector<std::set<unsigned int>> conflicts;
};

// Exact support, without amplitude threshold: one explicit Q1 element stencil.
SupportPrediction predict_support(const Identity &identity);
}

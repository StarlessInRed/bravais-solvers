// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include "tracking.h"

namespace bravais::grains
{
inline constexpr double reused_field_tail_cutoff = 1e-12;
// Numerical implementation settings, not material parameters or user Q choices.
struct SupportPolicy
{
  double cutoff = 0;
  std::size_t search_budget = 100000;
};

struct TailRemoval
{
  std::size_t count = 0;
  double maximum = 0;
  double absolute_sum = 0;
};

struct Remapping
{
  std::vector<unsigned int> channels;
  Owners predicted_owners;
  Fields fields;
  unsigned int changed_grains;
  TailRemoval removed;
};

// Construct everything before mutating the native solution or identity.
Remapping prepare_remapping(const Identity &identity, const Fields &fields,
                           SupportPolicy policy = {});

// Bounded deterministic DSATUR search; exhaustion is not an impossibility proof.
std::vector<unsigned int> allocate_channels(const std::vector<std::set<unsigned int>> &conflicts,
                                          const std::vector<unsigned int> &preferred,
                                          unsigned int count, std::size_t search_budget);
}

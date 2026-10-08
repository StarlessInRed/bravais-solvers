// SPDX-License-Identifier: LGPL-2.1-or-later
#include "remapping.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>

namespace bravais::grains
{
std::vector<unsigned int> allocate_channels(const std::vector<std::set<unsigned int>> &conflicts,
                                          const std::vector<unsigned int> &preferred,
                                          const unsigned int count, std::size_t search_budget)
{
  const auto unassigned = std::numeric_limits<unsigned int>::max();
  if (count == 0 || conflicts.size() != preferred.size())
    throw std::invalid_argument("Invalid channel allocation dimensions");
  for (unsigned int grain = 0; grain < conflicts.size(); ++grain)
    {
      if (preferred[grain] >= count) throw std::invalid_argument("Preferred channel outside fixed Q");
      for (const auto other : conflicts[grain])
        if (other >= conflicts.size() || other == grain || !conflicts[other].contains(grain))
          throw std::invalid_argument("Channel conflict graph must be symmetric and loop-free");
    }
  std::vector<unsigned int> channels(preferred.size(), unassigned);
  struct Choice { unsigned int grain; std::vector<unsigned int> candidates; std::size_t next = 0; };
  std::vector<Choice> stack;
  // Iterative backtracking keeps the C++ call stack independent of grain count.
  while (true)
    {
      if (stack.size() == preferred.size()) return channels;
      unsigned int selected = unassigned;
      std::tuple<std::size_t, std::size_t> best{0, 0};
      for (unsigned int grain = 0; grain < preferred.size(); ++grain)
        if (channels[grain] == unassigned)
          {
            std::set<unsigned int> saturation;
            for (const auto other : conflicts[grain])
              if (channels[other] != unassigned) saturation.insert(channels[other]);
            const auto rank = std::make_tuple(saturation.size(), conflicts[grain].size());
            if (selected == unassigned || rank > best) { selected = grain; best = rank; }
          }
      std::set<unsigned int> forbidden;
      for (const auto other : conflicts[selected])
        if (channels[other] != unassigned) forbidden.insert(channels[other]);
      Choice choice{selected, {}};
      if (!forbidden.contains(preferred[selected])) choice.candidates.push_back(preferred[selected]);
      for (unsigned int channel = 0; channel < count; ++channel)
        if (channel != preferred[selected] && !forbidden.contains(channel))
          choice.candidates.push_back(channel);
      // Unused colours are equivalent for feasibility. Try only the preferred
      // first representative, rather than enumerating factorial permutations.
      const std::set<unsigned int> used(channels.begin(), channels.end());
      bool kept_unused = false;
      std::erase_if(choice.candidates, [&](const unsigned int channel) {
        if (used.contains(channel)) return false;
        if (kept_unused) return true;
        kept_unused = true;
        return false;
      });
      stack.push_back(std::move(choice));
      while (!stack.empty())
        {
          auto &last = stack.back();
          channels[last.grain] = unassigned;
          if (last.next < last.candidates.size())
            {
              if (search_budget == 0)
                throw std::runtime_error("Fixed Q allocation search budget exhausted; feasibility is unknown");
              --search_budget;
              channels[last.grain] = last.candidates[last.next++];
              break;
            }
          stack.pop_back();
        }
      if (stack.empty())
        throw std::runtime_error("Fixed Q has no coloring for this predicted support graph; Q was not increased");
    }
}

Remapping prepare_remapping(const Identity &identity, const Fields &fields, const SupportPolicy policy)
{
  if (!std::isfinite(policy.cutoff) || policy.cutoff < 0 || policy.cutoff >= 1 || policy.search_budget == 0)
    throw std::invalid_argument("Invalid grain support policy");
  auto active = identity;
  active.observe(fields);
  auto retained = fields;
  TailRemoval removed;
  for (auto &field : retained)
    for (auto &value : field)
      if (value != 0 && std::abs(value) <= policy.cutoff)
        {
          ++removed.count;
          removed.maximum = std::max(removed.maximum, std::abs(value));
          removed.absolute_sum += std::abs(value);
          value = 0;
        }
  active.observe(retained);
  const auto prediction = predict_support(active);
  const auto count = static_cast<unsigned int>(fields.size());
  const auto &old = identity.channels();
  std::vector<unsigned int> channels = old;
  bool valid = true;
  for (unsigned int grain = 0; grain < old.size(); ++grain)
    for (const auto other : prediction.conflicts[grain])
      if (old[grain] == old[other]) valid = false;
  if (!valid) channels = allocate_channels(prediction.conflicts, old, count, policy.search_budget);
  Owners owners(count, std::vector<unsigned int>(identity.grid().node_count()));
  for (unsigned int grain = 0; grain < old.size(); ++grain)
    for (const auto node : prediction.nodes_by_grain[grain])
      {
        auto &owner = owners[channels[grain]][node];
        if (owner != 0 && owner != grain + 1)
          throw std::logic_error("Conflicting remapping prediction");
        owner = grain + 1;
      }
  Fields moved(count, std::vector<double>(identity.grid().node_count()));
  for (unsigned int channel = 0; channel < count; ++channel)
    for (std::size_t node = 0; node < identity.grid().node_count(); ++node)
      {
        const auto grain = active.owners()[channel][node];
        if (grain != 0) moved[channels[grain - 1]][node] = retained[channel][node];
      }
  unsigned int changed = 0;
  for (unsigned int grain = 0; grain < old.size(); ++grain)
    if (channels[grain] != old[grain]) ++changed;
  return {std::move(channels), std::move(owners), std::move(moved), changed, removed};
}

}

// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <istream>
#include <memory>

namespace bravais::defects
{
// Isothermal branch evolution from immutable initial density, not an OP fraction.
class Recovery
{
public:
  virtual ~Recovery() = default;
  virtual double density(double initial_density, double time_s) const = 0;
};

class NoRecovery final : public Recovery
{
public:
  double density(double initial_density, double) const override { return initial_density; }
};

class FirstOrder final : public Recovery
{
public:
  FirstOrder(double rate_s_inv, double target_density_m2);
  double density(double initial_density, double time_s) const override;

private:
  double rate, target;
};

std::shared_ptr<const Recovery> read_recovery(std::istream &stream);
}

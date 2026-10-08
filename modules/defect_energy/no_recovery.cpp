// SPDX-License-Identifier: LGPL-2.1-or-later
#include "recovery.h"
#include <stdexcept>
#include <string>

namespace bravais::defects
{
std::shared_ptr<const Recovery> read_recovery(std::istream &stream)
{
  std::string type;
  stream >> type;
  if (type != "none") throw std::invalid_argument("This app uses fixed initial GND density (recovery=none)");
  return std::make_shared<NoRecovery>();
}
}

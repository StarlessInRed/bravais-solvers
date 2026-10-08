// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include "identity.h"
#include "remapping.h"
#include <filesystem>

namespace bravais::grains
{
InitialIdentity read_identity(const std::filesystem::path &path);
void write_identity(const std::filesystem::path &path, const Identity &identity,
                    const Fields &fields, unsigned int cycle, double time,
                    SupportPolicy policy = {}, TailRemoval removed = {});
}

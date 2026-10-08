// SPDX-License-Identifier: LGPL-2.1-or-later
#include "io.h"
#include <fstream>
#include <iomanip>
#include <charconv>
#include <stdexcept>

namespace bravais::grains
{
namespace
{
unsigned int read_unsigned(std::istream &file)
{
  std::string token;
  file >> token;
  unsigned int value = 0;
  const auto parsed = std::from_chars(token.data(), token.data() + token.size(), value);
  if (!file || parsed.ec != std::errc{} || parsed.ptr != token.data() + token.size())
    throw std::invalid_argument("Grain identity expects an unsigned decimal integer");
  return value;
}
}

InitialIdentity read_identity(const std::filesystem::path &path)
{
  std::ifstream file(path);
  InitialIdentity input{};
  std::string format;
  file >> format;
  input.grid.dimension = read_unsigned(file);
  for (auto &axis : input.grid.cells) axis = read_unsigned(file);
  const auto grains = read_unsigned(file);
  input.count = read_unsigned(file);
  if (!file || format != "bravais_grain_identity_v1" || grains == 0)
    throw std::invalid_argument("Invalid initial grain identity header");
  input.grid.validate();
  // Even one-character values need separators. Reject impossible/truncated
  // declared shapes before allocating arrays, without a machine-specific cap.
  const auto bytes = std::filesystem::file_size(path);
  if (grains > input.grid.cell_count() || input.grid.cell_count() > bytes / 2)
    throw std::invalid_argument("Initial identity shape exceeds its payload");
  input.channels.resize(grains);
  input.cell_grains.resize(input.grid.cell_count());
  for (auto &channel : input.channels) channel = read_unsigned(file);
  for (auto &grain : input.cell_grains) grain = read_unsigned(file);
  if (!file) throw std::invalid_argument("Incomplete initial grain identity");
  std::string extra;
  if (file >> extra) throw std::invalid_argument("Unexpected initial grain identity data");
  input.validate();
  return input;
}

void write_identity(const std::filesystem::path &path, const Identity &identity,
                    const Fields &fields, const unsigned int cycle, const double time,
                    const SupportPolicy policy, const TailRemoval removed)
{
  const auto grains = identity.cell_grains(fields);
  std::filesystem::create_directories(path.parent_path());
  std::ofstream file(path);
  file << std::setprecision(17) << "{\"schema_version\":1,\"cycle\":" << cycle
       << ",\"time\":" << time << ",\"dimension\":" << identity.grid().dimension
       << ",\"grid\":[" << identity.grid().cells[0] << ',' << identity.grid().cells[1]
       << ',' << identity.grid().cells[2] << "],\"support\":\"bounded_tail_q1\","
       << "\"tail_cutoff\":" << policy.cutoff
       << ",\"removed_value_count\":" << removed.count
       << ",\"removed_maximum\":" << removed.maximum
       << ",\"removed_absolute_sum\":" << removed.absolute_sum << ','
       << "\"channel_by_grain\":[";
  for (std::size_t i = 0; i < identity.channels().size(); ++i)
    file << (i ? "," : "") << identity.channels()[i];
  file << "],\"grain_id\":[";
  for (std::size_t i = 0; i < grains.size(); ++i) file << (i ? "," : "") << grains[i];
  file << "]}\n";
  file.close();
  if (!file) throw std::runtime_error("Cannot write frame grain identity");
}
}

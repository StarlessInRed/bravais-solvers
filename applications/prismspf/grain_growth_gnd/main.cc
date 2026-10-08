// SPDX-License-Identifier: LGPL-2.1-or-later
#include <defect_energy/prismspf/operator.h>

template <unsigned int dim>
class GrainGrowth : public StoredEnergyGrowth<dim>
{
public:
  GrainGrowth(const UserInputParameters<dim> &input, PhaseFieldTools<dim> &tools,
              const unsigned int count)
    : StoredEnergyGrowth<dim>(input, tools, count)
  {
    AssertThrow(this->stored_energy.has_value(), dealii::ExcMessage("GND recrystallization requires stored_energy"));
  }
private:
  // All nuclei already exist in the supplied initial state.
  void update_grains(SolveContext<dim, 1, double> &) override {}
};

template <unsigned int dim>
void run(const std::string &filename)
{
  UserInputParameters<dim> input(filename);
  const int count = input.user_constants.get_int("order_parameter_count");
  AssertThrow(count >= 1, dealii::ExcMessage("Expected at least one field"));
  std::vector<FieldAttributes> fields;
  std::set<std::string> dependencies;
  SolveBlock block;
  block.id = 0;
  block.solve_type = Explicit;
  block.solve_timing = Primary;
  for (int i = 0; i < count; ++i)
    {
      const auto name = "n" + std::to_string(i);
      fields.emplace_back(name, Scalar);
      block.field_indices.insert(i);
      dependencies.insert("old_1(" + name + ")");
      dependencies.insert("grad(old_1(" + name + "))");
    }
  block.dependencies_rhs = make_dependency_set(fields, dependencies);
  PhaseFieldTools<dim> tools;
  GrainGrowth<dim> equations(input, tools, count);
  std::cout << "Bravais PRISMS-PF v4 grain_growth_gnd dimension=" << dim << " degree=1 fields="
            << count << std::endl;
  Problem<dim, 1, double> problem(fields, {block}, input, tools, equations);
  problem.solve();
}

int main(int argc, char **argv)
{
  MPIInitFinalize mpi(argc, argv);
  try
    {
      AssertThrow(dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1,
                  dealii::ExcMessage("This application currently supports one MPI rank"));
      ParseCMDOptions options(argc, argv);
      const auto filename = options.get_parameters_filename();
      UserInputParameters<3> selection(filename);
      const auto dimension = selection.user_constants.get_int("dimension");
      AssertThrow(dimension == 2 || dimension == 3, dealii::ExcMessage("Expected dimension 2 or 3"));
      if (dimension == 2) run<2>(filename);
      else run<3>(filename);
    }
  catch (const std::exception &error)
    {
      std::cerr << "Bravais grain_growth_gnd failed: " << error.what() << '\n';
      // Upstream manual timer sections can remain active on an exception.
      // Clear them while MPI is alive, before static TimerOutput destruction.
      Timer::serial_timer().reset();
      return 1;
    }
  return 0;
}

// SPDX-License-Identifier: LGPL-2.1-or-later
// v4 runtime-count app, using PRISMS v2.4.2 grain algorithms (see README).
#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/pde_operator_base.h>
#include <prismspf/core/problem.h>
#include <prismspf/core/timer.h>
#include <grain_structure/native.h>
#include <grain_growth/model.h>
#include <iostream>
#include <set>

using namespace prismspf;
template <unsigned int dim>
class GrainGrowth : public PDEOperatorBase<dim, 1, double>
{
public:
  GrainGrowth(const UserInputParameters<dim> &input, PhaseFieldTools<dim> &tools, unsigned int n)
    : PDEOperatorBase<dim, 1, double>(input, tools), count(n), structure(input, n),
      mobility(input.user_constants.get_double("mobility")),
      gradient(input.user_constants.get_double("gradient_energy")),
      interaction(input.user_constants.get_double("interaction"))
  {
    AssertThrow(std::isfinite(mobility) && mobility > 0 && std::isfinite(gradient) && gradient > 0 &&
                std::isfinite(interaction) && interaction > 0, dealii::ExcMessage("Invalid growth constants"));
    const auto names = UserConstants<dim>::get_names(input.user_constants.file_name);
    AssertThrow(!names.contains("stored_energy") && !names.contains("birth"),
                dealii::ExcMessage("Basic app does not implement stored energy or birth"));
  }
private:
  void pre_solve_block(SolveContext<dim, 1, double> &context, unsigned int) override
  { structure.before_step(context); }
  void post_solve_block(SolveContext<dim, 1, double> &context, unsigned int) override
  {
    if (context.get_simulation_timer().get_increment() == 0) structure.initialize(context);
    structure.write_frame(context);
  }
  void set_initial_condition(const unsigned int &, const unsigned int &,
                            const dealii::Point<dim> &, double &value, double &) const override
  { value = 0; } // The cycle-zero hook assigns FeatureIds and smooths all channels together.
  void compute_rhs(FieldContainer<dim, 1, double> &fields,
                   const SimulationTimer &timer, unsigned int) const override
  {
    using Value = dealii::VectorizedArray<double>;
    Value sum = {};
    std::vector<Value> values(count);
    std::vector<dealii::Tensor<1, dim, Value>> gradients(count);
    for (unsigned int i = 0; i < count; ++i)
      { const auto n = fields.template get_value<Scalar, OldOne>(i); sum += n * n; }
    for (unsigned int i = 0; i < count; ++i)
      {
        const auto n = fields.template get_value<Scalar, OldOne>(i);
        Value derivative;
        for (unsigned int lane = 0; lane < Value::size(); ++lane)
          derivative[lane] = phasefield::growth::potential_derivative(n[lane], sum[lane], interaction);
        values[i] = n - timer.get_timestep() * mobility * derivative;
        gradients[i] = -timer.get_timestep() * mobility * gradient * fields.template get_gradient<Scalar, OldOne>(i);
      }
    for (unsigned int i = 0; i < count; ++i)
      { fields.set_value_term(i, values[i]); fields.set_gradient_term(i, gradients[i]); }
  }
  unsigned int count;
  phasefield::grains::NativeStructure<dim> structure;
  double mobility, gradient, interaction;
};

template <unsigned int dim>
void run(const std::string &filename)
{
  UserInputParameters<dim> input(filename);
  const int count = input.user_constants.get_int("order_parameter_count");
  AssertThrow(count > 0, dealii::ExcMessage("Expected positive order parameter count"));
  std::vector<FieldAttributes> fields;
  std::set<std::string> dependencies;
  SolveBlock block; block.id = 0; block.solve_type = Explicit; block.solve_timing = Primary;
  for (int i = 0; i < count; ++i)
    {
      const auto name = "n" + std::to_string(i);
      fields.emplace_back(name, Scalar); block.field_indices.insert(i);
      dependencies.insert("old_1(" + name + ")"); dependencies.insert("grad(old_1(" + name + "))");
    }
  block.dependencies_rhs = make_dependency_set(fields, dependencies);
  PhaseFieldTools<dim> tools;
  GrainGrowth<dim> equations(input, tools, count);
  std::cout << "PRISMS-PF v4 grain_growth_dream3d dimension=" << dim << " degree=1 fields=" << count << '\n';
  Problem<dim, 1, double> problem(fields, {block}, input, tools, equations);
  problem.solve();
}
int main(int argc, char **argv)
{
  MPIInitFinalize mpi(argc, argv);
  try
    {
      AssertThrow(dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1,
                  dealii::ExcMessage("The v4 grain adapter supports one MPI rank"));
      ParseCMDOptions options(argc, argv);
      const auto filename = options.get_parameters_filename();
      UserInputParameters<3> selection(filename);
      const auto dim = selection.user_constants.get_int("dimension");
      AssertThrow(dim == 2 || dim == 3, dealii::ExcMessage("Expected dimension 2 or 3"));
      if (dim == 2) run<2>(filename); else run<3>(filename);
    }
  catch (const std::exception &error)
    { std::cerr << "Native grain growth failed: " << error.what() << '\n'; Timer::serial_timer().reset(); return 1; }
  return 0;
}

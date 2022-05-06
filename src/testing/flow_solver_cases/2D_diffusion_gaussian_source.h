#ifndef __2D_DIFFUSION_GAUSSIAN_SOURCE__
#define __2D_DIFFUSION_GAUSSIAN_SOURCE__


// for FlowSolver class:
#include "physics/initial_conditions/initial_condition.h"
#include "dg/dg.h"
#include "physics/physics.h"
#include "parameters/all_parameters.h"
#include <deal.II/base/table_handler.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>
#include "flow_solver_case_base.h"

namespace PHiLiP {
namespace Tests {

template <int dim, int nstate>
class DiffusionGaussian : public FlowSolverCaseBase<dim,nstate>
{
#if PHILIP_DIM==1
    using Triangulation = dealii::Triangulation<PHILIP_DIM>;
#else
    using Triangulation = dealii::parallel::distributed::Triangulation<PHILIP_DIM>;
#endif
public:
    /// Constructor.
    DiffusionGaussian(const Parameters::AllParameters *const parameters_input);

    /// Destructor
    ~DiffusionGaussian() {};

protected:
    /// Displays the flow setup parameters
    void display_flow_solver_setup(std::shared_ptr<InitialConditionFunction<dim,nstate,double>> initial_condition) const override;

    /// Function to generate the grid
    std::shared_ptr<Triangulation> generate_grid() const override;

};

} // Tests namespace
} // PHiLiP namespace

#endif

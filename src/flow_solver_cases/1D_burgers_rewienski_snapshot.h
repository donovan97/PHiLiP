#ifndef __1D_BURGERS_REWIENSKI_SNAPSHOT__
#define __1D_BURGERS_REWIENSKI_SNAPSHOT__

// for FlowSolver class:
#include "physics/initial_conditions/initial_condition.h"
#include "dg/dg.h"
#include "physics/physics.h"
#include "parameters/all_parameters.h"
#include <deal.II/base/table_handler.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>
#include "flow_solver_case_base.h"
#include <deal.II/base/function.h>
#include <stdlib.h>
#include <iostream>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_values.h>
#include "physics/physics_factory.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/vector.h>
#include "linear_solver/linear_solver.h"

namespace PHiLiP {
namespace FlowSolverCases {

#if PHILIP_DIM==1
    using Triangulation = dealii::Triangulation<PHILIP_DIM>;
#else
    using Triangulation = dealii::parallel::distributed::Triangulation<PHILIP_DIM>;
#endif

template <int dim, int nstate>
class BurgersRewienskiSnapshot: public FlowSolverCaseBase<dim, nstate>
{
public:
    /// Constructor.
    BurgersRewienskiSnapshot(const Parameters::AllParameters *const parameters_input);

    /// Destructor
    ~BurgersRewienskiSnapshot() {};

protected:
    const int number_of_refinements; ///< Number of cells per direction for the grid
    const double domain_left; ///< Domain left-boundary value for generating the grid
    const double domain_right; ///< Domain right-boundary value for generating the grid

    /// Function to generate the grid
    std::shared_ptr<Triangulation> generate_grid() const override;

    /// Function to write unsteady snapshot data to table
    void compute_unsteady_data_and_write_to_table(
            const unsigned int current_iteration,
            const double current_time,
            const std::shared_ptr <DGBase<dim, double>> dg,
            const std::shared_ptr<dealii::TableHandler> unsteady_data_table) const override;

    /// Function for postprocessing when solving for steady state
    void steady_state_postprocessing(std::shared_ptr <DGBase<dim, double>> dg) const override;

};

} // Tests namespace
} // PHiLiP namespace

#endif

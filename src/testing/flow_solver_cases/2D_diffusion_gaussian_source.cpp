#include "2D_diffusion_gaussian_source.h"
#include <deal.II/base/function.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <stdlib.h>
#include <iostream>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_values.h>
#include "physics/physics_factory.h"
#include "dg/dg.h"
#include <deal.II/base/table_handler.h>
#include "mesh/grids/naca_airfoil_grid.hpp"
#include "mesh/gmsh_reader.hpp"
#include "functional/lift_drag.hpp"

namespace PHiLiP {
namespace Tests {
//=========================================================
// NACA0012
//=========================================================
template <int dim, int nstate>
DiffusionGaussian<dim, nstate>::DiffusionGaussian(const PHiLiP::Parameters::AllParameters *const parameters_input)
    : FlowSolverCaseBase<dim, nstate>(parameters_input)
{
}

template <int dim, int nstate>
void DiffusionGaussian<dim,nstate>::display_flow_solver_setup(std::shared_ptr<InitialConditionFunction<dim,nstate,double>> /*initial_condition*/) const
{
    using PDE_enum = Parameters::AllParameters::PartialDifferentialEquation;
    const PDE_enum pde_type = this->all_param.pde_type;
    std::string pde_string;
    if (pde_type == PDE_enum::diffusion)                {pde_string = "diffusion";}
    this->pcout << "- PDE Type: " << pde_string << std::endl;
    this->pcout << "- Polynomial degree: " << this->all_param.grid_refinement_study_param.poly_degree << std::endl;
    this->pcout << "- Diffusion coefficient: " << this->all_param.manufactured_convergence_study_param.manufactured_solution_param.diffusion_coefficient << std::endl;
    this->pcout << "- Gaussian source height: " << this->all_param.manufactured_convergence_study_param.manufactured_solution_param.gaussian_source_height << std::endl;
}

template <int dim, int nstate>
std::shared_ptr<Triangulation> DiffusionGaussian<dim,nstate>::generate_grid() const
{
    std::shared_ptr<Triangulation> grid = std::make_shared<Triangulation>(
#if PHILIP_DIM!=1
            this->mpi_communicator
#endif
    );
    dealii::GridGenerator::subdivided_hyper_cube(*grid, this->all_param.grid_refinement_study_param.grid_size);
    return grid;
}

#if PHILIP_DIM==2
template class DiffusionGaussian<PHILIP_DIM,1>;
#endif

} // Tests namespace
} // PHiLiP namespace



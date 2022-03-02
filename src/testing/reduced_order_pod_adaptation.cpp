#include <fstream>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include "reduced_order_pod_adaptation.h"
#include "parameters/all_parameters.h"
#include "dg/dg_factory.hpp"
#include "ode_solver/ode_solver_factory.h"
#include "reduced_order/pod_adaptation.h"
#include "reduced_order/pod_basis_sensitivity.h"
#include "reduced_order/pod_basis_sensitivity_types.h"
#include "flow_solver.h"


namespace PHiLiP {
namespace Tests {

template <int dim, int nstate>
ReducedOrderPODAdaptation<dim, nstate>::ReducedOrderPODAdaptation(const PHiLiP::Parameters::AllParameters *const parameters_input)
        : TestsBase::TestsBase(parameters_input)
{}

template <int dim, int nstate>
int ReducedOrderPODAdaptation<dim, nstate>::run_test() const
{
    /*
    const Parameters::AllParameters param = *(TestsBase::all_parameters);

    pcout << "Running Burgers Rewienski with parameter a: "
          << param.burgers_param.rewienski_a
          << " and parameter b: "
          << param.burgers_param.rewienski_b
          << std::endl;

    std::shared_ptr<dealii::Triangulation<dim>> grid = std::make_shared<dealii::Triangulation<dim>>();

    double left = param.grid_refinement_study_param.grid_left;
    double right = param.grid_refinement_study_param.grid_right;
    const bool colorize = true;
    int n_refinements = param.grid_refinement_study_param.num_refinements;
    unsigned int poly_degree = param.grid_refinement_study_param.poly_degree;
    dealii::GridGenerator::hyper_cube(*grid, left, right, colorize);

    grid->refine_global(n_refinements);
    pcout << "Grid generated and refined" << std::endl;

    std::shared_ptr < PHiLiP::DGBase<dim, double> > dg = PHiLiP::DGFactory<dim,double>::create_discontinuous_galerkin(all_parameters, poly_degree, grid);
    pcout << "dg created" <<std::endl;
    dg->allocate_system ();

    // casting to dg state
    std::shared_ptr< DGBaseState<dim,nstate,double> > dg_state = std::dynamic_pointer_cast< DGBaseState<dim,nstate,double> >(dg);

    pcout << "Implement initial conditions" << std::endl;
    dealii::FunctionParser<1> initial_condition;
    std::string variables = "x";
    std::map<std::string,double> constants;
    constants["pi"] = dealii::numbers::PI;
    std::string expression = "1";
    initial_condition.initialize(variables, expression, constants);
    dealii::VectorTools::interpolate(dg->dof_handler,initial_condition,dg->solution);

    pcout << "Dimension: " << dim
          << "\t Polynomial degree p: " << poly_degree
          << std::endl
          << ". Number of active cells: " << grid->n_global_active_cells()
          << ". Number of degrees of freedom: " << dg->dof_handler.n_dofs()
          << std::endl;

    // Create functional
    auto burgers_functional = BurgersRewienskiFunctional<dim,nstate,double>(dg,dg_state->pde_physics_fad_fad,true,false);

    //POD adaptation
    std::shared_ptr<ProperOrthogonalDecomposition::PODAdaptation<dim, nstate>> pod_adapt = std::make_shared<ProperOrthogonalDecomposition::PODAdaptation<dim, nstate>>(dg, burgers_functional);
    pod_adapt->progressivePODAdaptation();

    //Evaluate functional on fine space to compare
    std::shared_ptr < PHiLiP::DGBase<dim, double> > dg_fine = PHiLiP::DGFactory<dim,double>::create_discontinuous_galerkin(all_parameters, poly_degree, grid);
    std::shared_ptr< DGBaseState<dim,nstate,double> > dg_state_fine = std::dynamic_pointer_cast< DGBaseState<dim,nstate,double> >(dg_fine);
    dg_fine->allocate_system ();
    dealii::VectorTools::interpolate(dg_fine->dof_handler,initial_condition,dg_fine->solution);
    //std::shared_ptr<ProperOrthogonalDecomposition::FineStatePOD<dim>> finePOD = std::make_shared<ProperOrthogonalDecomposition::FineStatePOD<dim>>(dg_fine);
    //std::shared_ptr<PHiLiP::ODE::ODESolverBase<dim, double>> ode_solver_fine = ODE::ODESolverFactory<dim, double>::create_ODESolver(dg_fine, finePOD);
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::implicit_solver;
    std::shared_ptr<PHiLiP::ODE::ODESolverBase<dim, double>> ode_solver_fine =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, dg_fine);
    ode_solver_fine->steady_state();
    auto functional_fine = BurgersRewienskiFunctional<dim,nstate,double>(dg_fine,dg_state_fine->pde_physics_fad_fad,true,false);

    pcout << "Fine functional: " << std::setprecision(15)  << functional_fine.evaluate_functional(false,false) << std::setprecision(6) << std::endl;
    pcout << "Coarse functional: " << std::setprecision(15)  << pod_adapt->getCoarseFunctional() << std::setprecision(6) << std::endl;

    if(abs(pod_adapt->getCoarseFunctional() - functional_fine.evaluate_functional(false,false)) > all_parameters->reduced_order_param.adaptation_tolerance){
        pcout << "Adaptation tolerance not reached." << std::endl;
        return -1;
    }
    else{
        pcout << "Adaptation tolerance reached." << std::endl;
        return 0;
    }
    */
    /*

    std::unique_ptr<FlowSolver<dim,nstate>> flow_solver_implicit = FlowSolverFactory<dim,nstate>::create_FlowSolver(all_parameters);
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::implicit_solver;
    flow_solver_implicit->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver_implicit->dg);
    flow_solver_implicit->ode_solver->allocate_ode_system();

    std::unique_ptr<FlowSolver<dim,nstate>> flow_solver_standard = FlowSolverFactory<dim,nstate>::create_FlowSolver(all_parameters);
    ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::pod_galerkin_solver;
    std::shared_ptr<ProperOrthogonalDecomposition::CoarseStatePOD<dim>> pod_standard = std::make_shared<ProperOrthogonalDecomposition::CoarseStatePOD<dim>>(flow_solver_standard->dg);
    flow_solver_standard->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver_standard->dg, pod_standard);
    flow_solver_standard->ode_solver->allocate_ode_system();

    std::unique_ptr<FlowSolver<dim,nstate>> flow_solver_expanded = FlowSolverFactory<dim,nstate>::create_FlowSolver(all_parameters);
    ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::pod_galerkin_solver;
    std::shared_ptr<ProperOrthogonalDecomposition::CoarseExpandedPOD<dim>> pod_expanded = std::make_shared<ProperOrthogonalDecomposition::CoarseExpandedPOD<dim>>(flow_solver_expanded->dg);
    flow_solver_expanded->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver_expanded->dg, pod_expanded);
    flow_solver_expanded->ode_solver->allocate_ode_system();

    std::unique_ptr<FlowSolver<dim,nstate>> flow_solver_extrapolated = FlowSolverFactory<dim,nstate>::create_FlowSolver(all_parameters);
    ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::pod_galerkin_solver;
    std::shared_ptr<ProperOrthogonalDecomposition::ExtrapolatedPOD<dim>> pod_extrapolated = std::make_shared<ProperOrthogonalDecomposition::ExtrapolatedPOD<dim>>(flow_solver_extrapolated->dg);
    flow_solver_extrapolated->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver_extrapolated->dg, pod_extrapolated);
    flow_solver_extrapolated->ode_solver->allocate_ode_system();
    */
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    /*Time-averaged relative error, E = 1/n_t * sum_{n=1}^{n_t} (||U_FOM(t^{n}) - U_ROM(t^{n})||_L2 / ||U_FOM(t^{n})||_L2 )
     *Refer to section 6.1 in "The GNAT method for nonlinear model reduction: Effective implementation and application to computational ﬂuid dynamics and turbulent ﬂows"
     *Authors: Kevin Carlberg, Charbel Farhat, ulien Cortial,  David Amsallem
     *Journal of Computational Physics, 2013
     */
    /*
    double finalTime = all_parameters->flow_solver_param.final_time;

    const unsigned int number_of_time_steps = static_cast<int>(ceil(finalTime/all_parameters->ode_solver_param.initial_time_step));
    const double constant_time_step = finalTime/number_of_time_steps;

    pcout << " Advancing solution by " << finalTime << " time units, using "
          << number_of_time_steps << " iterations of size dt=" << constant_time_step << " ... " << std::endl;

    double standard_error_norm_sum = 0;
    double expanded_error_norm_sum = 0;
    double extrapolated_error_norm_sum = 0;

    unsigned int current_iteration = 0;

    std::ofstream out_file("expandedbasis.txt");
    unsigned int precision = 7;
    pod_expanded->getPODBasis()->print(out_file, precision);

    std::ofstream out_file2("standardbasis.txt");
    pod_standard->getPODBasis()->print(out_file2, precision);


    while (current_iteration < number_of_time_steps)
    {
        pcout << " ********************************************************** "
              << std::endl
              << " Iteration: " << current_iteration + 1
              << " out of: " << number_of_time_steps
              << std::endl;

        const bool pseudotime = false;
        flow_solver_implicit->ode_solver->step_in_time(constant_time_step, pseudotime);
        flow_solver_standard->ode_solver->step_in_time(constant_time_step, pseudotime);
        flow_solver_expanded->ode_solver->step_in_time(constant_time_step, pseudotime);
        flow_solver_extrapolated->ode_solver->step_in_time(constant_time_step, pseudotime);

        dealii::LinearAlgebra::distributed::Vector<double> standard_solution(flow_solver_standard->dg->solution);
        dealii::LinearAlgebra::distributed::Vector<double> expanded_solution(flow_solver_expanded->dg->solution);
        dealii::LinearAlgebra::distributed::Vector<double> extrapolated_solution(flow_solver_extrapolated->dg->solution);
        dealii::LinearAlgebra::distributed::Vector<double> implicit_solution(flow_solver_implicit->dg->solution);

        standard_error_norm_sum = standard_error_norm_sum + ((standard_solution-=implicit_solution).l2_norm()/implicit_solution.l2_norm());
        expanded_error_norm_sum = expanded_error_norm_sum + (((expanded_solution-=implicit_solution).l2_norm())/implicit_solution.l2_norm());
        extrapolated_error_norm_sum = extrapolated_error_norm_sum + (((extrapolated_solution-=implicit_solution).l2_norm())/implicit_solution.l2_norm());

        pcout << (double)((standard_solution).l2_norm()/implicit_solution.l2_norm()) << std::endl;
        pcout << (double)((expanded_solution).l2_norm()/implicit_solution.l2_norm()) << std::endl;
        pcout << (double)((extrapolated_solution).l2_norm()/implicit_solution.l2_norm()) << std::endl;
        current_iteration++;
    }

    double standard_error = (1/(double)number_of_time_steps) * standard_error_norm_sum;
    double expanded_error = (1/(double)number_of_time_steps) * expanded_error_norm_sum;
    double extrapolated_error = (1/(double)number_of_time_steps) * extrapolated_error_norm_sum;

    pcout << "Standard error: " << standard_error << std::endl;
    pcout << "Expanded error: " << expanded_error << std::endl;
    pcout << "Extrapolated error: " << extrapolated_error << std::endl;

    return 0;
     */

    std::unique_ptr<FlowSolver<dim,nstate>> flow_solver_implicit = FlowSolverFactory<dim,nstate>::create_FlowSolver(all_parameters);
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::implicit_solver;
    flow_solver_implicit->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver_implicit->dg);
    flow_solver_implicit->ode_solver->allocate_ode_system();
    std::shared_ptr<DGBaseState<dim,nstate,double>> dg_state_implicit = std::dynamic_pointer_cast<DGBaseState<dim,nstate,double>>(flow_solver_implicit->dg);
    auto functional_implicit = BurgersRewienskiFunctional<dim,nstate,double>(flow_solver_implicit->dg, dg_state_implicit->pde_physics_fad_fad, true, false);

    std::unique_ptr<FlowSolver<dim,nstate>> flow_solver_standard = FlowSolverFactory<dim,nstate>::create_FlowSolver(all_parameters);
    ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::pod_petrov_galerkin_solver;
    std::shared_ptr<ProperOrthogonalDecomposition::CoarseStatePOD<dim>> pod_standard = std::make_shared<ProperOrthogonalDecomposition::CoarseStatePOD<dim>>(flow_solver_standard->dg);
    flow_solver_standard->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver_standard->dg, pod_standard);
    flow_solver_standard->ode_solver->allocate_ode_system();
    std::shared_ptr<DGBaseState<dim,nstate,double> > dg_state_standard = std::dynamic_pointer_cast< DGBaseState<dim,nstate,double>>(flow_solver_standard->dg);
    auto functional_standard = BurgersRewienskiFunctional<dim,nstate,double>(flow_solver_standard->dg, dg_state_standard->pde_physics_fad_fad, true, false);

    std::unique_ptr<FlowSolver<dim,nstate>> flow_solver_expanded = FlowSolverFactory<dim,nstate>::create_FlowSolver(all_parameters);
    ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::pod_petrov_galerkin_solver;
    std::shared_ptr<ProperOrthogonalDecomposition::CoarseExpandedPOD<dim>> pod_expanded = std::make_shared<ProperOrthogonalDecomposition::CoarseExpandedPOD<dim>>(flow_solver_expanded->dg);
    flow_solver_expanded->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver_expanded->dg, pod_expanded);
    flow_solver_expanded->ode_solver->allocate_ode_system();
    std::shared_ptr<DGBaseState<dim,nstate,double> > dg_state_expanded = std::dynamic_pointer_cast< DGBaseState<dim,nstate,double>>(flow_solver_expanded->dg);
    auto functional_expanded = BurgersRewienskiFunctional<dim,nstate,double>(flow_solver_expanded->dg, dg_state_expanded->pde_physics_fad_fad, true, false);

    flow_solver_implicit->ode_solver->steady_state();
    flow_solver_standard->ode_solver->steady_state();
    flow_solver_expanded->ode_solver->steady_state();

    dealii::LinearAlgebra::distributed::Vector<double> standard_solution(flow_solver_standard->dg->solution);
    dealii::LinearAlgebra::distributed::Vector<double> expanded_solution(flow_solver_expanded->dg->solution);
    dealii::LinearAlgebra::distributed::Vector<double> implicit_solution(flow_solver_implicit->dg->solution);

    double standard_error = ((standard_solution-=implicit_solution).l2_norm()/implicit_solution.l2_norm());
    double expanded_error = (((expanded_solution-=implicit_solution).l2_norm())/implicit_solution.l2_norm());
    double standard_func_error = abs(functional_implicit.evaluate_functional(false,false) - functional_standard.evaluate_functional(false,false));
    double expanded_func_error = abs(functional_implicit.evaluate_functional(false,false) - functional_expanded.evaluate_functional(false,false));

    pcout << "Standard error: " << standard_error << std::endl;
    pcout << "Expanded error: " << expanded_error << std::endl;
    pcout << "Standard func error: " << std::setprecision(15)  << standard_func_error << std::setprecision(6) << std::endl;
    pcout << "Expanded func error: " << std::setprecision(15)  << expanded_func_error << std::setprecision(6) << std::endl;

    return 0;
}

template <int dim, int nstate, typename real>
template <typename real2>
real2 BurgersRewienskiFunctional<dim,nstate,real>::evaluate_volume_integrand(
const PHiLiP::Physics::PhysicsBase<dim,nstate,real2> &/*physics*/,
const dealii::Point<dim,real2> &/*phys_coord*/,
const std::array<real2,nstate> &soln_at_q,
const std::array<dealii::Tensor<1,dim,real2>,nstate> &/*soln_grad_at_q*/) const
{
real2 val = 0;

// integrating over the domain
for (int istate=0; istate<nstate; ++istate) {
    val += soln_at_q[istate];
}

return val;
}


#if PHILIP_DIM==1
template class ReducedOrderPODAdaptation<PHILIP_DIM,PHILIP_DIM>;
template class BurgersRewienskiFunctional<PHILIP_DIM, PHILIP_DIM, double>;
#endif
} // Tests namespace
} // PHiLiP namespace

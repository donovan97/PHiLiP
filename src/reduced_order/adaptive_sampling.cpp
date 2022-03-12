#include "adaptive_sampling.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template <int dim, int nstate>
AdaptiveSampling<dim, nstate>::AdaptiveSampling(std::vector<double> parameter_space, double tolerance, const PHiLiP::Parameters::AllParameters *const parameters_input)
        : all_parameters(parameters_input)
        , parameter_space(parameter_space)
        , tolerance(tolerance)
{
}

template <int dim, int nstate>
void AdaptiveSampling<dim, nstate>::start(){
    std::cout << "Starting adaptive sampling process" << std::endl;

    //Generate trial locations
    int n = 10;
    generateTrialLocations(n);
    /*
    for(int i = 0 ; i < n ; i++){
        std::cout << trial_locations.front() << std::endl;
        trial_locations.pop();
    }
    */
    initializeSampling();


}

template <int dim, int nstate>
void AdaptiveSampling<dim, nstate>::generateTrialLocations(int n){
    std::cout << "Generating trial locations" << std::endl;
    for(int m = 1 ; m <= n ; m++){
        int base = 2;
        double q=0;
        double bk=(double)1/base;
        int i = m;

        while (i > 0) {
          q += (i % base)*bk;
          i /= base;
          bk /= base;
        }

        trial_locations.emplace(q);
    }
}

template <int dim, int nstate>
void AdaptiveSampling<dim, nstate>::initializeSampling(){
    int initialSnapshots = 2;
    for(int idx = 0 ; idx < initialSnapshots ; idx++){
        std::cout << idx << std::endl;
        std::shared_ptr<FOMSolution<dim,nstate>> fom_solution = solveSnapshotFOM(trial_locations.front());
        Snapshot<dim,nstate> snapshot = Snapshot<dim,nstate>(trial_locations.front(), fom_solution);
        std::shared_ptr<ROMSolution<dim,nstate>> rom_solution = solveSnapshotROM(trial_locations.front());
        snapshot.add_ROM(rom_solution);
        snapshots.push_back(snapshot);
        trial_locations.pop();
    }
}

template <int dim, int nstate>
std::shared_ptr<FOMSolution<dim,nstate>> AdaptiveSampling<dim, nstate>::solveSnapshotFOM(double parameter){

    Parameters::AllParameters params = reinitParams(parameter);

    std::shared_ptr<Tests::BurgersRewienskiSnapshot<dim, nstate>> flow_solver_case = std::make_shared<Tests::BurgersRewienskiSnapshot<dim,nstate>>(all_parameters);
    std::unique_ptr<Tests::FlowSolver<dim,nstate>> flow_solver = std::make_unique<Tests::FlowSolver<dim,nstate>>(all_parameters, flow_solver_case);

    // Solve implicit solution
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::implicit_solver;
    flow_solver->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver->dg);
    flow_solver->ode_solver->steady_state();

    // Casting to dg state
    std::shared_ptr< DGBaseState<dim,nstate,double>> dg_state = std::dynamic_pointer_cast< DGBaseState<dim,nstate,double> >(flow_solver->dg);

    // Create functional
    auto functional = BurgersRewienskiFunctional<dim,nstate,double>(flow_solver->dg,dg_state->pde_physics_fad_fad,true,false);

    //Get sensitivity from FlowSolver
    //std::shared_ptr<Tests::BurgersRewienskiSnapshot<dim,nstate>> burgers_rewienski_snapshot = std::dynamic_pointer_cast<Tests::BurgersRewienskiSnapshot<dim,nstate>>(flow_solver->flow_solver_case);

    std::shared_ptr<FOMSolution<dim,nstate>> fom_solution = std::make_shared<FOMSolution<dim, nstate>>(flow_solver->dg, functional, *flow_solver_case->sensitivity_dWdb);

    return fom_solution;
}

template <int dim, int nstate>
std::shared_ptr<ROMSolution<dim,nstate>> AdaptiveSampling<dim, nstate>::solveSnapshotROM(double parameter){

    Parameters::AllParameters params = reinitParams(parameter);

    std::shared_ptr<Tests::BurgersRewienskiSnapshot<dim, nstate>> flow_solver_case = std::make_shared<Tests::BurgersRewienskiSnapshot<dim,nstate>>(all_parameters);
    std::unique_ptr<Tests::FlowSolver<dim,nstate>> flow_solver = std::make_unique<Tests::FlowSolver<dim,nstate>>(all_parameters, flow_solver_case);

    // Solve implicit solution
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::pod_galerkin_solver;
    flow_solver->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver->dg, current_pod);
    flow_solver->ode_solver->steady_state();

    // Casting to dg state
    std::shared_ptr< DGBaseState<dim,nstate,double>> dg_state = std::dynamic_pointer_cast< DGBaseState<dim,nstate,double>>(flow_solver->dg);

    // Create functional
    auto functional = BurgersRewienskiFunctional<dim,nstate,double>(flow_solver->dg,dg_state->pde_physics_fad_fad,true,false);

    std::shared_ptr<ROMSolution<dim,nstate>> rom_solution = std::make_shared<ROMSolution<dim, nstate>>(flow_solver->dg, functional, current_pod);

    return rom_solution;
}

template <int dim, int nstate>
Parameters::AllParameters AdaptiveSampling<dim, nstate>::reinitParams(double parameter){
    dealii::ParameterHandler parameter_handler;
    PHiLiP::Parameters::AllParameters::declare_parameters (parameter_handler);
    PHiLiP::Parameters::AllParameters parameters;
    parameters.parse_parameters(parameter_handler);

    // Copy all parameters
    parameters.manufactured_convergence_study_param = all_parameters->manufactured_convergence_study_param;
    parameters.ode_solver_param = all_parameters->ode_solver_param;
    parameters.linear_solver_param = all_parameters->linear_solver_param;
    parameters.euler_param = all_parameters->euler_param;
    parameters.navier_stokes_param = all_parameters->navier_stokes_param;
    parameters.reduced_order_param= all_parameters->reduced_order_param;
    parameters.burgers_param = all_parameters->burgers_param;
    parameters.grid_refinement_study_param = all_parameters->grid_refinement_study_param;
    parameters.artificial_dissipation_param = all_parameters->artificial_dissipation_param;
    parameters.flow_solver_param = all_parameters->flow_solver_param;
    parameters.mesh_adaptation_param = all_parameters->mesh_adaptation_param;
    parameters.artificial_dissipation_param = all_parameters->artificial_dissipation_param;
    parameters.artificial_dissipation_param = all_parameters->artificial_dissipation_param;
    parameters.artificial_dissipation_param = all_parameters->artificial_dissipation_param;
    parameters.dimension = all_parameters->dimension;
    parameters.pde_type = all_parameters->pde_type;
    parameters.use_weak_form = all_parameters->use_weak_form;
    parameters.use_collocated_nodes = all_parameters->use_collocated_nodes;

    using FlowCaseEnum = Parameters::FlowSolverParam::FlowCaseType;
    const FlowCaseEnum flow_type = this->all_parameters->flow_solver_param.flow_case_type;
    if (flow_type == FlowCaseEnum::burgers_rewienski_snapshot){
        parameters.burgers_param.rewienski_b = parameter;
    }
    else if (flow_type == FlowCaseEnum::burgers_viscous_snapshot){
        parameters.burgers_param.diffusion_coefficient = parameter;
    }
    else{
        std::cout << "Invalid flow case. You probably forgot to add it to the list of flow cases in finite_difference_sensitivity.cpp" << std::endl;
        std::abort();
    }

    return parameters;
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
        template class BurgersRewienskiFunctional<PHILIP_DIM, PHILIP_DIM, double>;
#endif


template class AdaptiveSampling<PHILIP_DIM, 1>;
template class AdaptiveSampling<PHILIP_DIM, 2>;
template class AdaptiveSampling<PHILIP_DIM, 3>;
template class AdaptiveSampling<PHILIP_DIM, 4>;
template class AdaptiveSampling<PHILIP_DIM, 5>;

}
}
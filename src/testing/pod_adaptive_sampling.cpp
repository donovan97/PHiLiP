#include "pod_adaptive_sampling.h"

namespace PHiLiP {
namespace Tests {

template<int dim, int nstate>
AdaptiveSampling<dim, nstate>::AdaptiveSampling(const PHiLiP::Parameters::AllParameters *const parameters_input)
    : TestsBase::TestsBase(parameters_input)
    , current_pod(std::make_shared<ProperOrthogonalDecomposition::OnlinePOD<dim>>())
    {
        tolerance = 0.001;
        std::vector<double> parameter_range = {0, 1};
        parameter_space = parameter_range;
    }

template <int dim, int nstate>
int AdaptiveSampling<dim, nstate>::run_test() const
{
    std::cout << "Starting adaptive sampling process" << std::endl;

    //Generate trial locations
    int n = 10;
    generateTrialLocations(n);
    initializeSampling();


    return 0;
}

template <int dim, int nstate>
void AdaptiveSampling<dim, nstate>::generateTrialLocations(int n) const{
    std::cout << "Generating trial locations:" << std::endl;
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
        std::cout << q << std::endl;
        trial_locations.push_back(q);
    }
}

template <int dim, int nstate>
void AdaptiveSampling<dim, nstate>::initializeSampling() const{
    int initialSnapshots = 2;
    for(int idx = 0 ; idx < initialSnapshots ; idx++){
        std::shared_ptr<ProperOrthogonalDecomposition::FOMSolution<dim,nstate>> fom_solution = solveSnapshotFOM(trial_locations[idx]);
        ProperOrthogonalDecomposition::Snapshot<dim,nstate> snapshot = ProperOrthogonalDecomposition::Snapshot<dim,nstate>(trial_locations[idx], fom_solution);
        snapshots.push_back(snapshot);
        current_pod->addSnapshot(fom_solution->state);
    }

    current_pod->computeBasis();

    for(int idx = 0 ; idx < initialSnapshots ; idx++){
        std::shared_ptr<ProperOrthogonalDecomposition::ROMSolution<dim,nstate>> rom_solution = solveSnapshotROM(trial_locations[idx]);
        snapshots[idx].add_ROM(rom_solution);
    }

}

template <int dim, int nstate>
void AdaptiveSampling<dim, nstate>::updatePODBasis() const{

}

template <int dim, int nstate>
std::shared_ptr<ProperOrthogonalDecomposition::FOMSolution<dim,nstate>> AdaptiveSampling<dim, nstate>::solveSnapshotFOM(double parameter) const{

    Parameters::AllParameters params = reinitParams(parameter);

    std::shared_ptr<Tests::BurgersRewienskiSnapshot<dim, nstate>> flow_solver_case = std::make_shared<Tests::BurgersRewienskiSnapshot<dim,nstate>>(all_parameters);
    std::unique_ptr<Tests::FlowSolver<dim,nstate>> flow_solver = std::make_unique<Tests::FlowSolver<dim,nstate>>(all_parameters, flow_solver_case);

    // Solve implicit solution
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::implicit_solver;
    flow_solver->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver->dg);
    flow_solver->ode_solver->allocate_ode_system();
    flow_solver->ode_solver->steady_state();
    flow_solver->flow_solver_case->steady_state_postprocessing(flow_solver->dg);
    // Casting to dg state
    std::shared_ptr< DGBaseState<dim,nstate,double>> dg_state = std::dynamic_pointer_cast< DGBaseState<dim,nstate,double> >(flow_solver->dg);
    // Create functional
    auto functional = BurgersRewienskiFunctional2<dim,nstate,double>(flow_solver->dg,dg_state->pde_physics_fad_fad,true,false);
    //Get sensitivity from FlowSolver
    std::shared_ptr<ProperOrthogonalDecomposition::FOMSolution<dim,nstate>> fom_solution = std::make_shared<ProperOrthogonalDecomposition::FOMSolution<dim, nstate>>(flow_solver->dg, functional, *flow_solver_case->sensitivity_dWdb);
    return fom_solution;
}

template <int dim, int nstate>
std::shared_ptr<ProperOrthogonalDecomposition::ROMSolution<dim,nstate>> AdaptiveSampling<dim, nstate>::solveSnapshotROM(double parameter) const{

    Parameters::AllParameters params = reinitParams(parameter);

    std::shared_ptr<Tests::BurgersRewienskiSnapshot<dim, nstate>> flow_solver_case = std::make_shared<Tests::BurgersRewienskiSnapshot<dim,nstate>>(all_parameters);
    std::unique_ptr<Tests::FlowSolver<dim,nstate>> flow_solver = std::make_unique<Tests::FlowSolver<dim,nstate>>(all_parameters, flow_solver_case);

    // Solve implicit solution
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::pod_galerkin_solver;
    flow_solver->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver->dg, current_pod);
    flow_solver->ode_solver->allocate_ode_system();
    flow_solver->ode_solver->steady_state();

    // Casting to dg state
    std::shared_ptr< DGBaseState<dim,nstate,double>> dg_state = std::dynamic_pointer_cast< DGBaseState<dim,nstate,double>>(flow_solver->dg);

    // Create functional
    auto functional = BurgersRewienskiFunctional2<dim,nstate,double>(flow_solver->dg,dg_state->pde_physics_fad_fad,true,false);

    std::shared_ptr<ProperOrthogonalDecomposition::ROMSolution<dim,nstate>> rom_solution = std::make_shared<ProperOrthogonalDecomposition::ROMSolution<dim, nstate>>(flow_solver->dg, functional, current_pod);

    return rom_solution;
}

template <int dim, int nstate>
Parameters::AllParameters AdaptiveSampling<dim, nstate>::reinitParams(double parameter) const{
    dealii::ParameterHandler parameter_handler;
    PHiLiP::Parameters::AllParameters::declare_parameters (parameter_handler);
    PHiLiP::Parameters::AllParameters parameters;
    parameters.parse_parameters(parameter_handler);

    // Copy all parameters
    parameters.manufactured_convergence_study_param = this->all_parameters->manufactured_convergence_study_param;
    parameters.ode_solver_param = this->all_parameters->ode_solver_param;
    parameters.linear_solver_param = this->all_parameters->linear_solver_param;
    parameters.euler_param = this->all_parameters->euler_param;
    parameters.navier_stokes_param = this->all_parameters->navier_stokes_param;
    parameters.reduced_order_param= this->all_parameters->reduced_order_param;
    parameters.burgers_param = this->all_parameters->burgers_param;
    parameters.grid_refinement_study_param = this->all_parameters->grid_refinement_study_param;
    parameters.artificial_dissipation_param = this->all_parameters->artificial_dissipation_param;
    parameters.flow_solver_param = this->all_parameters->flow_solver_param;
    parameters.mesh_adaptation_param = this->all_parameters->mesh_adaptation_param;
    parameters.artificial_dissipation_param = this->all_parameters->artificial_dissipation_param;
    parameters.artificial_dissipation_param = this->all_parameters->artificial_dissipation_param;
    parameters.artificial_dissipation_param = this->all_parameters->artificial_dissipation_param;
    parameters.dimension = this->all_parameters->dimension;
    parameters.pde_type = this->all_parameters->pde_type;
    parameters.use_weak_form = this->all_parameters->use_weak_form;
    parameters.use_collocated_nodes = this->all_parameters->use_collocated_nodes;

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
real2 BurgersRewienskiFunctional2<dim,nstate,real>::evaluate_volume_integrand(
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
    template class AdaptiveSampling<PHILIP_DIM, PHILIP_DIM>;
    template class BurgersRewienskiFunctional2<PHILIP_DIM, PHILIP_DIM, double>;
#endif


}
}
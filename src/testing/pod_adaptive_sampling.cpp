#include "pod_adaptive_sampling.h"

namespace PHiLiP {
namespace Tests {

template<int dim, int nstate>
AdaptiveSampling<dim, nstate>::AdaptiveSampling(const PHiLiP::Parameters::AllParameters *const parameters_input)
    : TestsBase::TestsBase(parameters_input)
    {
        std::vector<double> parameter_range = {0.01, 0.1};
        parameter_space = parameter_range;

        std::shared_ptr<Tests::BurgersRewienskiSnapshot<dim, nstate>> flow_solver_case = std::make_shared<Tests::BurgersRewienskiSnapshot<dim,nstate>>(all_parameters);
        std::unique_ptr<Tests::FlowSolver<dim,nstate>> flow_solver = std::make_unique<Tests::FlowSolver<dim,nstate>>(all_parameters, flow_solver_case);
        current_pod = std::make_shared<ProperOrthogonalDecomposition::OnlinePOD<dim>>(flow_solver->dg);
    }

template <int dim, int nstate>
int AdaptiveSampling<dim, nstate>::run_test() const
{
    std::cout << "Starting adaptive sampling process" << std::endl;

    int n = 5;
    initializeSampling(n);

    double max_error = 1;
    double tolerance = 1E-03;
    double max_parameter = getMaxErrorROM();

    double parameter = parameter_space[0];
    Parameters::AllParameters params = reinitParams(parameter);
    std::shared_ptr<Tests::BurgersRewienskiSnapshot<dim, nstate>> flow_solver_case = std::make_shared<Tests::BurgersRewienskiSnapshot<dim,nstate>>(&params);
    std::unique_ptr<Tests::FlowSolver<dim,nstate>> flow_solver = std::make_unique<Tests::FlowSolver<dim,nstate>>(&params, flow_solver_case);
    const dealii::LinearAlgebra::distributed::Vector<double> initial_conditions = flow_solver->dg->solution;

    while(abs(max_error) > tolerance){

        std::cout << "Sampling snapshot at " << max_parameter << std::endl;

        std::shared_ptr<ProperOrthogonalDecomposition::FOMSolution<dim,nstate>> fom_solution = solveSnapshotFOM(max_parameter);
        dealii::LinearAlgebra::distributed::Vector<double> state_tmp = fom_solution->state;
        current_pod->addSnapshot(state_tmp -= initial_conditions);
        current_pod->computeBasis();

        //Find adjacent FOMS
        double parameter_before = parameter_space[0];
        double parameter_after = parameter_space[1];
        for(unsigned int i = 0 ; i < sampled_locations.size() ; i++){
            if(sampled_locations[i] >= parameter_before && sampled_locations[i] < max_parameter){
                std::cout << "FOM sampled before is " << sampled_locations[i] << std::endl;
                parameter_before = sampled_locations[i];
            }
            if(sampled_locations[i] <= parameter_after && sampled_locations[i] > max_parameter) {
                std::cout << "FOM sampled after is " << sampled_locations[i] << std::endl;
                parameter_after = sampled_locations[i];
            }
        }

        //Get new ROM locations
        double new_location1 = (parameter_before + max_parameter)/2;
        double new_location2 = (parameter_after + max_parameter)/2;
        //Remove location from ROM sampled locations and add to FOM sampled locations
        sampled_locations.push_back(max_parameter);
        rom_locations.erase(max_parameter);

        std::cout << "Removed ROM at " << max_parameter << " and updating error at other ROM locations." << std::endl;

        //Update previous ROM errors
        for(auto& [key, value] : rom_locations){
            value.compute_initial_rom_to_final_rom_error(current_pod);
            value.compute_total_error();
        }

        std::cout << "Adding ROMs at " << new_location1 << " and " << new_location2 <<std::endl;

        //Sample ROMs at new locations
        std::shared_ptr<ProperOrthogonalDecomposition::ROMSolution<dim,nstate>> rom_solution1 = solveSnapshotROM(new_location1);
        ProperOrthogonalDecomposition::ROMTestLocation<dim,nstate> rom_location1 = ProperOrthogonalDecomposition::ROMTestLocation<dim,nstate>(new_location1, rom_solution1);
        rom_locations.emplace(new_location1, rom_location1);

        std::shared_ptr<ProperOrthogonalDecomposition::ROMSolution<dim,nstate>> rom_solution2 = solveSnapshotROM(new_location2);
        ProperOrthogonalDecomposition::ROMTestLocation<dim,nstate> rom_location2 = ProperOrthogonalDecomposition::ROMTestLocation<dim,nstate>(new_location2, rom_solution2);
        rom_locations.emplace(new_location2, rom_location2);

        max_parameter = getMaxErrorROM();
        max_error = rom_locations.at(max_parameter).total_error;
        std::cout << "Max error is: " << max_error << std::endl;
    }

    std::shared_ptr<dealii::TableHandler> snapshot_table = std::make_shared<dealii::TableHandler>();

    for(auto& location : sampled_locations){
        snapshot_table->add_value("Snapshots", location);
        snapshot_table->set_precision("Snapshots", 16);
    }

    std::ofstream snapshot_table_file("adaptive_sampling_snapshot_table.txt");
    snapshot_table->write_text(snapshot_table_file);

    std::shared_ptr<dealii::TableHandler> rom_table = std::make_shared<dealii::TableHandler>();

    for(auto& [key, value] : rom_locations){
        rom_table->add_value("ROM params", value.parameter);
        rom_table->set_precision("ROM params", 16);

        rom_table->add_value("ROM errors", value.total_error);
        rom_table->set_precision("ROM errors", 16);
    }

    std::ofstream rom_table_file("adaptive_sampling_rom_table.txt");
    rom_table->write_text(rom_table_file);

    return 0;
}


template <int dim, int nstate>
double AdaptiveSampling<dim, nstate>::getMaxErrorROM() const{
    std::cout << "Updating error curve fit..." << std::endl;

    //Set errors at ROM locations
    double max_parameter = 0;
    double max_error = 0;
    for(auto& [key, value] : rom_locations){
        if(abs(max_error) < abs(value.total_error)){
            max_parameter = value.parameter;
            max_error = value.total_error;
        }
        std::cout << "Parameter: " << value.parameter << " Error: " << value.total_error << std::endl;
    }
    return max_parameter;
}

template <int dim, int nstate>
void AdaptiveSampling<dim, nstate>::initializeSampling(int n) const{

    double dx = (parameter_space[1] - parameter_space[0])/(n - 1);

    //Get initial conditions
    double parameter = parameter_space[0];
    Parameters::AllParameters params = reinitParams(parameter);
    std::shared_ptr<Tests::BurgersRewienskiSnapshot<dim, nstate>> flow_solver_case = std::make_shared<Tests::BurgersRewienskiSnapshot<dim,nstate>>(&params);
    std::unique_ptr<Tests::FlowSolver<dim,nstate>> flow_solver = std::make_unique<Tests::FlowSolver<dim,nstate>>(&params, flow_solver_case);
    const dealii::LinearAlgebra::distributed::Vector<double> initial_conditions = flow_solver->dg->solution;
    //current_pod->addSnapshot(initial_conditions);

    for(int i = 0 ; i < n ; i++){
        parameter = i*dx + parameter_space[0];
        std::cout << "Sampling initial snapshot at " << parameter << std::endl;
        sampled_locations.push_back(parameter);
        std::shared_ptr<ProperOrthogonalDecomposition::FOMSolution<dim,nstate>> fom_solution = solveSnapshotFOM(parameter);
        dealii::LinearAlgebra::distributed::Vector<double> state_tmp = fom_solution->state;
        current_pod->addSnapshot(state_tmp -= initial_conditions);
        //rom_locations.push_back(*fom_solution);
    }

    current_pod->computeBasis();

    for(int i = 0 ; i < n-1 ; i++) {
        parameter = i * dx + dx / 2 + parameter_space[0];
        std::cout << "Computing ROM at " << parameter << std::endl;
        std::shared_ptr<ProperOrthogonalDecomposition::ROMSolution<dim, nstate>> rom_solution = solveSnapshotROM(parameter);
        ProperOrthogonalDecomposition::ROMTestLocation < dim,nstate > rom_location = ProperOrthogonalDecomposition::ROMTestLocation < dim, nstate>(parameter, rom_solution);
        rom_locations.emplace(parameter, rom_location);
    }

}

template <int dim, int nstate>
dealii::Vector<double> AdaptiveSampling<dim, nstate>::polyFit(dealii::Vector<double> x, dealii::Vector<double> y, unsigned int n) const{

    dealii::FullMatrix<double> vandermonde(x.size(), n+1);

    for(unsigned int i = 0 ; i < x.size() ; i++){
        for(unsigned int j = 0 ; j < n+1 ; j++){
            vandermonde.set(i, j, std::pow(x(i), j));
        }
    }

    dealii::Householder<double> householder(vandermonde);

    dealii::Vector<double> coefficients(n+1);
    householder.least_squares(coefficients, y);

    return coefficients;
}

template <int dim, int nstate>
dealii::Vector<double> AdaptiveSampling<dim, nstate>::polyVal(dealii::Vector<double> polynomial, dealii::Vector<double> x) const{

    dealii::FullMatrix<double> vandermonde(x.size(), polynomial.size());

    for(unsigned int i = 0 ; i < x.size() ; i++){
        for(unsigned int j = 0 ; j < polynomial.size() ; j++){
            vandermonde.set(i, j, std::pow(x(i), j));
        }
    }
    dealii::Vector<double> y(x.size());
    vandermonde.vmult(y, polynomial);

    return y;

}

template <int dim, int nstate>
std::shared_ptr<ProperOrthogonalDecomposition::FOMSolution<dim,nstate>> AdaptiveSampling<dim, nstate>::solveSnapshotFOM(double parameter) const{
    std::cout << "Solving FOM at " << parameter << std::endl;
    Parameters::AllParameters params = reinitParams(parameter);

    std::shared_ptr<Tests::BurgersRewienskiSnapshot<dim, nstate>> flow_solver_case = std::make_shared<Tests::BurgersRewienskiSnapshot<dim,nstate>>(&params);
    std::unique_ptr<Tests::FlowSolver<dim,nstate>> flow_solver = std::make_unique<Tests::FlowSolver<dim,nstate>>(&params, flow_solver_case);

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
    std::shared_ptr<ProperOrthogonalDecomposition::FOMSolution<dim,nstate>> fom_solution = std::make_shared<ProperOrthogonalDecomposition::FOMSolution<dim, nstate>>(flow_solver->dg, functional, flow_solver_case->sensitivity_dWdb_l2norm);

    std::cout << "Done solving FOM." << std::endl;
    return fom_solution;
}

template <int dim, int nstate>
std::shared_ptr<ProperOrthogonalDecomposition::ROMSolution<dim,nstate>> AdaptiveSampling<dim, nstate>::solveSnapshotROM(double parameter) const{
    std::cout << "Solving ROM at " << parameter << std::endl;
    Parameters::AllParameters params = reinitParams(parameter);

    std::shared_ptr<Tests::BurgersRewienskiSnapshot<dim, nstate>> flow_solver_case = std::make_shared<Tests::BurgersRewienskiSnapshot<dim,nstate>>(&params);
    std::unique_ptr<Tests::FlowSolver<dim,nstate>> flow_solver = std::make_unique<Tests::FlowSolver<dim,nstate>>(&params, flow_solver_case);

    // Solve implicit solution
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::pod_petrov_galerkin_solver;
    flow_solver->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver->dg, current_pod);
    flow_solver->ode_solver->allocate_ode_system();
    flow_solver->ode_solver->steady_state();

    // Casting to dg state
    std::shared_ptr< DGBaseState<dim,nstate,double>> dg_state = std::dynamic_pointer_cast< DGBaseState<dim,nstate,double>>(flow_solver->dg);

    // Create functional
    auto functional = BurgersRewienskiFunctional2<dim,nstate,double>(flow_solver->dg,dg_state->pde_physics_fad_fad,true,false);

    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> system_matrix_transpose = std::make_shared<dealii::TrilinosWrappers::SparseMatrix>();
    system_matrix_transpose->copy_from(flow_solver->dg->system_matrix_transpose);

    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> pod_basis = std::make_shared<dealii::TrilinosWrappers::SparseMatrix>();
    pod_basis->copy_from(*current_pod->getPODBasis());

    std::shared_ptr<ProperOrthogonalDecomposition::ROMSolution<dim,nstate>> rom_solution = std::make_shared<ProperOrthogonalDecomposition::ROMSolution<dim, nstate>>(flow_solver->dg, system_matrix_transpose,functional, pod_basis);
    std::cout << "Done solving ROM." << std::endl;
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
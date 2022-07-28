#include "pod_adaptive_sampling.h"
#include <iostream>
#include <filesystem>
#include "functional/functional.h"
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector_operation.h>
#include "reduced_order/reduced_order_solution.h"
#include "flow_solver/flow_solver.h"
#include "flow_solver/flow_solver_factory.h"
#include <cmath>
#include "reduced_order/rbf_interpolation.h"
#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Stream.hpp"
#include "ROL_Bounds.hpp"
#include "reduced_order/halton.h"
#include "reduced_order/min_max_scaler.h"
#include "tests.h"

namespace PHiLiP {
namespace Tests {

template<int dim, int nstate>
AdaptiveSampling<dim, nstate>::AdaptiveSampling(const PHiLiP::Parameters::AllParameters *const parameters_input,
                                                const dealii::ParameterHandler &parameter_handler_input)
        : TestsBase::TestsBase(parameters_input)
        , parameter_handler(parameter_handler_input)
{
    configureParameterSpace();
    std::unique_ptr<FlowSolver::FlowSolver<dim,nstate>> flow_solver = FlowSolver::FlowSolverFactory<dim,nstate>::select_flow_case(all_parameters, parameter_handler);
    const bool compute_dRdW = true;
    flow_solver->dg->assemble_residual(compute_dRdW);
    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> system_matrix = std::make_shared<dealii::TrilinosWrappers::SparseMatrix>();
    system_matrix->copy_from(flow_solver->dg->system_matrix);
    current_pod = std::make_shared<ProperOrthogonalDecomposition::OnlinePOD<dim>>(system_matrix);
    nearest_neighbors = std::make_shared<ProperOrthogonalDecomposition::NearestNeighbors>();
}

template <int dim, int nstate>
int AdaptiveSampling<dim, nstate>::run_test() const
{
    this->pcout << "Starting adaptive sampling process" << std::endl;

    placeInitialSnapshots();
    current_pod->computeBasis();

    MatrixXd rom_points = nearest_neighbors->kPairwiseNearestNeighborsMidpoint();
    pcout << rom_points << std::endl;

    placeTriangulationROMs(rom_points);

    RowVectorXd max_error_params = getMaxErrorROM();
    int iteration = 0;

    while(max_error > all_parameters->reduced_order_param.adaptation_tolerance){

        outputErrors(iteration);

        this->pcout << "Sampling snapshot at " << max_error_params << std::endl;
        dealii::LinearAlgebra::distributed::Vector<double> fom_solution = solveSnapshotFOM(max_error_params);
        snapshot_parameters.conservativeResize(snapshot_parameters.rows()+1, snapshot_parameters.cols());
        snapshot_parameters.row(snapshot_parameters.rows()-1) = max_error_params;
        nearest_neighbors->updateSnapshots(snapshot_parameters, fom_solution);
        current_pod->addSnapshot(fom_solution);
        current_pod->computeBasis();

        //Update previous ROM errors with updated current_pod
        for(auto it = rom_locations.begin(); it != rom_locations.end(); ++it){
            it->second->compute_initial_rom_to_final_rom_error(current_pod);
            it->second->compute_total_error();
        }

        //updateNearestExistingROMs(max_error_params);

        rom_points = nearest_neighbors->kNearestNeighborsMidpoint(max_error_params);
        pcout << rom_points << std::endl;

        bool error_greater_than_tolerance = placeTriangulationROMs(rom_points);


        if(!error_greater_than_tolerance){
            updateNearestExistingROMs(max_error_params);
        }


        //updateNearestExistingROMs(max_error_params);

        //Update max error
        max_error_params = getMaxErrorROM();

        this->pcout << "Max error is: " << max_error << std::endl;
        iteration++;
    }

    outputErrors(iteration);

    return 0;
}

template <int dim, int nstate>
void AdaptiveSampling<dim, nstate>::outputErrors(int iteration) const{
    std::unique_ptr<dealii::TableHandler> snapshot_table = std::make_unique<dealii::TableHandler>();

    std::ofstream solution_out_file("solution_snapshots_iteration_" +  std::to_string(iteration) + ".txt");
    unsigned int precision = 16;
    current_pod->dealiiSnapshotMatrix.print_formatted(solution_out_file, precision);
    solution_out_file.close();

    for(auto parameters : snapshot_parameters.rowwise()){
        for(int i = 0 ; i < snapshot_parameters.cols() ; i++){
            snapshot_table->add_value(all_parameters->reduced_order_param.parameter_names[i], parameters(i));
            snapshot_table->set_precision(all_parameters->reduced_order_param.parameter_names[i], 16);
        }
    }

    std::ofstream snapshot_table_file("snapshot_table_iteration_" + std::to_string(iteration) + ".txt");
    snapshot_table->write_text(snapshot_table_file, dealii::TableHandler::TextOutputFormat::org_mode_table);
    snapshot_table_file.close();

    std::unique_ptr<dealii::TableHandler> rom_table = std::make_unique<dealii::TableHandler>();

    for(auto it = rom_locations.begin(); it != rom_locations.end(); ++it){
        for(int i = 0 ; i < snapshot_parameters.cols() ; i++){
            rom_table->add_value(all_parameters->reduced_order_param.parameter_names[i], it->first(i));
            rom_table->set_precision(all_parameters->reduced_order_param.parameter_names[i], 16);
        }
        rom_table->add_value("ROM_errors", it->second->total_error);
        rom_table->set_precision("ROM_errors", 16);
    }

    std::ofstream rom_table_file("rom_table_iteration_" + std::to_string(iteration) + ".txt");
    rom_table->write_text(rom_table_file, dealii::TableHandler::TextOutputFormat::org_mode_table);
    rom_table_file.close();
}

template <int dim, int nstate>
RowVectorXd AdaptiveSampling<dim, nstate>::getMaxErrorROM() const{
    this->pcout << "Updating RBF interpolation..." << std::endl;

    int n_rows = snapshot_parameters.rows() + rom_locations.size();
    MatrixXd parameters(n_rows, snapshot_parameters.cols());
    VectorXd errors(n_rows);

    int i;
    for(i = 0 ; i < snapshot_parameters.rows() ; i++){
        errors(i) = 0;
        parameters.row(i) = snapshot_parameters.row(i);
    }
    this->pcout << i << std::endl;
    for(auto it = rom_locations.begin(); it != rom_locations.end(); ++it){
        parameters.row(i) = it->first.array();
        errors(i) = it->second->total_error;
        i++;
    }

    //Must scale both axes between [0,1] for the 2d rbf interpolation to work optimally
    ProperOrthogonalDecomposition::MinMaxScaler scaler;
    MatrixXd parameters_scaled = scaler.fit_transform(parameters);

    //Construct radial basis function
    std::string kernel = "thin_plate_spline";
    ProperOrthogonalDecomposition::RBFInterpolation rbf = ProperOrthogonalDecomposition::RBFInterpolation(parameters_scaled, errors, kernel);

    // Set parameters.
    ROL::ParameterList parlist;
    //parlist.sublist("General").set("Recompute Objective Function",false);
    //parlist.sublist("Step").sublist("Line Search").set("Initial Step Size", 1.0);
    //parlist.sublist("Step").sublist("Line Search").set("User Defined Initial Step Size",true);
    parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").set("Type","Backtracking");
    parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type","Newton-Krylov");
    parlist.sublist("Step").sublist("Line Search").sublist("Curvature Condition").set("Type","Null Curvature Condition");
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-10);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist.sublist("Status Test").set("Iteration Limit",100);

    //Find max error and parameters by minimizing function starting at each ROM location
    RowVectorXd max_error_params(parameters.cols());
    max_error = 0;

    for(auto it = rom_locations.begin(); it != rom_locations.end(); ++it){

        Eigen::RowVectorXd rom_unscaled = it->first;
        Eigen::RowVectorXd rom_scaled = scaler.transform(rom_unscaled);

        //start bounds
        int dimension = parameters.cols();
        ROL::Ptr<std::vector<double>> l_ptr = ROL::makePtr<std::vector<double>>(dimension,0.0);
        ROL::Ptr<std::vector<double>> u_ptr = ROL::makePtr<std::vector<double>>(dimension,1.0);
        ROL::Ptr<ROL::Vector<double>> lo = ROL::makePtr<ROL::StdVector<double>>(l_ptr);
        ROL::Ptr<ROL::Vector<double>> up = ROL::makePtr<ROL::StdVector<double>>(u_ptr);
        ROL::Bounds<double> bcon(lo,up);
        //end bounds

        ROL::Ptr<ROL::Step<double>> step = ROL::makePtr<ROL::LineSearchStep<double>>(parlist);
        ROL::Ptr<ROL::StatusTest<double>> status = ROL::makePtr<ROL::StatusTest<double>>(parlist);
        ROL::Algorithm<double> algo(step,status,false);

        // Iteration Vector
        ROL::Ptr<std::vector<double>> x_ptr = ROL::makePtr<std::vector<double>>(dimension, 0.0);

        // Set Initial Guess
        for(int j = 0 ; j < dimension ; j++){
            (*x_ptr)[j] = rom_scaled(j);
        }

        this->pcout << "Unscaled parameter: " << rom_unscaled << std::endl;
        this->pcout << "Scaled parameter: ";
        for(int j = 0 ; j < dimension ; j++){
            this->pcout << (*x_ptr)[j] << " ";
        }
        this->pcout << std::endl;

        ROL::StdVector<double> x(x_ptr);

        // Run Algorithm
        //ROL::Ptr<std::ostream> outStream;
        //outStream = ROL::makePtrFromRef(std::cout);
        //algo.run(x, rbf, bcon,true, *outStream);
        algo.run(x, rbf, bcon,false);

        ROL::Ptr<std::vector<double>> x_min = x.getVector();

        for(int j = 0 ; j < dimension ; j++){
            rom_scaled(j) = (*x_min)[j];
        }

        RowVectorXd rom_unscaled_optim = scaler.inverse_transform(rom_scaled);

        this->pcout << "Parameters of optimization convergence: " << rom_unscaled_optim << std::endl;

        double error = std::abs(rbf.evaluate(rom_scaled));
        this->pcout << "RBF error at optimization convergence: " << error << std::endl;
        if(error > max_error){
            this->pcout << "RBF error is greater than current max error. Updating max error." << std::endl;
            max_error = error;
            max_error_params = rom_unscaled_optim;
            this->pcout << "RBF Max error: " << max_error << std::endl;
        }
    }

    //Check if max_error_params is a ROM point
    for(auto it = rom_locations.begin(); it != rom_locations.end(); ++it){
        if(max_error_params.isApprox(it->first)){
            this->pcout << "Max error location is approximately the same as a ROM location. Removing ROM location." << std::endl;
            rom_locations.erase(it);
            break;
        }
    }

    return max_error_params;
}

template <int dim, int nstate>
void AdaptiveSampling<dim, nstate>::placeInitialSnapshots() const{
    for(auto snap_param : snapshot_parameters.rowwise()){
        this->pcout << "Sampling initial snapshot at " << snap_param << std::endl;
        dealii::LinearAlgebra::distributed::Vector<double> fom_solution = solveSnapshotFOM(snap_param);
        nearest_neighbors->updateSnapshots(snapshot_parameters, fom_solution);
        current_pod->addSnapshot(fom_solution);
    }
}

template <int dim, int nstate>
bool AdaptiveSampling<dim, nstate>::placeTriangulationROMs(const MatrixXd& rom_points) const{
    bool error_greater_than_tolerance = false;
    for(auto midpoint : rom_points.rowwise()){

        //Check if ROM point already exists as another ROM point
        auto element = std::find_if(rom_locations.begin(), rom_locations.end(), [&midpoint](std::pair<RowVectorXd, std::shared_ptr<ProperOrthogonalDecomposition::ROMTestLocation<dim,nstate>>>& location){ return location.first.isApprox(midpoint);} );

        //Check if ROM point already exists as a snapshot
        bool snapshot_exists = false;
        for(auto snap_param : snapshot_parameters.rowwise()){
            if(snap_param.isApprox(midpoint)){
                snapshot_exists = true;
            }
        }

        if(element == rom_locations.end() && snapshot_exists == false){
            ProperOrthogonalDecomposition::ROMSolution<dim, nstate> rom_solution = solveSnapshotROM(midpoint);
            std::shared_ptr<ProperOrthogonalDecomposition::ROMTestLocation < dim,nstate >> rom_location = std::make_shared<ProperOrthogonalDecomposition::ROMTestLocation < dim, nstate>>(midpoint, std::move(rom_solution));
            rom_locations.emplace_back(std::make_pair(midpoint, rom_location));
            if(abs(rom_location->total_error) > all_parameters->reduced_order_param.adaptation_tolerance){
                error_greater_than_tolerance = true;
            }
        }
        else{
            this->pcout << "ROM already computed." << std::endl;
        }
    }
    return error_greater_than_tolerance;
}

template <int dim, int nstate>
void AdaptiveSampling<dim, nstate>::updateNearestExistingROMs(const RowVectorXd& parameter) const{

    /*

    for(auto it = rom_locations.begin(); it != rom_locations.end(); ++it){
        if(std::abs(it->second->total_error) > all_parameters->reduced_order_param.adaptation_tolerance && std::abs(it->second->total_error) > std::abs(it->second->minimum_error) && it->second->iteration_count > 4){
            this->pcout << "Recomputing point: " << it->first << std::endl;
            ProperOrthogonalDecomposition::ROMSolution<dim, nstate> rom_solution = solveSnapshotROM(it->first);
            std::shared_ptr<ProperOrthogonalDecomposition::ROMTestLocation < dim,nstate >> rom_location = std::make_shared<ProperOrthogonalDecomposition::ROMTestLocation < dim, nstate>>(it->first, std::move(rom_solution));
            *it = std::make_pair(it->first, rom_location);
        }
    }
     */

    /*
    this->pcout << "Computing mean error: " << std::endl;
    //Compute mean error
    double error_sum = 0;
    for(auto it = rom_locations.begin(); it != rom_locations.end(); ++it){
        error_sum = error_sum + std::abs(it->second->total_error);
    }
    double error_mean = error_sum / rom_locations.size();
    this->pcout << "Mean error: " << error_mean << std::endl;

    for(auto it = rom_locations.begin(); it != rom_locations.end(); ++it){
        if(std::abs(it->second->total_error) > error_mean && it->second->initial_rom_to_final_rom_error < -1*all_parameters->reduced_order_param.adaptation_tolerance){
            this->pcout << "Recomputing point: " << it->first << std::endl;
            ProperOrthogonalDecomposition::ROMSolution<dim, nstate> rom_solution = solveSnapshotROM(it->first);
            std::shared_ptr<ProperOrthogonalDecomposition::ROMTestLocation < dim,nstate >> rom_location = std::make_shared<ProperOrthogonalDecomposition::ROMTestLocation < dim, nstate>>(it->first, std::move(rom_solution));
            *it = std::make_pair(it->first, rom_location);
        }
    }

    */



    //Assemble ROM points in a matrix
    MatrixXd rom_points(0,0);
    for(auto it = rom_locations.begin(); it != rom_locations.end(); ++it){
        rom_points.conservativeResize(rom_points.rows()+1, 2);
        rom_points.row(rom_points.rows()-1) = it->first;
    }

    //Get distances between newly added snapshots and ROM points
    ProperOrthogonalDecomposition::MinMaxScaler scaler;
    MatrixXd scaled_rom_params = scaler.fit_transform(rom_points);
    RowVectorXd scaled_point = scaler.transform(parameter);

    VectorXd distances = (scaled_rom_params.rowwise() - scaled_point).rowwise().squaredNorm();

    std::vector<int> index(distances.size());
    std::iota(index.begin(), index.end(), 0);

    std::sort(index.begin(), index.end(),
              [&](const int& a, const int& b) {
                  return distances[a] < distances[b];
              });

    //For 2(n+1) nearest points, if error>tolerance, recompute adjoint 1
    pcout << "Searching ROM points near: " << parameter << std::endl;
    for(int i = 0 ; i < 6 ; i++){
        pcout << "ROM point: " << rom_locations[index[i]].first << " Error: " << rom_locations[index[i]].second->total_error << std::endl;
        if(std::abs(rom_locations[index[i]].second->total_error) > all_parameters->reduced_order_param.adaptation_tolerance){
            pcout << "Total error greater than tolerance. Recomputing ROM solution" << std::endl;
            ProperOrthogonalDecomposition::ROMSolution<dim, nstate> rom_solution = solveSnapshotROM(rom_locations[index[i]].first);
            std::shared_ptr<ProperOrthogonalDecomposition::ROMTestLocation < dim,nstate >> rom_location = std::make_shared<ProperOrthogonalDecomposition::ROMTestLocation < dim, nstate>>(rom_locations[index[i]].first, rom_solution);
            rom_locations[index[i]] = std::make_pair(rom_locations[index[i]].first, rom_location);
        }
    }

}

template <int dim, int nstate>
dealii::LinearAlgebra::distributed::Vector<double> AdaptiveSampling<dim, nstate>::solveSnapshotFOM(const RowVectorXd& parameter) const{
    this->pcout << "Solving FOM at " << parameter << std::endl;
    Parameters::AllParameters params = reinitParams(parameter);

    std::unique_ptr<FlowSolver::FlowSolver<dim,nstate>> flow_solver = FlowSolver::FlowSolverFactory<dim,nstate>::select_flow_case(&params, parameter_handler);

    // Solve implicit solution
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::implicit_solver;
    flow_solver->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver->dg);
    flow_solver->ode_solver->allocate_ode_system();
    flow_solver->run();

    this->pcout << "Done solving FOM." << std::endl;
    return flow_solver->dg->solution;
}

template <int dim, int nstate>
ProperOrthogonalDecomposition::ROMSolution<dim,nstate> AdaptiveSampling<dim, nstate>::solveSnapshotROM(const RowVectorXd& parameter) const{
    this->pcout << "Solving ROM at " << parameter << std::endl;
    Parameters::AllParameters params = reinitParams(parameter);

    std::unique_ptr<FlowSolver::FlowSolver<dim,nstate>> flow_solver = FlowSolver::FlowSolverFactory<dim,nstate>::select_flow_case(&params, parameter_handler);

    // Solve implicit solution
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::pod_petrov_galerkin_solver;
    flow_solver->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, flow_solver->dg, current_pod);
    //flow_solver->dg->solution = nearest_neighbors->nearestNeighborMidpointSolution(parameter);
    flow_solver->ode_solver->allocate_ode_system();
    flow_solver->ode_solver->steady_state();

    // Create functional
    std::shared_ptr<Functional<dim,nstate,double>> functional = FunctionalFactory<dim,nstate,double>::create_Functional(params.functional_param, flow_solver->dg);
    functional->evaluate_functional( true, false, false);

    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> system_matrix_transpose = std::make_shared<dealii::TrilinosWrappers::SparseMatrix>();
    system_matrix_transpose->copy_from(flow_solver->dg->system_matrix_transpose);
    dealii::LinearAlgebra::distributed::Vector<double> right_hand_side(flow_solver->dg->right_hand_side);

    dealii::LinearAlgebra::distributed::Vector<double> gradient(functional->dIdw);
    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> pod_basis = std::make_shared<dealii::TrilinosWrappers::SparseMatrix>();
    pod_basis->copy_from(*current_pod->getPODBasis());

    ProperOrthogonalDecomposition::ROMSolution<dim,nstate> rom_solution = ProperOrthogonalDecomposition::ROMSolution<dim, nstate>(system_matrix_transpose, right_hand_side, pod_basis, gradient);
    this->pcout << "Done solving ROM." << std::endl;
    return rom_solution;
}

template <int dim, int nstate>
Parameters::AllParameters AdaptiveSampling<dim, nstate>::reinitParams(const RowVectorXd& parameter) const{
    // Copy all parameters
    PHiLiP::Parameters::AllParameters parameters = *(this->all_parameters);

    using FlowCaseEnum = Parameters::FlowSolverParam::FlowCaseType;
    const FlowCaseEnum flow_type = this->all_parameters->flow_solver_param.flow_case_type;

    if (flow_type == FlowCaseEnum::burgers_rewienski_snapshot){
        if(all_parameters->reduced_order_param.parameter_names.size() == 1){
            if(all_parameters->reduced_order_param.parameter_names[0] == "rewienski_a"){
                parameters.burgers_param.rewienski_a = parameter(0);
            }
            else if(all_parameters->reduced_order_param.parameter_names[0] == "rewienski_b"){
                parameters.burgers_param.rewienski_b = parameter(0);
            }
        }
        else{
            parameters.burgers_param.rewienski_a = parameter(0);
            parameters.burgers_param.rewienski_b = parameter(1);
        }
    }
    else if (flow_type == FlowCaseEnum::naca0012){
        if(all_parameters->reduced_order_param.parameter_names.size() == 1){
            if(all_parameters->reduced_order_param.parameter_names[0] == "mach"){
                parameters.euler_param.mach_inf = parameter(0);
            }
            else if(all_parameters->reduced_order_param.parameter_names[0] == "alpha"){
                parameters.euler_param.angle_of_attack = parameter(0); //radians!
            }
        }
        else{
            parameters.euler_param.mach_inf = parameter(0);
            parameters.euler_param.angle_of_attack = parameter(1); //radians!
        }
    }
    else if (flow_type == FlowCaseEnum::gaussian_bump){
        if(all_parameters->reduced_order_param.parameter_names.size() == 1){
            if(all_parameters->reduced_order_param.parameter_names[0] == "mach"){
                parameters.euler_param.mach_inf = parameter(0);
            }
        }
    }
    else{
        this->pcout << "Invalid flow case. You probably forgot to specify a flow case in the prm file." << std::endl;
        std::abort();
    }
    return parameters;
}

template <int dim, int nstate>
void AdaptiveSampling<dim, nstate>::configureParameterSpace() const
{
    const double pi = atan(1.0) * 4.0;

    int n_halton = all_parameters->reduced_order_param.num_halton;

    if(all_parameters->reduced_order_param.parameter_names.size() == 1){
        RowVectorXd parameter1_range;
        parameter1_range.resize(2);
        parameter1_range << all_parameters->reduced_order_param.parameter_min_values[0], all_parameters->reduced_order_param.parameter_max_values[0];
        if(all_parameters->reduced_order_param.parameter_names[0] == "alpha"){
            parameter1_range *= pi/180; //convert to radians
        }

        //Place parameters at 2 ends and center
        snapshot_parameters.resize(3,1);
        snapshot_parameters  << parameter1_range[0],
                                parameter1_range[1],
                                (parameter1_range[0]+parameter1_range[1])/2;

        snapshot_parameters.conservativeResize(snapshot_parameters.rows() + n_halton, snapshot_parameters.cols());

        double *seq = nullptr;
        for (int i = 0; i < n_halton; i++)
        {
            seq = ProperOrthogonalDecomposition::halton(i+2, 1); //ignore the first two Halton point as they are the left end and center
                snapshot_parameters(i+3) = seq[0]*(parameter1_range[1] - parameter1_range[0])+parameter1_range[0];
        }
        delete [] seq;

        this->pcout << snapshot_parameters << std::endl;
    }
    else if(all_parameters->reduced_order_param.parameter_names.size() == 2){
        RowVectorXd parameter1_range;
        parameter1_range.resize(2);
        parameter1_range << all_parameters->reduced_order_param.parameter_min_values[0], all_parameters->reduced_order_param.parameter_max_values[0];

        RowVectorXd parameter2_range;
        parameter2_range.resize(2);
        parameter2_range << all_parameters->reduced_order_param.parameter_min_values[1], all_parameters->reduced_order_param.parameter_max_values[1];
        if(all_parameters->reduced_order_param.parameter_names[1] == "alpha"){
            parameter2_range *= pi/180; //convert to radians
        }

        //Place parameters at 4 corners and center
        snapshot_parameters.resize(5,2);
        snapshot_parameters  << parameter1_range[0], parameter2_range[0],
                                parameter1_range[0], parameter2_range[1],
                                parameter1_range[1], parameter2_range[1],
                                parameter1_range[1], parameter2_range[0],
                                0.5*(parameter1_range[1] - parameter1_range[0])+parameter1_range[0], 0.5*(parameter2_range[1] - parameter2_range[0])+parameter2_range[0];

        /*
        //Place 9 parameters in a grid
        snapshot_parameters.resize(9,2);
        snapshot_parameters  << parameter1_range[0], parameter2_range[0],
                                parameter1_range[0], parameter2_range[1],
                                parameter1_range[1], parameter2_range[1],
                                parameter1_range[1], parameter2_range[0],
                                0.5*(parameter1_range[1] - parameter1_range[0])+parameter1_range[0], parameter2_range[1],
                                0.5*(parameter1_range[1] - parameter1_range[0])+parameter1_range[0], parameter2_range[0],
                                parameter1_range[0], 0.5*(parameter2_range[1] - parameter2_range[0])+parameter2_range[0],
                                parameter1_range[1], 0.5*(parameter2_range[1] - parameter2_range[0])+parameter2_range[0],
                                0.5*(parameter1_range[1] - parameter1_range[0])+parameter1_range[0], 0.5*(parameter2_range[1] - parameter2_range[0])+parameter2_range[0];
        snapshot_parameters.conservativeResize(snapshot_parameters.rows() + n_halton, snapshot_parameters.cols());
        */
        double *seq = nullptr;
        for (int i = 0; i < n_halton; i++)
        {
            seq = ProperOrthogonalDecomposition::halton(i+1, 2); //ignore the first Halton point as it is one of the corners
            for (int j = 0; j < 2; j++)
            {
                if(j == 0){
                    snapshot_parameters(i+5, j) = seq[j]*(parameter1_range[1] - parameter1_range[0])+parameter1_range[0];
                }
                else if(j == 1){
                    snapshot_parameters(i+5, j) = seq[j]*(parameter2_range[1] - parameter2_range[0])+parameter2_range[0];
                }
            }
        }
        delete [] seq;

        this->pcout << snapshot_parameters << std::endl;
    }
}

#if PHILIP_DIM==1
        template class AdaptiveSampling<PHILIP_DIM, PHILIP_DIM>;
#endif

#if PHILIP_DIM!=1
        template class AdaptiveSampling<PHILIP_DIM, PHILIP_DIM+2>;
#endif

}
}
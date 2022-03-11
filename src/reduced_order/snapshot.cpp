#include "snapshot.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template <int dim, int nstate>
Snapshot<dim, nstate>::Snapshot(double parameter, std::shared_ptr<ROMSolution<dim, nstate>> rom_solution)
        : parameter(parameter)
        , rom_solution(rom_solution)
{
}

template <int dim, int nstate>
Snapshot<dim, nstate>::Snapshot(double parameter, std::shared_ptr<FOMSolution<dim, nstate>> fom_solution)
        : parameter(parameter)
        , fom_solution(fom_solution)
{
}

template <int dim, int nstate>
void Snapshot<dim, nstate>::add_FOM(std::shared_ptr<FOMSolution<dim, nstate>> fom_solution_input){
    fom_solution = fom_solution_input;
    compute_FOM_to_initial_ROM_error();
    initial_rom_to_final_rom_error = 0;
    compute_total_error();
}

template <int dim, int nstate>
void Snapshot<dim, nstate>::add_ROM(std::shared_ptr<ROMSolution<dim, nstate>> rom_solution_input){
    rom_solution = rom_solution_input;
    compute_FOM_to_initial_ROM_error();
    initial_rom_to_final_rom_error = 0;
    compute_total_error();
}

template <int dim, int nstate>
void Snapshot<dim, nstate>::compute_FOM_to_initial_ROM_error(){
    fom_to_initial_rom_error = rom_solution->functional_value - fom_solution->functional_value;
}

template <int dim, int nstate>
void Snapshot<dim, nstate>::compute_initial_rom_to_final_rom_error(std::shared_ptr<ProperOrthogonalDecomposition::POD<dim>> pod_updated){
    using DealiiVector = dealii::LinearAlgebra::distributed::Vector<double>;
    //Initialize
    DealiiVector fineGradient(pod_updated->getPODBasis()->n());
    DealiiVector fineAdjoint(pod_updated->getPODBasis()->n());
    DealiiVector fineResidual(pod_updated->getPODBasis()->n());
    DealiiVector dualWeightedResidual(pod_updated->getPODBasis()->n());

    pod_updated->getPODBasis()->Tvmult(fineGradient, rom_solution->gradient);

    dealii::TrilinosWrappers::SparseMatrix tmp;
    dealii::TrilinosWrappers::SparseMatrix fineJacobianTranspose;
    pod_updated->getPODBasis()->Tmmult(tmp, *rom_solution->system_matrix_transpose); //tmp = pod_basis^T * dg->system_matrix_transpose
    tmp.mmult(fineJacobianTranspose, *pod_updated->getPODBasis()); // reducedJacobianTranspose= tmp*pod_basis

    dealii::ParameterHandler parameter_handler;
    Parameters::LinearSolverParam linear_solver_param;
    Parameters::LinearSolverParam::declare_parameters (parameter_handler);
    linear_solver_param.parse_parameters (parameter_handler);
    linear_solver_param.linear_solver_type = Parameters::LinearSolverParam::direct;
    solve_linear(fineJacobianTranspose, fineGradient*=-1.0, fineAdjoint, linear_solver_param);

    //Compute fine residual
    pod_updated->getPODBasis()->Tvmult(fineResidual, rom_solution->right_hand_side);

    //Compute dual weighted residual
    initial_rom_to_final_rom_error = 0;
    std::cout << std::setw(10) << std::left << "Index" << std::setw(20) << std::left << "Reduced Adjoint" << std::setw(20) << std::left << "Reduced Residual" << std::setw(20) << std::left << "Dual Weighted Residual" << std::endl;
    for(unsigned int i = 0; i < fineAdjoint.size(); i++){
        dualWeightedResidual[i] = -(fineAdjoint[i] * fineResidual[i]);
        initial_rom_to_final_rom_error = initial_rom_to_final_rom_error + dualWeightedResidual[i];
        std::cout << std::setw(10) << std::left << i << std::setw(20) << std::left << fineAdjoint[i] << std::setw(20) << std::left << fineResidual[i] << std::setw(20) << std::left << dualWeightedResidual[i] << std::endl;
    }
    std::cout << std::endl << "initial_rom_to_final_rom_error: " << initial_rom_to_final_rom_error << std::endl;
}

template <int dim, int nstate>
void Snapshot<dim, nstate>::compute_total_error(){
    total_error = fom_to_initial_rom_error + initial_rom_to_final_rom_error;
}


template class Snapshot <PHILIP_DIM, 1>;
template class Snapshot <PHILIP_DIM, 2>;
template class Snapshot <PHILIP_DIM, 3>;
template class Snapshot <PHILIP_DIM, 4>;
template class Snapshot <PHILIP_DIM, 5>;

}
}
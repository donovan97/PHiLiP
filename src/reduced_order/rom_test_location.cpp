#include "rom_test_location.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template <int dim, int nstate>
ROMTestLocation<dim, nstate>::ROMTestLocation(RowVector2d parameter, std::shared_ptr<ROMSolution<dim, nstate>> rom_solution)
        : parameter(parameter)
        , rom_solution(rom_solution)
{
    std::cout << "Creating ROM test location..." << std::endl;
    compute_FOM_to_initial_ROM_error();
    initial_rom_to_final_rom_error = 0;
    compute_total_error();
    std::cout << "ROM test location created. Error estimate updated." << std::endl;
}

template <int dim, int nstate>
void ROMTestLocation<dim, nstate>::compute_FOM_to_initial_ROM_error(){
    std::cout << "Computing adjoint-based error estimate between ROM and FOM..." << std::endl;
    dealii::LinearAlgebra::distributed::Vector<double> gradient(rom_solution->right_hand_side);
    dealii::LinearAlgebra::distributed::Vector<double> adjoint(rom_solution->right_hand_side);
    adjoint.update_ghost_values();
    gradient.update_ghost_values();

    gradient = rom_solution->gradient;

    Parameters::LinearSolverParam linear_solver_param;
    linear_solver_param.linear_solver_type = Parameters::LinearSolverParam::direct;
    solve_linear(*rom_solution->system_matrix_transpose, gradient*=-1.0, adjoint, linear_solver_param);

    //Compute dual weighted residual
    fom_to_initial_rom_error = 0;
    fom_to_initial_rom_error = -(adjoint * rom_solution->right_hand_side);
    fom_to_initial_rom_error  = dealii::Utilities::MPI::sum(fom_to_initial_rom_error, MPI_COMM_WORLD);
    std::cout << "Parameter: " << parameter << ". Error estimate between ROM and FOM: " << fom_to_initial_rom_error << std::endl;
}

template <int dim, int nstate>
void ROMTestLocation<dim, nstate>::compute_initial_rom_to_final_rom_error(std::shared_ptr<ProperOrthogonalDecomposition::POD<dim>> pod_updated){
    std::cout << "Computing adjoint-based error estimate between initial ROM and updated ROM..." << std::endl;

    dealii::LinearAlgebra::distributed::Vector<double> fineGradient(pod_updated->getPODBasis()->locally_owned_range_indices(), MPI_COMM_WORLD);
    dealii::LinearAlgebra::distributed::Vector<double> fineResidual(pod_updated->getPODBasis()->locally_owned_range_indices(), MPI_COMM_WORLD);

    pod_updated->getPODBasis()->Tvmult(fineGradient, rom_solution->gradient);

    //Compute fine jacobian transpose (Petrov-Galerkin)
    dealii::TrilinosWrappers::SparseMatrix petrov_galerkin_basis;
    dealii::TrilinosWrappers::SparseMatrix fineJacobianTranspose;
    rom_solution->system_matrix_transpose->Tmmult(petrov_galerkin_basis, *pod_updated->getPODBasis()); // petrov_galerkin_basis = system_matrix * pod_basis. Note, use transpose in subsequent multiplications
    petrov_galerkin_basis.Tmmult(fineJacobianTranspose, petrov_galerkin_basis); //reduced_lhs = petrov_galerkin_basis^T * petrov_galerkin_basis , equivalent to V^T*J^T*J*V

    dealii::TrilinosWrappers::SparseMatrix::iterator it = fineJacobianTranspose.begin();
    dealii::TrilinosWrappers::SparseMatrix::iterator it_end = fineJacobianTranspose.end();
    Eigen::SparseMatrix<double, Eigen::RowMajor> eigenFineJacobianTranspose(fineJacobianTranspose.n(), fineJacobianTranspose.m());

    unsigned int row, col;
    double val;
    int approximate_elements_per_row = fineJacobianTranspose.n_nonzero_elements()/fineJacobianTranspose.m();
    eigenFineJacobianTranspose.reserve(Eigen::VectorXi::Constant(fineJacobianTranspose.m(), approximate_elements_per_row));
    for (; it!=it_end; ++it)
    {
        row = it->row();
        col = it->column();
        val = it->value();
        eigenFineJacobianTranspose.insert(row, col) = val;
    }
    eigenFineJacobianTranspose.makeCompressed();

    std::cout << eigenFineJacobianTranspose << std::endl;

    dealii::LinearAlgebra::ReadWriteVector<double> readFineGradient(fineGradient.size());
    readFineGradient.import(fineGradient, dealii::VectorOperation::values::insert);
    Eigen::VectorXd eigenFineGradient(fineGradient.size());
    for(unsigned int i = 0 ; i < fineGradient.size() ; i++){
        eigenFineGradient(i) = readFineGradient(i);
    }

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> householder(eigenFineJacobianTranspose);
    Eigen::VectorXd eigenFineAdjoint = householder.solve(-eigenFineGradient);

    std::cout << eigenFineAdjoint << std::endl;

    //Compute fine residual (Petrov_Galerkin)
    petrov_galerkin_basis.Tvmult(fineResidual, rom_solution->right_hand_side);

    dealii::LinearAlgebra::ReadWriteVector<double> readFineResidual(fineResidual.size());
    readFineResidual.import(fineResidual, dealii::VectorOperation::values::insert);
    Eigen::VectorXd eigenFineResidual(fineResidual.size());
    for(unsigned int i = 0 ; i < fineResidual.size() ; i++){
        eigenFineResidual(i) = readFineResidual(i);
    }

    std::cout << eigenFineResidual << std::endl;

    //Compute dual weighted residual
    initial_rom_to_final_rom_error = 0;
    initial_rom_to_final_rom_error = -(eigenFineAdjoint.dot(eigenFineResidual));

    std::cout << "Parameter: " << parameter << ". Error estimate between initial ROM and updated ROM: " << initial_rom_to_final_rom_error << std::endl;
}

template <int dim, int nstate>
void ROMTestLocation<dim, nstate>::compute_total_error(){
    std::cout << "Computing total error estimate between FOM and updated ROM..." << std::endl;
    total_error = fom_to_initial_rom_error - initial_rom_to_final_rom_error;
    std::cout << "Parameter: " << parameter <<  ". Total error estimate between FOM and updated ROM: " << total_error << std::endl;
}


template class ROMTestLocation <PHILIP_DIM, 1>;
template class ROMTestLocation <PHILIP_DIM, 2>;
template class ROMTestLocation <PHILIP_DIM, 3>;
template class ROMTestLocation <PHILIP_DIM, 4>;
template class ROMTestLocation <PHILIP_DIM, 5>;

}
}
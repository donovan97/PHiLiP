#include "pod_petrov_galerkin_ode_solver.h"
#include <deal.II/lac/householder.h>
#include <deal.II/lac/la_parallel_vector.h>

namespace PHiLiP {
namespace ODE {

template <int dim, typename real, typename MeshType>
PODPetrovGalerkinODESolver<dim,real,MeshType>::PODPetrovGalerkinODESolver(std::shared_ptr< DGBase<dim, real, MeshType> > dg_input, std::shared_ptr<ProperOrthogonalDecomposition::POD<dim>> pod)
        : ImplicitODESolver<dim,real,MeshType>(dg_input)
        , pod(pod)
{}

template <int dim, typename real, typename MeshType>
void PODPetrovGalerkinODESolver<dim,real,MeshType>::step_in_time (real dt, const bool pseudotime)
{
    const bool compute_dRdW = true;
    this->dg->assemble_residual(compute_dRdW);
    this->current_time += dt;
    // Solve (M/dt - dRdW) dw = R
    // w = w + dw

    this->dg->system_matrix *= -1.0;

    if (pseudotime) {
        const double CFL = dt;
        this->dg->time_scaled_mass_matrices(CFL);
        this->dg->add_time_scaled_mass_matrices();
    } else {
        this->dg->add_mass_matrices(1.0/dt);
    }


    if ((this->ode_param.ode_output) == Parameters::OutputEnum::verbose &&
        (this->current_iteration%this->ode_param.print_iteration_modulo) == 0 ) {
        this->pcout << " Evaluating system update... " << std::endl;
    }

    /* Reference for Petrov-Galerkin projection: Refer to Equation (23) in the following reference:
    "Efficient non-linear model reduction via a least-squares Petrov–Galerkin projection and compressive tensor approximations"
    Kevin Carlberg, Charbel Bou-Mosleh, Charbel Farhat
    International Journal for Numerical Methods in Engineering, 2011
    */
    //Petrov-Galerkin projection, petrov_galerkin_basis = V^T*J^T, pod basis V, system matrix J
    //V^T*J*V*p = -V^T*R
    this->pcout << "here1" << std::endl;

    dealii::TrilinosWrappers::SparseMatrix::iterator system_matrix_it = this->dg->system_matrix.begin();
    dealii::TrilinosWrappers::SparseMatrix::iterator system_matrix_it_end = this->dg->system_matrix.end();

    Eigen::SparseMatrix<double, Eigen::RowMajor> eigen_system_matrix(this->dg->system_matrix.n(), this->dg->system_matrix.m());

    unsigned int row,col;
    double val;
    int approximate_elements_per_row = this->dg->system_matrix.n_nonzero_elements()/this->dg->system_matrix.m();
    eigen_system_matrix.reserve(Eigen::VectorXi::Constant(this->dg->system_matrix.m(), approximate_elements_per_row));
    for (; system_matrix_it!=system_matrix_it_end; ++system_matrix_it)
    {
        row = system_matrix_it->row();
        col = system_matrix_it->column();
        val = system_matrix_it->value();
        eigen_system_matrix.insert(row, col) = val;
    }
    eigen_system_matrix.makeCompressed();

    Eigen::MatrixXd LHS = eigen_system_matrix * pod->getEigenPODBasis();

    Eigen::MatrixXd LHS2 = LHS.transpose() * LHS;

    this->pcout << LHS2 << std::endl;

    dealii::LinearAlgebra::ReadWriteVector<double> read_right_hand_side(this->dg->right_hand_side.size());
    read_right_hand_side.import(this->dg->right_hand_side, dealii::VectorOperation::values::insert);
    //snapshotVectors.push_back(read_snapshot);
    Eigen::VectorXd eigen_right_hand_side(this->dg->right_hand_side.size());
    for(unsigned int i = 0 ; i < this->dg->right_hand_side.size() ; i++){
        eigen_right_hand_side(i) = read_right_hand_side(i);
    }

    Eigen::MatrixXd eigen_rhs = -pod->getEigenPODBasis().transpose() * eigen_right_hand_side;
    this->pcout << eigen_rhs << std::endl;

    Eigen::HouseholderQR<Eigen::MatrixXd> householder(LHS);
    reduced_solution_update = householder.solve(eigen_right_hand_side);

    this->pcout << reduced_solution_update << std::endl;

    this->pcout << "here2" << std::endl;

    linesearch();

    //FIX THIS this->update_norm = this->solution_update.l2_norm();
    ++(this->current_iteration);
}

template <int dim, typename real, typename MeshType>
double PODPetrovGalerkinODESolver<dim,real,MeshType>::linesearch()
{
    double step_length = 1.0;
    const double step_reduction = 0.5;
    const int maxline = 10;
    const double reduction_tolerance_1 = 1.0;

    const auto old_reduced_solution = reduced_solution;
    const double initial_residual = this->dg->get_residual_l2norm();

    reduced_solution += step_length*reduced_solution_update;

    Eigen::VectorXd new_solution = reference_solution + pod->getEigenPODBasis()*reduced_solution;

    dealii::LinearAlgebra::ReadWriteVector<double> write_solution(this->dg->solution.size());
    for(unsigned int i = 0 ; i < new_solution.size() ; i++){
        write_solution(i) = new_solution(i);
    }
    this->dg->solution.import(write_solution, dealii::VectorOperation::values::insert);

    for(unsigned int i = 0 ; i < pod->getEigenPODBasis().cols() ; i++){
        this->pcout << this->dg->solution(i) << std::endl;
    }

    this->dg->assemble_residual ();
    double new_residual = this->dg->get_residual_l2norm();
    this->pcout << " Step length " << step_length << ". Old residual: " << initial_residual << " New residual: " << new_residual << std::endl;

    int iline = 0;
    for (iline = 0; iline < maxline && new_residual > initial_residual * reduction_tolerance_1; ++iline) {
        step_length = step_length * step_reduction;
        reduced_solution = old_reduced_solution;

        reduced_solution += step_length*reduced_solution_update;
        new_solution = reference_solution + pod->getEigenPODBasis()*reduced_solution;

        write_solution.reinit(this->dg->solution.size());
        for(unsigned int i = 0 ; i < new_solution.size() ; i++){
            write_solution(i) = new_solution(i);
        }
        this->dg->solution.import(write_solution, dealii::VectorOperation::values::insert);

        this->dg->assemble_residual();
        new_residual = this->dg->get_residual_l2norm();
        this->pcout << " Step length " << step_length << " . Old residual: " << initial_residual << " New residual: " << new_residual << std::endl;
    }
    if (iline == 0) this->CFL_factor *= 1.1;

    return step_length;
}

template <int dim, typename real, typename MeshType>
void PODPetrovGalerkinODESolver<dim,real,MeshType>::allocate_ode_system ()
{
    this->pcout << "Allocating ODE system and evaluating mass matrix..." << std::endl;
    reference_solution = this->pod->getReferenceState();

    dealii::LinearAlgebra::ReadWriteVector<double> read_solution(this->dg->solution.size());
    read_solution.import(this->dg->solution, dealii::VectorOperation::values::insert);
    Eigen::VectorXd initial_condition(this->dg->solution.size());
    for(unsigned int i = 0 ; i < this->dg->solution.size() ; i++){
        initial_condition(i) = read_solution(i);
    }

    reduced_solution = pod->getEigenPODBasis().transpose()*(initial_condition - reference_solution);
    Eigen::VectorXd projected_initial_condition = reference_solution + pod->getEigenPODBasis()*reduced_solution;

    dealii::LinearAlgebra::ReadWriteVector<double> write_solution(this->dg->solution.size());
    for(unsigned int i = 0 ; i < this->dg->solution.size() ; i++){
        write_solution(i) = projected_initial_condition(i);
    }
    this->dg->solution.import(write_solution, dealii::VectorOperation::values::insert);

    const bool do_inverse_mass_matrix = false;
    this->dg->evaluate_mass_matrices(do_inverse_mass_matrix);
}

template class PODPetrovGalerkinODESolver<PHILIP_DIM, double, dealii::Triangulation<PHILIP_DIM>>;
template class PODPetrovGalerkinODESolver<PHILIP_DIM, double, dealii::parallel::shared::Triangulation<PHILIP_DIM>>;
#if PHILIP_DIM != 1
template class PODPetrovGalerkinODESolver<PHILIP_DIM, double, dealii::parallel::distributed::Triangulation<PHILIP_DIM>>;
#endif

} // ODE namespace
} // PHiLiP namespace//
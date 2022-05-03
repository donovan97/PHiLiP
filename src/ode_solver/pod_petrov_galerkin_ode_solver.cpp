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
void PODPetrovGalerkinODESolver<dim,real,MeshType>::step_in_time (real dt, const bool /*pseudotime*/)
{
    const bool compute_dRdW = true;
    this->dg->assemble_residual(compute_dRdW);
    this->current_time += dt;
    // Solve (M/dt - dRdW) dw = R
    // w = w + dw

    this->dg->system_matrix *= -1.0;

    this->dg->add_mass_matrices(1.0/dt);

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

    dealii::LinearAlgebra::ReadWriteVector<double> read_right_hand_side(this->dg->right_hand_side.size());
    read_right_hand_side.import(this->dg->right_hand_side, dealii::VectorOperation::values::insert);
    //snapshotVectors.push_back(read_snapshot);
    Eigen::VectorXd eigen_right_hand_side(this->dg->right_hand_side.size());
    //for(unsigned int i = 0 ; i < this->dg->right_hand_side.size() ; i++){
    //    eigen_right_hand_side(i) = read_right_hand_side(i);
    //}

    Eigen::HouseholderQR<Eigen::MatrixXd> householder(LHS);
    Eigen::VectorXd eigen_reduced_solution_update = householder.solve(eigen_right_hand_side);

    dealii::LinearAlgebra::distributed::Vector<double> solution_read_write(this->reduced_solution_update->size());
    for(unsigned int i = 0 ; i < pod->getEigenPODBasis().cols() ; i++){
        solution_read_write(i) = eigen_reduced_solution_update(i);
        //this->pcout << eigen_reduced_solution_update(i) << std::endl;
    }

    *reduced_solution_update = solution_read_write;
    reduced_solution_update->compress(dealii::VectorOperation::insert);

    this->pcout << "here2" << std::endl;

    linesearch();

    this->update_norm = this->solution_update.l2_norm();
    ++(this->current_iteration);
}

template <int dim, typename real, typename MeshType>
double PODPetrovGalerkinODESolver<dim,real,MeshType>::linesearch()
{
    this->pcout << "here3" << std::endl;
    const auto old_reduced_solution = reduced_solution;
    double step_length = 1.0;
    this->pcout << "here4" << std::endl;
    const double step_reduction = 0.5;
    const int maxline = 10;
    const double reduction_tolerance_1 = 1.0;
    this->pcout << "here5" << std::endl;
    const double initial_residual = this->dg->get_residual_l2norm();
    this->pcout << "here6" << std::endl;
    //reduced_solution = old_reduced_solution;
    reduced_solution.add(step_length, *this->reduced_solution_update);
    //for(unsigned int i = 0 ; i < reduced_solution.size() ; i++){
    //    this->pcout << reduced_solution(i) << std::endl;
    //}
    this->pcout << "here7" << std::endl;
    dealii::LinearAlgebra::distributed::Vector<double> new_solution(this->dg->solution.size());
    pod->getPODBasis()->vmult(new_solution, reduced_solution);
    //for(unsigned int i = 0 ; i < new_solution.size() ; i++){
    //    this->pcout << new_solution(i) << std::endl;
    //}
    this->pcout << "here7.2" << std::endl;
    this->dg->solution.copy_locally_owned_data_from(new_solution);
    //this->dg->solution = new_solution;
    this->pcout << "here7.3" << std::endl;

    this->dg->solution.add(1, reference_solution);
    this->pcout << "here8" << std::endl;
    this->dg->assemble_residual ();
    double new_residual = this->dg->get_residual_l2norm();
    this->pcout << " Step length " << step_length << ". Old residual: " << initial_residual << " New residual: " << new_residual << std::endl;
    this->pcout << "here9" << std::endl;
    int iline = 0;
    for (iline = 0; iline < maxline && new_residual > initial_residual * reduction_tolerance_1; ++iline) {
        step_length = step_length * step_reduction;
        reduced_solution = old_reduced_solution;
        reduced_solution.add(step_length, *this->reduced_solution_update);
        new_solution.reinit(this->dg->solution.size());
        pod->getPODBasis()->vmult(new_solution, reduced_solution);
        this->dg->solution.copy_locally_owned_data_from(new_solution);
        //this->dg->solution = new_solution;
        this->dg->solution.add(1, reference_solution);
        this->dg->assemble_residual();
        new_residual = this->dg->get_residual_l2norm();
        this->pcout << " Step length " << step_length << " . Old residual: " << initial_residual << " New residual: " << new_residual << std::endl;
    }
    if (iline == 0) this->CFL_factor *= 2.0;

    return step_length;
}

template <int dim, typename real, typename MeshType>
void PODPetrovGalerkinODESolver<dim,real,MeshType>::allocate_ode_system ()
{
    this->pcout << "Allocating ODE system and evaluating mass matrix..." << std::endl;
    const bool do_inverse_mass_matrix = false;
    this->dg->evaluate_mass_matrices(do_inverse_mass_matrix);

    this->solution_update.reinit(this->dg->right_hand_side);

    reduced_solution_update = std::make_unique<dealii::LinearAlgebra::distributed::Vector<double>>(pod->getEigenPODBasis().cols());
    //reduced_rhs = std::make_unique<dealii::LinearAlgebra::distributed::Vector<double>>(pod->getPODBasis()->n());
    petrov_galerkin_basis = std::make_unique<dealii::TrilinosWrappers::SparseMatrix>(pod->getEigenPODBasis().rows(), pod->getEigenPODBasis().cols(), pod->getPODBasis()->n());
    reduced_lhs = std::make_unique<dealii::TrilinosWrappers::SparseMatrix>();
    reference_solution = this->dg->solution; //Set reference solution to initial conditions
    reduced_solution = dealii::LinearAlgebra::distributed::Vector<double>(pod->getEigenPODBasis().cols()); //Zero if reference solution is the initial conditions
    //reduced_solution.update_ghost_values();
}

template class PODPetrovGalerkinODESolver<PHILIP_DIM, double, dealii::Triangulation<PHILIP_DIM>>;
template class PODPetrovGalerkinODESolver<PHILIP_DIM, double, dealii::parallel::shared::Triangulation<PHILIP_DIM>>;
#if PHILIP_DIM != 1
template class PODPetrovGalerkinODESolver<PHILIP_DIM, double, dealii::parallel::distributed::Triangulation<PHILIP_DIM>>;
#endif

} // ODE namespace
} // PHiLiP namespace//
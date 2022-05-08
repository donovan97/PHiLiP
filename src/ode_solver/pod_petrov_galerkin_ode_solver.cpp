#include "pod_petrov_galerkin_ode_solver.h"
#include <deal.II/lac/la_parallel_vector.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#include "Amesos.h"
#include "Amesos_BaseSolver.h"

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


    Epetra_CrsMatrix *epetra_system_matrix = const_cast<Epetra_CrsMatrix *>(&(this->dg->system_matrix.trilinos_matrix()));
    Epetra_Map system_matrix_rowmap = epetra_system_matrix->RowMap();

    Epetra_CrsMatrix *epetra_pod_basis = const_cast<Epetra_CrsMatrix *>(&(pod->getPODBasis()->trilinos_matrix()));

    Epetra_CrsMatrix epetra_petrov_galerkin_basis(Epetra_DataAccess::Copy, system_matrix_rowmap, pod->getPODBasis()->n());

    EpetraExt::MatrixMatrix::Multiply(*epetra_system_matrix, false, *epetra_pod_basis, false, epetra_petrov_galerkin_basis, false);

    Epetra_MpiComm epetra_comm(MPI_COMM_WORLD);
    Epetra_Map domain_map((int)pod->getPODBasis()->n(), 0, epetra_comm);
    epetra_petrov_galerkin_basis.FillComplete(domain_map, system_matrix_rowmap);

    std::cout << "here" << std::endl;

    Epetra_Vector epetra_right_hand_side(Epetra_DataAccess::View, epetra_system_matrix->RowMap(), this->dg->right_hand_side.begin());

    Epetra_Vector epetra_reduced_rhs(epetra_petrov_galerkin_basis.DomainMap());

    epetra_petrov_galerkin_basis.Multiply(true, epetra_right_hand_side, epetra_reduced_rhs);

    Epetra_CrsMatrix epetra_reduced_lhs(Epetra_DataAccess::Copy, epetra_petrov_galerkin_basis.DomainMap(), pod->getPODBasis()->n());

    EpetraExt::MatrixMatrix::Multiply(epetra_petrov_galerkin_basis, true, epetra_petrov_galerkin_basis, false, epetra_reduced_lhs);

    Epetra_Vector epetra_reduced_solution_update(epetra_reduced_lhs.DomainMap());

    Epetra_LinearProblem linearProblem(&epetra_reduced_lhs, &epetra_reduced_solution_update, &epetra_reduced_rhs);

    Amesos_BaseSolver* Solver;
    Amesos Factory;
    std::string SolverType = "Klu";
    Solver = Factory.Create(SolverType, linearProblem);

    Teuchos::ParameterList List;
    List.set("PrintTiming", true);
    List.set("PrintStatus", true);
    Solver->SetParameters(List);

    this->pcout << "Starting symbolic factorization..." << std::endl;
    Solver->SymbolicFactorization();

    // you can change the matrix values here
    this->pcout << "Starting numeric factorization..." << std::endl;
    Solver->NumericFactorization();

    // you can change LHS and RHS here
    this->pcout << "Starting solution phase..." << std::endl;
    Solver->Solve();

    //delete Solver;
    //delete epetra_pod_basis;
    //delete epetra_system_matrix;

    this->pcout << "here4" << std::endl;

    linesearch();

    this->update_norm = this->solution_update.l2_norm();
    ++(this->current_iteration);
}

template <int dim, typename real, typename MeshType>
double PODPetrovGalerkinODESolver<dim,real,MeshType>::linesearch()
{
    const auto old_reduced_solution = reduced_solution;
    double step_length = 1.0;

    const double step_reduction = 0.5;
    const int maxline = 10;
    const double reduction_tolerance_1 = 1.0;

    const double initial_residual = this->dg->get_residual_l2norm();

    reduced_solution = old_reduced_solution;
    reduced_solution.add(step_length, *this->reduced_solution_update);
    pod->getPODBasis()->vmult(this->dg->solution, reduced_solution);
    this->dg->solution.add(1, reference_solution);

    this->dg->assemble_residual ();
    double new_residual = this->dg->get_residual_l2norm();
    this->pcout << " Step length " << step_length << ". Old residual: " << initial_residual << " New residual: " << new_residual << std::endl;

    int iline = 0;
    for (iline = 0; iline < maxline && new_residual > initial_residual * reduction_tolerance_1; ++iline) {
        step_length = step_length * step_reduction;
        reduced_solution = old_reduced_solution;
        reduced_solution.add(step_length, *this->reduced_solution_update);
        pod->getPODBasis()->vmult(this->dg->solution, reduced_solution);
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

    reduced_solution_update = std::make_unique<dealii::LinearAlgebra::distributed::Vector<double>>(pod->getPODBasis()->n());
    reduced_rhs = std::make_unique<dealii::LinearAlgebra::distributed::Vector<double>>(pod->getPODBasis()->n());
    petrov_galerkin_basis = std::make_unique<dealii::TrilinosWrappers::SparseMatrix>(pod->getPODBasis()->m(), pod->getPODBasis()->n(), pod->getPODBasis()->n());
    reduced_lhs = std::make_unique<dealii::TrilinosWrappers::SparseMatrix>();
    reference_solution = this->dg->solution; //Set reference solution to initial conditions
    reduced_solution = dealii::LinearAlgebra::distributed::Vector<double>(pod->getPODBasis()->n()); //Zero if reference solution is the initial conditions
}

template class PODPetrovGalerkinODESolver<PHILIP_DIM, double, dealii::Triangulation<PHILIP_DIM>>;
template class PODPetrovGalerkinODESolver<PHILIP_DIM, double, dealii::parallel::shared::Triangulation<PHILIP_DIM>>;
#if PHILIP_DIM != 1
template class PODPetrovGalerkinODESolver<PHILIP_DIM, double, dealii::parallel::distributed::Triangulation<PHILIP_DIM>>;
#endif

} // ODE namespace
} // PHiLiP namespace//
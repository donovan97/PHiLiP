#include "pod_petrov_galerkin_adaptation.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

using DealiiVector = dealii::LinearAlgebra::distributed::Vector<double>;

template<int dim, int nstate>
PODPetrovGalerkinAdaptation<dim, nstate>::PODPetrovGalerkinAdaptation(std::shared_ptr<DGBase<dim, double>> &dg_input, Functional<dim, nstate, double> &functional_input)
        : PODAdaptation<dim, nstate>(dg_input, functional_input)
        {}


template <int dim, int nstate>
void PODPetrovGalerkinAdaptation<dim, nstate>::getDualWeightedResidual()
{
    //Initialize
    DealiiVector fineGradient(this->finePOD->getPODBasis()->n());
    DealiiVector fineAdjoint(this->finePOD->getPODBasis()->n());
    std::vector<double> fineNotInCoarseAdjoint(this->fineNotInCoarsePOD->getPODBasis()->n());
    DealiiVector fineNotInCoarseResidual(this->fineNotInCoarsePOD->getPODBasis()->n());
    this->dualWeightedResidual.reinit(this->fineNotInCoarsePOD->getPODBasis()->n());

    //Compute coarse solution
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::pod_petrov_galerkin_solver;
    this->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, this->dg, this->coarsePOD);
    this->ode_solver->steady_state();

    //Output coarse functional
    this->pcout << "Coarse functional: " << std::setprecision(15) << this->functional.evaluate_functional(false,false) << std::setprecision(6) << std::endl;

    //Compute fine adjoint
    this->pcout << "HERE1" << std::endl;
    this->getReducedGradient(fineGradient);
    this->pcout << "HERE2" << std::endl;
    applyReducedJacobianTranspose(fineAdjoint, fineGradient);

    this->pcout << "HERE3" << std::endl;
    for(unsigned int i = 0 ; i < this->fineNotInCoarsePOD->getFullBasisIndices().size() ; i++){
        this->pcout << fineAdjoint[i] << std::endl;
        this->pcout << this->fineNotInCoarsePOD->getFullBasisIndices()[i] << std::endl;
    }
    //Extract fine not in coarse adjoint
    fineAdjoint.extract_subvector_to(this->fineNotInCoarsePOD->getFullBasisIndices(), fineNotInCoarseAdjoint);
    this->pcout << "HERE4" << std::endl;

    //Compute fine not in coarse residual
    dealii::TrilinosWrappers::SparseMatrix petrov_galerkin_basis;
    this->dg->system_matrix.mmult(petrov_galerkin_basis, *this->fineNotInCoarsePOD->getPODBasis()); // petrov_galerkin_basis= system_matrix * pod_basis.
    petrov_galerkin_basis.Tvmult(fineNotInCoarseResidual, this->dg->right_hand_side);

    //Compute dual weighted residual
    this->error = 0;
    this->pcout << std::setw(10) << std::left << "Index" << std::setw(20) << std::left << "Reduced Adjoint" << std::setw(20) << std::left << "Reduced Residual" << std::setw(20) << std::left << "Dual Weighted Residual" << std::endl;
    for(unsigned int i = 0; i < fineNotInCoarseAdjoint.size(); i++){
        this->dualWeightedResidual[i] = -(fineNotInCoarseAdjoint[i] * fineNotInCoarseResidual[i]);
        this->error = this->error + this->dualWeightedResidual[i];
        this->pcout << std::setw(10) << std::left << this->fineNotInCoarsePOD->getFullBasisIndices()[i] << std::setw(20) << std::left << fineNotInCoarseAdjoint[i] << std::setw(20) << std::left << fineNotInCoarseResidual[i] << std::setw(20) << std::left << this->dualWeightedResidual[i] << std::endl;
    }
    this->pcout << std::endl << "Total error: " << this->error << std::endl;
}

template <int dim, int nstate>
void PODPetrovGalerkinAdaptation<dim, nstate>::applyReducedJacobianTranspose(DealiiVector &reducedAdjoint, DealiiVector &reducedGradient)
{
    const bool compute_dRdW=true;
    const bool compute_dRdX=false;
    const bool compute_d2R=false;
    double flow_CFL_ = 0.0;
    this->dg->assemble_residual(compute_dRdW, compute_dRdX, compute_d2R, flow_CFL_);

    dealii::TrilinosWrappers::SparseMatrix petrov_galerkin_basis;
    dealii::TrilinosWrappers::SparseMatrix reducedJacobianTranspose;

    this->dg->system_matrix.mmult(petrov_galerkin_basis, *this->finePOD->getPODBasis()); // petrov_galerkin_basis = system_matrix * pod_basis. Note, use transpose in subsequent multiplications
    petrov_galerkin_basis.Tmmult(reducedJacobianTranspose, petrov_galerkin_basis); //reduced_lhs = petrov_galerkin_basis^T * petrov_galerkin_basis , equivalent to V^T*J^T*J*V

    solve_linear (reducedJacobianTranspose, reducedGradient*=-1.0, reducedAdjoint, this->linear_solver_param);
}

template class PODPetrovGalerkinAdaptation <PHILIP_DIM,1>;
template class PODPetrovGalerkinAdaptation <PHILIP_DIM,2>;
template class PODPetrovGalerkinAdaptation <PHILIP_DIM,3>;
template class PODPetrovGalerkinAdaptation <PHILIP_DIM,4>;
template class PODPetrovGalerkinAdaptation <PHILIP_DIM,5>;

}
}

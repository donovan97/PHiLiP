#include "pod_full_dim_adaptation.h"


namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

using DealiiVector = dealii::LinearAlgebra::distributed::Vector<double>;

/******************************************************
* PODFullDimAdaptation
******************************************************/

template<int dim, int nstate>
PODFullDimAdaptation<dim, nstate>::PODFullDimAdaptation(std::shared_ptr<DGBase<dim, double>> &dg_input,
                                                            Functional<dim, nstate, double> &functional_input)
        : PODAdaptation<dim, nstate>(dg_input, functional_input)
        {}

template<int dim, int nstate>
std::vector<unsigned int> PODFullDimAdaptation<dim, nstate>::getPODBasisColumnsToAdd() {
    std::map<double, unsigned int> dualWeightedResidualToIndex;
    std::vector<unsigned int> PODBasisColumnsToAdd;
    std::map<double, unsigned int>::iterator element;

    this->pcout << std::setw(10) << std::left << "Index" << std::setw(20) << std::left << "Dual-weighted residual" << std::endl;
    for(unsigned int i = 0; i < this->fineNotInCoarseDualWeightedResidual.size(); i++){
        dualWeightedResidualToIndex.emplace(abs(fineNotInCoarseDualWeightedResidual[i]), this->fineNotInCoarsePOD->getFullBasisIndices()[i]);
        this->pcout << std::setw(10) << std::left << this->fineNotInCoarsePOD->getFullBasisIndices()[i] << std::setw(20) << this->fineNotInCoarseDualWeightedResidual[i] << std::endl;
    }

    for (unsigned int i = 0; i < this->all_parameters->reduced_order_param.adapt_coarse_basis_constant; i++) { //Add user-specified number of basis vectors
        element = std::prev(dualWeightedResidualToIndex.end());
        PODBasisColumnsToAdd.push_back(element->second);
        this->pcout << "Adding POD basis: " << element->second << std::endl;
        dualWeightedResidualToIndex.erase(element);
    }

    return PODBasisColumnsToAdd;
}

template<int dim, int nstate>
void PODFullDimAdaptation<dim, nstate>::applyJacobianTranspose(DealiiVector &adjoint, DealiiVector &gradient) {
    const bool compute_dRdW=true;
    const bool compute_dRdX=false;
    const bool compute_d2R=false;
    double flow_CFL_ = 0.0;
    this->dg->assemble_residual(compute_dRdW, compute_dRdX, compute_d2R, flow_CFL_);

    solve_linear (this->dg->system_matrix_transpose, gradient*=-1.0, adjoint, this->linear_solver_param);
}

template<int dim, int nstate>
void PODFullDimAdaptation<dim, nstate>::getGradient(DealiiVector &gradient) {
    this->functional.set_state(this->dg->solution);
    this->functional.dg->high_order_grid->volume_nodes = this->dg->high_order_grid->volume_nodes;

    const bool compute_dIdW = true;
    const bool compute_dIdX = false;
    const bool compute_d2I = false;
    this->functional.evaluate_functional( compute_dIdW, compute_dIdX, compute_d2I );

    gradient = this->functional.dIdw;
}

/******************************************************
* PODPetrovGalerkinFullDimAdaptation
******************************************************/

template<int dim, int nstate>
PODPetrovGalerkinFullDimAdaptation<dim, nstate>::PODPetrovGalerkinFullDimAdaptation(
        std::shared_ptr<DGBase<dim, double>> &dg_input, Functional<dim, nstate, double> &functional_input)
        : PODFullDimAdaptation<dim, nstate>(dg_input, functional_input) {}


template<int dim, int nstate>
void PODPetrovGalerkinFullDimAdaptation<dim, nstate>::getDualWeightedResidual() {
    //Initialize
    DealiiVector gradient(this->fineNotInCoarsePOD->getPODBasis()->m());
    DealiiVector adjoint(this->fineNotInCoarsePOD->getPODBasis()->m());
    DealiiVector fineResidual(this->fineNotInCoarsePOD->getPODBasis()->m());
    this->dualWeightedResidual.reinit(this->fineNotInCoarsePOD->getPODBasis()->m());

    //Compute coarse solution
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::pod_petrov_galerkin_solver;
    this->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, this->dg, this->coarsePOD);
    this->ode_solver->steady_state();

    //Output coarse functional
    this->pcout << "Coarse functional: " << std::setprecision(15) << this->functional.evaluate_functional(false,false) << std::setprecision(6) << std::endl;

    //Compute fine adjoint
    this->pcout << "HERE1" << std::endl;
    this->getGradient(gradient);
    this->pcout << "HERE2" << std::endl;
    this->applyJacobianTranspose(adjoint, gradient);

    this->pcout << "HERE3" << std::endl;
    for(unsigned int i = 0 ; i < adjoint.size() ; i++){
        //pcout << adjoint[i] << std::endl;
    }
    //Extract fine not in coarse adjoint
    //fineAdjoint.extract_subvector_to(fineNotInCoarsePOD->getFullBasisIndices(), fineNotInCoarseAdjoint);
    //pcout << "HERE4" << std::endl;

    //Compute fine residual
    //coarsePOD->getPODBasis()->vmult(fineResidual, dg->right_hand_side);

    //Compute dual weighted residual
    this->error = 0;
    for(unsigned int i = 0; i < adjoint.size(); i++){
        this->dualWeightedResidual[i] = -(adjoint[i] * this->dg->right_hand_side[i]);
        this->error = this->error + this->dualWeightedResidual[i];
        //this->pcout << std::setw(10) << std::left << this->fineNotInCoarsePOD->getFullBasisIndices()[i] << std::setw(20) << std::left << adjoint[i] << std::setw(20) << std::left << this->dg->right_hand_side[i] << std::setw(20) << std::left << this->dualWeightedResidual[i] << std::endl;
    }
    this->pcout << std::endl << "Total error: " << this->error << std::endl;

    this->fineNotInCoarseDualWeightedResidual.reinit(this->fineNotInCoarsePOD->getPODBasis()->n());
    this->fineNotInCoarsePOD->getPODBasis()->Tvmult(this->fineNotInCoarseDualWeightedResidual, this->dualWeightedResidual);
}

/******************************************************
* PODGalerkinFullDimAdaptation
******************************************************/

template<int dim, int nstate>
PODGalerkinFullDimAdaptation<dim, nstate>::PODGalerkinFullDimAdaptation(
        std::shared_ptr<DGBase<dim, double>> &dg_input, Functional<dim, nstate, double> &functional_input)
        : PODFullDimAdaptation<dim, nstate>(dg_input, functional_input) {}


template<int dim, int nstate>
void PODGalerkinFullDimAdaptation<dim, nstate>::getDualWeightedResidual() {
    //Initialize
    DealiiVector gradient(this->fineNotInCoarsePOD->getPODBasis()->m());
    DealiiVector adjoint(this->fineNotInCoarsePOD->getPODBasis()->m());
    //std::vector<double> fineNotInCoarseAdjoint(fineNotInCoarsePOD->getPODBasis()->n());
    DealiiVector fineResidual(this->fineNotInCoarsePOD->getPODBasis()->m());
    this->dualWeightedResidual.reinit(this->fineNotInCoarsePOD->getPODBasis()->m());

    //Compute coarse solution
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::pod_galerkin_solver;
    this->ode_solver =  PHiLiP::ODE::ODESolverFactory<dim, double>::create_ODESolver_manual(ode_solver_type, this->dg, this->coarsePOD);
    this->ode_solver->steady_state();

    //Output coarse functional
    this->pcout << "Coarse functional: " << std::setprecision(15) << this->functional.evaluate_functional(false,false) << std::setprecision(6) << std::endl;

    //Compute fine adjoint
    this->pcout << "HERE1" << std::endl;
    this->getGradient(gradient);
    this->pcout << "HERE2" << std::endl;
    this->applyJacobianTranspose(adjoint, gradient);

    this->pcout << "HERE3" << std::endl;
    for(unsigned int i = 0 ; i < adjoint.size() ; i++){
        //pcout << adjoint[i] << std::endl;
    }
    //Extract fine not in coarse adjoint
    //fineAdjoint.extract_subvector_to(fineNotInCoarsePOD->getFullBasisIndices(), fineNotInCoarseAdjoint);
    //pcout << "HERE4" << std::endl;

    //Compute fine residual
    //coarsePOD->getPODBasis()->vmult(fineResidual, dg->right_hand_side);

    //Compute dual weighted residual
    this->error = 0;
    //pcout << std::setw(10) << std::left << "Index" << std::setw(20) << std::left << "Reduced Adjoint" << std::setw(20) << std::left << "Reduced Residual" << std::setw(20) << std::left << "Dual Weighted Residual" << std::endl;
    for(unsigned int i = 0; i < adjoint.size(); i++){
        this->dualWeightedResidual[i] = -(adjoint[i] * this->dg->right_hand_side[i]);
        this->error = this->error + this->dualWeightedResidual[i];
        //pcout << std::setw(10) << std::left << fineNotInCoarsePOD->getFullBasisIndices()[i] << std::setw(20) << std::left << fineNotInCoarseAdjoint[i] << std::setw(20) << std::left << fineNotInCoarseResidual[i] << std::setw(20) << std::left << dualWeightedResidual[i] << std::endl;
    }
    this->pcout << std::endl << "Total error: " << this->error << std::endl;

    dealii::TrilinosWrappers::SparseMatrix petrov_galerkin_basis;
    this->fineNotInCoarseDualWeightedResidual.reinit(this->fineNotInCoarsePOD->getPODBasis()->n());

    this->dg->system_matrix.mmult(petrov_galerkin_basis, *this->fineNotInCoarsePOD->getPODBasis()); // petrov_galerkin_basis = system_matrix * pod_basis. Note, use transpose in subsequent multiplications
    petrov_galerkin_basis.Tvmult(this->fineNotInCoarseDualWeightedResidual, this->dualWeightedResidual);
}

template class PODFullDimAdaptation <PHILIP_DIM,1>;
template class PODFullDimAdaptation <PHILIP_DIM,2>;
template class PODFullDimAdaptation <PHILIP_DIM,3>;
template class PODFullDimAdaptation <PHILIP_DIM,4>;
template class PODFullDimAdaptation <PHILIP_DIM,5>;

template class PODPetrovGalerkinFullDimAdaptation <PHILIP_DIM,1>;
template class PODPetrovGalerkinFullDimAdaptation <PHILIP_DIM,2>;
template class PODPetrovGalerkinFullDimAdaptation <PHILIP_DIM,3>;
template class PODPetrovGalerkinFullDimAdaptation <PHILIP_DIM,4>;
template class PODPetrovGalerkinFullDimAdaptation <PHILIP_DIM,5>;

template class PODGalerkinFullDimAdaptation <PHILIP_DIM,1>;
template class PODGalerkinFullDimAdaptation <PHILIP_DIM,2>;
template class PODGalerkinFullDimAdaptation <PHILIP_DIM,3>;
template class PODGalerkinFullDimAdaptation <PHILIP_DIM,4>;
template class PODGalerkinFullDimAdaptation <PHILIP_DIM,5>;
}
}
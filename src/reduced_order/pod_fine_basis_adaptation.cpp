#include "pod_fine_basis_adaptation.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

using DealiiVector = dealii::LinearAlgebra::distributed::Vector<double>;

/******************************************************
 * PODFineBasisAdaptation
 ******************************************************/

template<int dim, int nstate>
PODFineBasisAdaptation<dim, nstate>::PODFineBasisAdaptation(std::shared_ptr<DGBase<dim, double>> &dg_input, Functional<dim, nstate, double> &functional_input)
        : PODAdaptation<dim, nstate>(dg_input, functional_input)
        , finePOD(std::make_shared<ProperOrthogonalDecomposition::FineStatePOD<dim>>(dg_input))
{}

template <int dim, int nstate>
std::vector<unsigned int> PODFineBasisAdaptation<dim, nstate>::getPODBasisColumnsToAdd()
{
    std::map<double, unsigned int> dualWeightedResidualToIndex;
    std::vector<unsigned int> PODBasisColumnsToAdd;
    std::map<double, unsigned int>::iterator element;

    if(this->all_parameters->reduced_order_param.consider_error_sign){ //if considering sign of error
        for(unsigned int i = 0; i < this->dualWeightedResidual.size(); i++){
            dualWeightedResidualToIndex.emplace(this->dualWeightedResidual[i], this->fineNotInCoarsePOD->getFullBasisIndices()[i]); //will automatically sort
        }
        double adaptationError = this->error;
        if(this->all_parameters->reduced_order_param.adapt_coarse_basis_constant == 0){ //Automatically choose how many basis vectors to add
            while(abs(adaptationError) > this->all_parameters->reduced_order_param.adaptation_tolerance){
                if(adaptationError > 0){
                    element = std::prev(dualWeightedResidualToIndex.end()); //Get largest negative value
                }
                else{
                    element = dualWeightedResidualToIndex.begin(); //Get largest positive value
                }
                PODBasisColumnsToAdd.push_back(element->second);
                this->pcout << "Adding POD basis: " << element->second << std::endl;
                adaptationError = adaptationError - element->first;
                this->pcout << "Estimated adaptation error: " << adaptationError << std::endl;
                dualWeightedResidualToIndex.erase(element);
            }
        }
        else{
            for (unsigned int i = 0; i < this->all_parameters->reduced_order_param.adapt_coarse_basis_constant; i++) { //Add user-specified number of basis vectors
                if(abs(adaptationError) < this->all_parameters->reduced_order_param.adaptation_tolerance){
                    break;
                }
                if(adaptationError > 0){
                    element = std::prev(dualWeightedResidualToIndex.end()); //Get largest negative value
                }
                else{
                    element = dualWeightedResidualToIndex.begin(); //Get largest positive value
                }
                PODBasisColumnsToAdd.push_back(element->second);
                this->pcout << "Adding POD basis: " << element->second << std::endl;
                adaptationError = adaptationError - element->first;
                this->pcout << "Estimated adaptation error: " << adaptationError << std::endl;
                dualWeightedResidualToIndex.erase(element);
            }
        }
    }
    else{ //If consdering only the absolute value of errors
        for(unsigned int i = 0; i < this->dualWeightedResidual.size(); i++){
            dualWeightedResidualToIndex.emplace(abs(this->dualWeightedResidual[i]), this->fineNotInCoarsePOD->getFullBasisIndices()[i]);
        }
        for (unsigned int i = 0; i < this->all_parameters->reduced_order_param.adapt_coarse_basis_constant; i++) { //Add user-specified number of basis vectors
            element = std::prev(dualWeightedResidualToIndex.end());
            PODBasisColumnsToAdd.push_back(element->second);
            this->pcout << "Adding POD basis: " << element->second << std::endl;
            dualWeightedResidualToIndex.erase(element);
        }
    }

    return PODBasisColumnsToAdd;
}

template <int dim, int nstate>
void PODFineBasisAdaptation<dim, nstate>::getReducedGradient(DealiiVector &reducedGradient)
{
    this->functional.set_state(this->dg->solution);
    this->functional.dg->high_order_grid->volume_nodes = this->dg->high_order_grid->volume_nodes;

    const bool compute_dIdW = true;
    const bool compute_dIdX = false;
    const bool compute_d2I = false;
    this->functional.evaluate_functional( compute_dIdW, compute_dIdX, compute_d2I );

    finePOD->getPODBasis()->Tvmult(reducedGradient, this->functional.dIdw); // reducedGradient= (pod_basis)^T * gradient
}

/******************************************************
 * PODPetrovGalerkinFineBasisAdaptation
 ******************************************************/

template<int dim, int nstate>
PODPetrovGalerkinFineBasisAdaptation<dim, nstate>::PODPetrovGalerkinFineBasisAdaptation(std::shared_ptr<DGBase<dim, double>> &dg_input, Functional<dim, nstate, double> &functional_input)
        : PODFineBasisAdaptation<dim, nstate>(dg_input, functional_input)
        {}


template <int dim, int nstate>
void PODPetrovGalerkinFineBasisAdaptation<dim, nstate>::getDualWeightedResidual()
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
void PODPetrovGalerkinFineBasisAdaptation<dim, nstate>::applyReducedJacobianTranspose(DealiiVector &reducedAdjoint, DealiiVector &reducedGradient)
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

/******************************************************
 * PODGalerkinFineBasisAdaptation
 ******************************************************/

template<int dim, int nstate>
PODGalerkinFineBasisAdaptation<dim, nstate>::PODGalerkinFineBasisAdaptation(std::shared_ptr<DGBase<dim, double>> &dg_input, Functional<dim, nstate, double> &functional_input)
        : PODFineBasisAdaptation<dim, nstate>(dg_input, functional_input)
        {}


template <int dim, int nstate>
void PODGalerkinFineBasisAdaptation<dim, nstate>::getDualWeightedResidual()
{
    //Initialize
    DealiiVector fineGradient(this->finePOD->getPODBasis()->n());
    DealiiVector fineAdjoint(this->finePOD->getPODBasis()->n());
    std::vector<double> fineNotInCoarseAdjoint(this->fineNotInCoarsePOD->getPODBasis()->n());
    DealiiVector fineNotInCoarseResidual(this->fineNotInCoarsePOD->getPODBasis()->n());
    this->dualWeightedResidual.reinit(this->fineNotInCoarsePOD->getPODBasis()->n());

    //Compute coarse solution
    auto ode_solver_type = Parameters::ODESolverParam::ODESolverEnum::pod_galerkin_solver;
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
    this->fineNotInCoarsePOD->getPODBasis()->Tvmult(fineNotInCoarseResidual, this->dg->right_hand_side);

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
void PODGalerkinFineBasisAdaptation<dim, nstate>::applyReducedJacobianTranspose(DealiiVector &reducedAdjoint, DealiiVector &reducedGradient)
{
    const bool compute_dRdW=true;
    const bool compute_dRdX=false;
    const bool compute_d2R=false;
    double flow_CFL_ = 0.0;
    this->dg->assemble_residual(compute_dRdW, compute_dRdX, compute_d2R, flow_CFL_);

    dealii::TrilinosWrappers::SparseMatrix tmp;
    dealii::TrilinosWrappers::SparseMatrix reducedJacobianTranspose;
    this->finePOD->getPODBasis()->Tmmult(tmp, this->dg->system_matrix_transpose); //tmp = pod_basis^T * dg->system_matrix_transpose
    tmp.mmult(reducedJacobianTranspose, *this->finePOD->getPODBasis()); // reducedJacobianTranspose= tmp*pod_basis

    solve_linear (reducedJacobianTranspose, reducedGradient*=-1.0, reducedAdjoint, this->linear_solver_param);
}

template class PODFineBasisAdaptation <PHILIP_DIM,1>;
template class PODFineBasisAdaptation <PHILIP_DIM,2>;
template class PODFineBasisAdaptation <PHILIP_DIM,3>;
template class PODFineBasisAdaptation <PHILIP_DIM,4>;
template class PODFineBasisAdaptation <PHILIP_DIM,5>;

template class PODPetrovGalerkinFineBasisAdaptation <PHILIP_DIM,1>;
template class PODPetrovGalerkinFineBasisAdaptation <PHILIP_DIM,2>;
template class PODPetrovGalerkinFineBasisAdaptation <PHILIP_DIM,3>;
template class PODPetrovGalerkinFineBasisAdaptation <PHILIP_DIM,4>;
template class PODPetrovGalerkinFineBasisAdaptation <PHILIP_DIM,5>;

template class PODGalerkinFineBasisAdaptation <PHILIP_DIM,1>;
template class PODGalerkinFineBasisAdaptation <PHILIP_DIM,2>;
template class PODGalerkinFineBasisAdaptation <PHILIP_DIM,3>;
template class PODGalerkinFineBasisAdaptation <PHILIP_DIM,4>;
template class PODGalerkinFineBasisAdaptation <PHILIP_DIM,5>;

}
}

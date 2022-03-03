#include "pod_adaptation_base.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

using DealiiVector = dealii::LinearAlgebra::distributed::Vector<double>;

template <int dim, int nstate>
PODAdaptation<dim, nstate>::PODAdaptation(std::shared_ptr<DGBase<dim,double>> &dg_input, Functional<dim,nstate,double> &functional_input)
        : functional(functional_input)
        , dg(dg_input)
        , all_parameters(dg->all_parameters)
        , coarsePOD(std::make_shared<ProperOrthogonalDecomposition::CoarseStatePOD<dim>>(dg))
        , fineNotInCoarsePOD(std::make_unique<ProperOrthogonalDecomposition::FineNotInCoarseStatePOD<dim>>(dg))
        , mpi_communicator(MPI_COMM_WORLD)
        , pcout(std::cout, dealii::Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
{
    dealii::ParameterHandler parameter_handler;
    Parameters::LinearSolverParam::declare_parameters (parameter_handler);
    this->linear_solver_param.parse_parameters (parameter_handler);
    linear_solver_param.linear_solver_type = Parameters::LinearSolverParam::direct;
}

template <int dim, int nstate>
void PODAdaptation<dim, nstate>::progressivePODAdaptation()
{
    getDualWeightedResidual();
    while(abs(error) > all_parameters->reduced_order_param.adaptation_tolerance){
        std::vector<unsigned int> newColumns = getPODBasisColumnsToAdd();
        coarsePOD->addPODBasisColumns(newColumns);
        fineNotInCoarsePOD->removePODBasisColumns(newColumns);
        getDualWeightedResidual();

        if(fineNotInCoarsePOD->getFullBasisIndices().empty()){
            pcout << "No basis vectors remaining to add!" << std::endl;
            break;
        }
    }
    pcout << "Error estimate is smaller than desired tolerance!" << std::endl;
}


template <int dim, int nstate>
void PODAdaptation<dim, nstate>::simplePODAdaptation()
{
    getDualWeightedResidual();
    if(abs(error) > all_parameters->reduced_order_param.adaptation_tolerance){
        std::vector<unsigned int> newColumns = getPODBasisColumnsToAdd();
        coarsePOD->addPODBasisColumns(newColumns);
        fineNotInCoarsePOD->removePODBasisColumns(newColumns);

        getDualWeightedResidual();
    }
    else{
        pcout << "Error estimate is smaller than desired tolerance!" << std::endl;
    }
}

template <int dim, int nstate>
double PODAdaptation<dim, nstate>::getCoarseFunctional()
{
    return functional.evaluate_functional(false,false);
}

template class PODAdaptation <PHILIP_DIM,1>;
template class PODAdaptation <PHILIP_DIM,2>;
template class PODAdaptation <PHILIP_DIM,3>;
template class PODAdaptation <PHILIP_DIM,4>;
template class PODAdaptation <PHILIP_DIM,5>;

}
}
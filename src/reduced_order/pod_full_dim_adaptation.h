#ifndef __POD_FULL_DIM_ADAPTATION__
#define __POD_FULL_DIM_ADAPTATION__

#include <fstream>
#include <iostream>
#include <filesystem>

#include "functional/functional.h"
#include "dg/dg.h"
#include "reduced_order/pod_basis_types.h"
#include <deal.II/base/function_parser.h>
#include "linear_solver/linear_solver.h"

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include "ode_solver/ode_solver_factory.h"
#include "pod_interfaces.h"
#include "pod_basis_sensitivity_types.h"
#include "pod_adaptation_base.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

using DealiiVector = dealii::LinearAlgebra::distributed::Vector<double>;

template <int dim, int nstate>
class PODFullDimAdaptation: public PODAdaptation<dim, nstate>
{

protected:

    DealiiVector fineNotInCoarseDualWeightedResidual;

    /// Constructor
    PODFullDimAdaptation(std::shared_ptr<DGBase < dim,double>> &dg_input, Functional<dim,nstate,double> &functional_input);

    /// Destructor
    virtual ~PODFullDimAdaptation () {};

    /// Compute reduced-order gradient
    void getGradient(DealiiVector &reducedGradient);

    /// Apply reduced-order Jacobian transpose to solve for reduced-order adjoint
    void applyJacobianTranspose(DealiiVector &reducedAdjoint, DealiiVector &reducedGradient);

    /// Compute dual-weighted residual
    virtual void getDualWeightedResidual() = 0;

    /// Determine which POD basis to add based on dual-weighted residual error
    std::vector<unsigned int> getPODBasisColumnsToAdd();
};

template <int dim, int nstate>
class PODPetrovGalerkinFullDimAdaptation: public PODFullDimAdaptation<dim, nstate>
{

public:
    /// Constructor
    PODPetrovGalerkinFullDimAdaptation(std::shared_ptr<DGBase < dim,double>> &dg_input, Functional<dim,nstate,double> &functional_input);

    /// Destructor
    ~PODPetrovGalerkinFullDimAdaptation () {};

    /// Compute dual-weighted residual
    void getDualWeightedResidual() override;
};

template <int dim, int nstate>
class PODGalerkinFullDimAdaptation: public PODFullDimAdaptation<dim, nstate>
{

public:
    /// Constructor
    PODGalerkinFullDimAdaptation(std::shared_ptr<DGBase < dim,double>> &dg_input, Functional<dim,nstate,double> &functional_input);

    /// Destructor
    ~PODGalerkinFullDimAdaptation () {};

    /// Compute dual-weighted residual
    void getDualWeightedResidual() override;
};

}
}

#endif

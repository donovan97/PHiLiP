#ifndef __POD_FINE_BASIS_ADAPTATION__
#define __POD_FINE_BASIS_ADAPTATION__

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
class PODFineBasisAdaptation: public PODAdaptation<dim, nstate>
{

public:
    /// Smart pointer to fine POD basis
    std::shared_ptr<ProperOrthogonalDecomposition::FineBasis<dim>> finePOD;

public:

    /// Constructor
    PODFineBasisAdaptation(std::shared_ptr<DGBase < dim,double>> &dg_input, Functional<dim,nstate,double> &functional_input);

    /// Destructor
    virtual ~PODFineBasisAdaptation () {};

    /// Compute reduced-order gradient
    void getReducedGradient(DealiiVector &reducedGradient);

    /// Apply reduced-order Jacobian transpose to solve for reduced-order adjoint
    virtual void applyReducedJacobianTranspose(DealiiVector &reducedAdjoint, DealiiVector &reducedGradient) = 0;

    /// Compute dual-weighted residual
    virtual void getDualWeightedResidual() = 0;

    /// Determine which POD basis to add based on dual-weighted residual error
    std::vector<unsigned int> getPODBasisColumnsToAdd();
};

template <int dim, int nstate>
class PODPetrovGalerkinFineBasisAdaptation: public PODFineBasisAdaptation<dim, nstate>
{

public:
    /// Constructor
    PODPetrovGalerkinFineBasisAdaptation(std::shared_ptr<DGBase < dim,double>> &dg_input, Functional<dim,nstate,double> &functional_input);

    /// Destructor
    ~PODPetrovGalerkinFineBasisAdaptation () {};

    /// Apply reduced-order Jacobian transpose to solve for reduced-order adjoint
    void applyReducedJacobianTranspose(DealiiVector &reducedAdjoint, DealiiVector &reducedGradient) override;

    /// Compute dual-weighted residual
    void getDualWeightedResidual() override;
};

template <int dim, int nstate>
class PODGalerkinFineBasisAdaptation: public PODFineBasisAdaptation<dim, nstate>
{

public:
    /// Constructor
    PODGalerkinFineBasisAdaptation(std::shared_ptr<DGBase < dim,double>> &dg_input, Functional<dim,nstate,double> &functional_input);

    /// Destructor
    ~PODGalerkinFineBasisAdaptation () {};

    /// Apply reduced-order Jacobian transpose to solve for reduced-order adjoint
    void applyReducedJacobianTranspose(DealiiVector &reducedAdjoint, DealiiVector &reducedGradient) override;

    /// Compute dual-weighted residual
    void getDualWeightedResidual() override;
};

}
}

#endif

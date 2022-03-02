#ifndef __POD_PETROV_GALERKIN_ADAPTATION__
#define __POD_PETROV_GALERKIN_ADAPTATION__

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
#include "pod_adaptation.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template <int dim, int nstate>
class PODPetrovGalerkinAdaptation: public PODAdaptation<dim, nstate>
{
    using DealiiVector = dealii::LinearAlgebra::distributed::Vector<double>;

public:
    /// Constructor
    PODPetrovGalerkinAdaptation(std::shared_ptr<DGBase<dim,double>> &dg_input, Functional<dim,nstate,double> &functional_input);

    /// Destructor
    ~PODPetrovGalerkinAdaptation () {};

    /// Apply reduced-order Jacobian transpose to solve for reduced-order adjoint
    void applyReducedJacobianTranspose(DealiiVector &reducedAdjoint, DealiiVector &reducedGradient) override;

    /// Compute dual-weighted residual
    void getDualWeightedResidual() override;
};

}
}

#endif

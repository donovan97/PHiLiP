#ifndef __ADAPTIVE_SAMPLING__
#define __ADAPTIVE_SAMPLING__

#include <fstream>
#include <iostream>
#include <filesystem>
#include "functional/functional.h"
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector_operation.h>
#include "parameters/all_parameters.h"
#include "dg/dg.h"
#include "pod_interfaces.h"
#include "reduced_order_solution.h"
#include "full_order_solution.h"
#include "linear_solver/linear_solver.h"
#include "testing/flow_solver.h"
#include "snapshot.h"
#include <cmath>
#include <iostream>

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

/// Class to hold information about the reduced-order solution
template <int dim, int nstate>
class AdaptiveSampling
{
public:
    /// Constructor
    AdaptiveSampling(std::vector<double> parameter_space, double tolerance, const PHiLiP::Parameters::AllParameters *const parameters_input);

    /// Destructor
    virtual ~AdaptiveSampling() {};

    const Parameters::AllParameters *const all_parameters;

    std::vector<double> parameter_space;

    double tolerance;

    std::queue<double> trial_locations;

    std::vector<Snapshot<dim,nstate>> snapshots;

    std::shared_ptr<POD<dim>> current_pod;

    void start();

    void generateTrialLocations(int n);

    void initializeSampling();

    void updateErrors();

    void updateSensitivityCurveFit();

    void updateErrorCurveFit();

    void addSnapshot();

    std::shared_ptr<FOMSolution<dim,nstate>> solveSnapshotFOM(double parameter);

    std::shared_ptr<ROMSolution<dim,nstate>> solveSnapshotROM(double parameter);

    Parameters::AllParameters reinitParams(double parameter);

};

///Functional to take the integral of the solution
template <int dim, int nstate, typename real>
class BurgersRewienskiFunctional : public Functional<dim, nstate, real>
{
public:
    using FadType = Sacado::Fad::DFad<real>; ///< Sacado AD type for first derivatives.
    using FadFadType = Sacado::Fad::DFad<FadType>; ///< Sacado AD type that allows 2nd derivatives.
public:
    /// Constructor
    BurgersRewienskiFunctional(
            std::shared_ptr<PHiLiP::DGBase<dim,real>> dg_input,
            std::shared_ptr<PHiLiP::Physics::PhysicsBase<dim,nstate,FadFadType>> _physics_fad_fad,
            const bool uses_solution_values = true,
            const bool uses_solution_gradient = false)
            : PHiLiP::Functional<dim,nstate,real>(dg_input,_physics_fad_fad,uses_solution_values,uses_solution_gradient)
    {}
    template <typename real2>
    /// Templated volume integrand
    real2 evaluate_volume_integrand(
            const PHiLiP::Physics::PhysicsBase<dim,nstate,real2> &physics,
            const dealii::Point<dim,real2> &phys_coord,
            const std::array<real2,nstate> &soln_at_q,
            const std::array<dealii::Tensor<1,dim,real2>,nstate> &soln_grad_at_q) const;

    /// Non-template functions to override the template classes
    real evaluate_volume_integrand(
            const PHiLiP::Physics::PhysicsBase<dim,nstate,real> &physics,
            const dealii::Point<dim,real> &phys_coord,
            const std::array<real,nstate> &soln_at_q,
            const std::array<dealii::Tensor<1,dim,real>,nstate> &soln_grad_at_q) const override
    {
        return evaluate_volume_integrand<>(physics, phys_coord, soln_at_q, soln_grad_at_q);
    }
    /// Non-template functions to override the template classes
    FadFadType evaluate_volume_integrand(
            const PHiLiP::Physics::PhysicsBase<dim,nstate,FadFadType> &physics,
            const dealii::Point<dim,FadFadType> &phys_coord,
            const std::array<FadFadType,nstate> &soln_at_q,
            const std::array<dealii::Tensor<1,dim,FadFadType>,nstate> &soln_grad_at_q) const override
    {
        return evaluate_volume_integrand<>(physics, phys_coord, soln_at_q, soln_grad_at_q);
    }
};

}
}


#endif

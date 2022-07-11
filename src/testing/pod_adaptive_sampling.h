#ifndef __POD_ADAPTIVE_SAMPLING__
#define __POD_ADAPTIVE_SAMPLING__

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
#include "reduced_order/pod_basis_online.h"
#include "reduced_order/reduced_order_solution.h"
#include "linear_solver/linear_solver.h"
#include "flow_solver/flow_solver.h"
#include "flow_solver/flow_solver_factory.h"
#include "reduced_order/rom_test_location.h"
#include <deal.II/lac/householder.h>
#include <cmath>
#include <iostream>
#include <deal.II/base/function_lib.h>
#include <eigen/Eigen/Dense>
#include "reduced_order/nearest_neighbors.h"
#include "reduced_order/rbf_interpolation.h"
#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Stream.hpp"
#include "ROL_Bounds.hpp"
#include "functional/lift_drag.hpp"
#include "functional/functional.h"
#include "reduced_order/halton.h"
#include "reduced_order/min_max_scaler.h"
#include "tests.h"

namespace PHiLiP {
namespace Tests {

using DealiiVector = dealii::LinearAlgebra::distributed::Vector<double>;
using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::VectorXd;

/// POD adaptive sampling
template <int dim, int nstate>
class AdaptiveSampling: public TestsBase
{
public:
    /// Constructor
    AdaptiveSampling(const PHiLiP::Parameters::AllParameters *const parameters_input,
                     const dealii::ParameterHandler &parameter_handler_input);

    /// Destructor
    ~AdaptiveSampling() {};

    /// Matrix of snapshot parameters
    mutable MatrixXd snapshot_parameters;

    /// Vector of parameter-ROMTestLocation pairs
    mutable std::vector<std::pair<RowVectorXd, std::shared_ptr<ProperOrthogonalDecomposition::ROMTestLocation<dim,nstate>>>> rom_locations;

    /// Maximum error
    mutable double max_error;

    /// Parameter handler for storing the .prm file being ran
    const dealii::ParameterHandler &parameter_handler;

    /// Most up to date POD basis
    std::shared_ptr<ProperOrthogonalDecomposition::OnlinePOD<dim>> current_pod;

    /// Nearest neighbors of snapshots
    std::shared_ptr<ProperOrthogonalDecomposition::NearestNeighbors> nearest_neighbors;

    /// Adaptation tolerance
    double tolerance;

    /// Run test
    int run_test () const override;

    /// Placement of initial snapshots
    void placeInitialSnapshots() const;

    /// Placement of ROMs
    bool placeTriangulationROMs(const MatrixXd& rom_points) const;

    /// Updates nearest ROM points to snapshot if error discrepancy is above tolerance
    void updateNearestExistingROMs(const RowVectorXd& parameter) const;

    /// Compute RBF and find max error
    RowVectorXd getMaxErrorROM() const;

    /// Solve full-order snapshot
    dealii::LinearAlgebra::distributed::Vector<double> solveSnapshotFOM(const RowVectorXd& parameter) const;

    /// Solve reduced-order solution
    std::shared_ptr<ProperOrthogonalDecomposition::ROMSolution<dim,nstate>> solveSnapshotROM(const RowVectorXd& parameter) const;

    /// Reinitialize parameters
    Parameters::AllParameters reinitParams(const RowVectorXd& parameter) const;

    /// Set up parameter space depending on test case
    void configureParameterSpace() const;

    /// Output for each iteration
    void outputErrors(int iteration) const;
};

}
}


#endif
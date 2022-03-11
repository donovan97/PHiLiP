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
#include "snapshot.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

/// Class to hold information about the reduced-order solution
template <int dim, int nstate>
class AdaptiveSampling
{
public:
    /// Constructor
    AdaptiveSampling(std::vector<double> parameter_space, double tolerance);

    /// Destructor
    virtual ~AdaptiveSampling() {};

    std::vector<double> parameter_space;

    double tolerance;

    std::vector<Snapshot<dim,nstate>> snapshots;

    std::shared_ptr<POD<dim>> current_pod;

    void start();

    void generateTrialLocations();

    void initializeSampling();

    void updateErrors();

    void updateSensitivityCurveFit();

    void updateErrorCurveFit();

    void addSnapshot();

    void solveSnapshotFOM();

    void solveSnapshotROM();

};

}
}


#endif

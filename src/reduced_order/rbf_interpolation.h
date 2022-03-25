#ifndef __RBF_INTERPOLATION__
#define __RBF_INTERPOLATION__

#include <Eigen/Dense>
#include <Eigen/LU>
#include <iostream>

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

/// Radial basis function interpolation
class RBFInterpolation
{
public:
    /// Constructor
    RBFInterpolation(MatrixXd data_coordinates, VectorXd data_values, std::string kernel);

    /// Destructor
    ~RBFInterpolation () {};

    void computeWeights();

    double radialBasisFunction(double r);

    VectorXd evaluate(RowVectorXd evaluate_coordinate);

    VectorXd weights;

    const MatrixXd data_coordinates;

    const VectorXd data_values;

    const std::string kernel;

};

}
}


#endif

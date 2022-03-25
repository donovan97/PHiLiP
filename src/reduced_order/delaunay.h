#ifndef __DELAUNAY__
#define __DELAUNAY__

#include <Eigen/Dense>
#include <Eigen/LU>
#include <iostream>

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

/// Radial basis function interpolation
class Delaunay
{
public:
    /// Constructor
    Delaunay();

    /// Destructor
    ~Delaunay() {};

};

}
}

#endif

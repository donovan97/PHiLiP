#ifndef __LEVENBERG_MARQUARDT__
#define __LEVENBERG_MARQUARDT__

#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include "rbf_interpolation.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

struct MyFunctor{
    typedef double Scalar;

    typedef Eigen::VectorXd InputType;
    typedef Eigen::VectorXd ValueType;
    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> JacobianType;

    enum {
        InputsAtCompileTime = Eigen::Dynamic,
        ValuesAtCompileTime = Eigen::Dynamic
    };

    RBFInterpolation rbf;

    MyFunctor(RBFInterpolation rbf);

    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const;
    int inputs() const;
    int values() const;
};

}
}

#endif

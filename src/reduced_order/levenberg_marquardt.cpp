#include "levenberg_marquardt.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

MyFunctor::MyFunctor(RBFInterpolation rbf)
    : rbf(rbf)
{}

int MyFunctor::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
{
    // IMPORTANT: The algorithm will square fvec internally. Therefore, for this rbf, taking the inverse
    // of the absolute value will give the right minimum
    RowVectorXd rbf_x = rbf.evaluate(x.transpose());
    fvec = rbf_x.transpose().cwiseAbs().cwiseInverse();
    return 0;
}

int MyFunctor::inputs() const { return 1; }// inputs is the dimension of x.
int MyFunctor::values() const { return 1; } // "values" is the number of f_i and


}
}
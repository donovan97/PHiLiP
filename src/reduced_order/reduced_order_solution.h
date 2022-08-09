#ifndef __REDUCED_ORDER_SOLUTION__
#define __REDUCED_ORDER_SOLUTION__

#include "functional/functional.h"
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector_operation.h>
#include "parameters/all_parameters.h"
#include "dg/dg.h"
#include "pod_basis_base.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

/// Class to hold information about the reduced-order solution
template<int dim, int nstate>
class ROMSolution
{
public:
    /// Constructor
    ROMSolution(dealii::TrilinosWrappers::SparseMatrix _system_matrix_transpose, dealii::LinearAlgebra::distributed::Vector<double> _right_hand_side, dealii::TrilinosWrappers::SparseMatrix _pod_basis, dealii::LinearAlgebra::distributed::Vector<double> _gradient);

    /// Destructor
    ~ROMSolution () {};

    /// Stores system matrix transpose
    dealii::TrilinosWrappers::SparseMatrix system_matrix_transpose;

    /// Stores residual
    dealii::LinearAlgebra::distributed::Vector<double> right_hand_side;

    /// Stores POD basis on which solution was computed
    dealii::TrilinosWrappers::SparseMatrix basis;

    /// Stores gradient
    dealii::LinearAlgebra::distributed::Vector<double> gradient;

};

}
}


#endif
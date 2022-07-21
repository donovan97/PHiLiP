#include "reduced_order_solution.h"

#include <utility>

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template <int dim, int nstate>
ROMSolution<dim, nstate>::ROMSolution(std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> _system_matrix_transpose, dealii::LinearAlgebra::distributed::Vector<double> _right_hand_side, std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> _pod_basis, dealii::LinearAlgebra::distributed::Vector<double> _gradient)
        : system_matrix_transpose(std::move(_system_matrix_transpose))
        , right_hand_side(_right_hand_side)
        , basis(std::move(_pod_basis))
        , gradient(_gradient)
{
}

template class ROMSolution <PHILIP_DIM, 1>;
template class ROMSolution <PHILIP_DIM, 2>;
template class ROMSolution <PHILIP_DIM, 3>;
template class ROMSolution <PHILIP_DIM, 4>;
template class ROMSolution <PHILIP_DIM, 5>;

}
}
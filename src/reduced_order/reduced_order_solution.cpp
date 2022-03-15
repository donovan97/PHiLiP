#include "reduced_order_solution.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template <int dim, int nstate>
ROMSolution<dim, nstate>::ROMSolution(std::shared_ptr<DGBase<dim,double>> &dg_input, Functional<dim,nstate,double> &functional_input, std::shared_ptr<ProperOrthogonalDecomposition::POD<dim>> pod_input)
        : system_matrix_transpose(dealii::TrilinosWrappers::SparseMatrix())
        , right_hand_side(dg_input->right_hand_side)
        , basis(pod_input->getPODBasis().get())
        , functional_value(functional_input.evaluate_functional( true, false, false))
        , gradient(functional_input.dIdw)
{
    system_matrix_transpose.copy_from(dg_input->system_matrix_transpose);
}

template class ROMSolution <PHILIP_DIM, 1>;
template class ROMSolution <PHILIP_DIM, 2>;
template class ROMSolution <PHILIP_DIM, 3>;
template class ROMSolution <PHILIP_DIM, 4>;
template class ROMSolution <PHILIP_DIM, 5>;

}
}

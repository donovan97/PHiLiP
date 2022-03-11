#include "reduced_order_solution.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template <int dim, int nstate>
ROMSolution<dim, nstate>::ROMSolution(std::shared_ptr<DGBase<dim,double>> &dg_input, Functional<dim,nstate,double> &functional_input, std::shared_ptr<ProperOrthogonalDecomposition::POD<dim>> pod_input)
        : system_matrix_transpose(std::make_shared<dealii::TrilinosWrappers::SparseMatrix>())
        , right_hand_side(dealii::LinearAlgebra::distributed::Vector<double>())
        , basis(std::make_shared<dealii::TrilinosWrappers::SparseMatrix>())
        , gradient(dealii::LinearAlgebra::distributed::Vector<double>())
{
    system_matrix_transpose->copy_from(dg_input->system_matrix_transpose);
    basis->copy_from(*pod_input->getPODBasis());
    right_hand_side.reinit(dg_input->right_hand_side);

    const bool compute_dIdW = true;
    const bool compute_dIdX = false;
    const bool compute_d2I = false;
    functional_value = functional_input.evaluate_functional( compute_dIdW, compute_dIdX, compute_d2I);
    gradient.reinit(functional_input.dIdw);
}

template class ROMSolution <PHILIP_DIM, 1>;
template class ROMSolution <PHILIP_DIM, 2>;
template class ROMSolution <PHILIP_DIM, 3>;
template class ROMSolution <PHILIP_DIM, 4>;
template class ROMSolution <PHILIP_DIM, 5>;

}
}

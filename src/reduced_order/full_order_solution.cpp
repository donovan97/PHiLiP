#include "full_order_solution.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template <int dim, int nstate>
FOMSolution<dim, nstate>::FOMSolution(std::shared_ptr<DGBase<dim,double>> &dg_input, Functional<dim,nstate,double> &functional_input, dealii::LinearAlgebra::distributed::Vector<double> sensitivity)
        : sensitivity(sensitivity)
{
    state = dg_input->solution;

    const bool compute_dIdW = false;
    const bool compute_dIdX = false;
    const bool compute_d2I = false;
    functional_value = functional_input.evaluate_functional(compute_dIdW, compute_dIdX, compute_d2I);
}

template class FOMSolution <PHILIP_DIM, 1>;
template class FOMSolution <PHILIP_DIM, 2>;
template class FOMSolution <PHILIP_DIM, 3>;
template class FOMSolution <PHILIP_DIM, 4>;
template class FOMSolution <PHILIP_DIM, 5>;

}
}

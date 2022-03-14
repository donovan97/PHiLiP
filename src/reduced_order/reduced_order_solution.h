#ifndef __REDUCED_ORDER_SOLUTION__
#define __REDUCED_ORDER_SOLUTION__

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

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

/// Class to hold information about the reduced-order solution
template<int dim, int nstate>
class ROMSolution
{
public:
    /// Constructor
    ROMSolution(std::shared_ptr<DGBase<dim,double>> &dg_input, Functional<dim,nstate,double> &functional_input, std::shared_ptr<ProperOrthogonalDecomposition::POD<dim>>);

    /// Destructor
    virtual ~ROMSolution () {};

    const std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> system_matrix_transpose;

    const dealii::LinearAlgebra::distributed::Vector<double> right_hand_side;

    const std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> basis;

    const double functional_value;

    const dealii::LinearAlgebra::distributed::Vector<double> gradient;

};

}
}


#endif

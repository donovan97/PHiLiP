#ifndef __POD_BASIS_SENSITIVITY_TYPES__
#define __POD_BASIS_SENSITIVITY_TYPES__

#include <fstream>
#include <iostream>
#include <filesystem>

#include "functional/functional.h"
#include "dg/dg.h"
#include "reduced_order/pod_basis.h"
#include "reduced_order/pod_basis_sensitivity.h"
#include "linear_solver/linear_solver.h"

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include "ode_solver/ode_solver_factory.h"
#include <algorithm>
#include "pod_basis_interface.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template<int dim>
/// Intermediary class that includes attributes common to all SensitivityPOD basis subtypes
class SpecificSensitivityPOD : public SensitivityPOD<dim> {
protected:
    /// Constructor
    SpecificSensitivityPOD(std::shared_ptr<DGBase<dim, double>> &dg_input);

    /// Destructor
    ~SpecificSensitivityPOD() {}

    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> basis; ///< pod basis
    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> basisTranspose; ///< Transpose of pod_basis

public:

    /// Function to add columns (basis functions) to POD basis. Used when building basis and refining when doing POD adaptation
    virtual void addPODBasisColumns(const std::vector<unsigned int> addColumns) = 0;

    /// Function to return basis
    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> getPODBasis() override;

    /// Function to return basisTranspose
    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> getPODBasisTranspose() override;

    /// Vector to store which indicies of the full basis are present in this basis
    std::vector<unsigned int> fullBasisIndices;

    /// Vector to store which indicies of the sensitivity basis are present in this basis
    std::vector<unsigned int> fullSensitivityBasisIndices;

};

/// Class for a coarse expanded POD basis, derived from SensitivityPOD and implementing CoarseBasis
template<int dim>
class CoarseExpandedPOD : public SensitivityPOD<dim>, public CoarseBasis<dim> {
private:

    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> coarseBasis; ///< pod basis
    std::vector<unsigned int> fullBasisIndices; ///< Vector to store which indicies of the full basis are present in this basis

public:
    /// Constructor
    CoarseExpandedPOD(std::shared_ptr<DGBase<dim, double>> &dg_input);

    /// Destructor
    ~CoarseExpandedPOD() {};

    /// Add columns to the basis
    void addPODBasisColumns(const std::vector<unsigned int> addColumns);

    /// Function to return basis
    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> getPODBasis();

    /// Return vector storing which indices of the full basis are present in this basis
    std::vector<unsigned int> getFullBasisIndices();
};

/*
/// Class for a fine not in coarse expanded POD basis, derived from SensitivityPOD and implementing FineNotInCoarseBasis
template<int dim>
class FineNotInCoarseExpandedPOD : public SensitivityPOD<dim>, public FineNotInCoarseBasis<dim> {
private:

    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> basis; ///< pod basis
    std::vector<unsigned int> fullBasisIndices; ///< Vector to store which indicies of the full basis are present in this basis

public:
    /// Constructor
    FineNotInCoarseExpandedPOD(std::shared_ptr<DGBase<dim, double>> &dg_input);

    /// Destructor
    ~FineNotInCoarseExpandedPOD() {};

    /// Removes columns of the basis, used during POD adaptation
    void removePODBasisColumns(const std::vector<unsigned int> removeColumns);

    /// Add columns to the basis
    void addPODBasisColumns(const std::vector<unsigned int> addColumns);

    /// Function to return basis
    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> getPODBasis();

    /// Return vector storing which indices of the full basis are present in this basis
    std::vector<unsigned int> getFullBasisIndices();
};

/// Class for a fine not in coarse expanded POD basis, derived from SensitivityPOD and implementing FineNotInCoarseBasis
template<int dim>
class FineExpandedPOD : public SensitivityPOD<dim>, public FineBasis<dim> {
private:

    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> basis; ///< pod basis
    std::vector<unsigned int> fullBasisIndices; ///< Vector to store which indicies of the full basis are present in this basis

public:
    /// Constructor
    FineExpandedPOD(std::shared_ptr<DGBase<dim, double>> &dg_input);

    /// Destructor
    ~FineExpandedPOD() {};

    /// Add columns to the basis
    void addPODBasisColumns(const std::vector<unsigned int> addColumns);

    /// Function to return basis
    std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> getPODBasis();

    /// Return vector storing which indices of the full basis are present in this basis
    std::vector<unsigned int> getFullBasisIndices();
};
*/
/// Class for Extrapolated POD basis, derived from SpecificSensitivityPOD
template<int dim>
class ExtrapolatedPOD : public SpecificSensitivityPOD<dim> {
public:
    /// Constructor
    ExtrapolatedPOD(std::shared_ptr<DGBase<dim, double>> &dg_input);

    /// Destructor
    ~ExtrapolatedPOD() {};

    void addPODBasisColumns(const std::vector<unsigned int> addColumns);
};

}
}

#endif

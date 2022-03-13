#include "pod_basis_online.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template<int dim>
OnlinePOD<dim>::OnlinePOD()
        : basis(std::make_shared<dealii::TrilinosWrappers::SparseMatrix>())
{}

template <int dim>
void OnlinePOD<dim>::addSnapshot(dealii::LinearAlgebra::distributed::Vector<double> snapshot) {
    snapshotVectors.push_back(snapshot);
}

template <int dim>
void OnlinePOD<dim>::computeBasis() {
    dealii::LAPACKFullMatrix<double> snapshot_matrix(snapshotVectors[0].size(), snapshotVectors.size());

    for(unsigned int n = 0 ; n < snapshotVectors.size() ; n++){
        for(unsigned int m = 0 ; m < snapshotVectors[0].size() ; m++){
            snapshot_matrix(m,n) = snapshotVectors[n][m];
        }
    }

    snapshot_matrix.compute_svd();
    dealii::LAPACKFullMatrix<double> svd_u = snapshot_matrix.get_svd_u();

    dealii::TrilinosWrappers::SparseMatrix basis_tmp(snapshotVectors[0].size(), snapshotVectors.size(), snapshotVectors.size());
    std::cout << snapshotVectors[0].size() << " " << snapshotVectors.size() << std::endl;
    for(unsigned int m = 0 ; m < snapshotVectors[0].size() ; m++){
        for(unsigned int n = 0 ; n < snapshotVectors.size() ; n++){
            basis_tmp.set(m, n, svd_u(m,n));
        }
    }

    basis_tmp.compress(dealii::VectorOperation::insert);

    basis->reinit(basis_tmp);
    basis->copy_from(basis_tmp);
}

template <int dim>
std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> OnlinePOD<dim>::getPODBasis() {
    return basis;
}

template class OnlinePOD <PHILIP_DIM>;

}
}


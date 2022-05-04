#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include "pod_basis_online.h"
#include <deal.II/base/index_set.h>

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template<int dim>
OnlinePOD<dim>::OnlinePOD(std::shared_ptr<DGBase<dim,double>> &dg_input)
        : basis(std::make_shared<dealii::TrilinosWrappers::SparseMatrix>())
        , dg(dg_input)
        , snapshotMatrix(0,0)
{}

template <int dim>
void OnlinePOD<dim>::addSnapshot(dealii::LinearAlgebra::distributed::Vector<double> snapshot) {
    std::cout << "Adding new snapshot to snapshot matrix..." << std::endl;
    dealii::LinearAlgebra::ReadWriteVector<double> read_snapshot(snapshot.size());
    read_snapshot.import(snapshot, dealii::VectorOperation::values::insert);
    VectorXd eigen_snapshot(snapshot.size());
    for(unsigned int i = 0 ; i < snapshot.size() ; i++){
        eigen_snapshot(i) = read_snapshot(i);
    }
    snapshotMatrix.conservativeResize(snapshot.size(), snapshotMatrix.cols()+1);
    snapshotMatrix.col(snapshotMatrix.cols()-1) = eigen_snapshot;
}

template <int dim>
void OnlinePOD<dim>::computeBasis() {
    std::cout << "Computing POD basis..." << std::endl;

    reference_state = snapshotMatrix.rowwise().mean();

    //for(unsigned int i = 0 ; i < reference_state.size() ; i++){
    //    reference_state(i) = 1;
    //}

    snapshotMatrix = snapshotMatrix.colwise() - reference_state;

    Eigen::BDCSVD<MatrixXd> svd(snapshotMatrix, Eigen::DecompositionOptions::ComputeThinU);
    pod_basis = svd.matrixU();

    dealii::TrilinosWrappers::SparseMatrix basis_tmp(snapshotMatrix.rows(), snapshotMatrix.cols(), snapshotMatrix.rows());
    for(unsigned int m = 0 ; m < snapshotMatrix.rows() ; m++){
        for(unsigned int n = 0 ; n < snapshotMatrix.cols() ; n++){
            //std::cout << pod_basis(m,n) << std::endl;
            basis_tmp.set(m, n, pod_basis(m,n));
        }
    }

    fullBasis.reinit(snapshotMatrix.rows(), snapshotMatrix.cols());
    for(unsigned int m = 0 ; m < snapshotMatrix.rows() ; m++){
        for(unsigned int n = 0 ; n < snapshotMatrix.cols() ; n++){
            fullBasis.set(m, n, pod_basis(m,n));
        }
    }

    basis_tmp.compress(dealii::VectorOperation::insert);

    basis->reinit(basis_tmp);
    basis->copy_from(basis_tmp);

    std::cout << "Done computing POD basis. Basis now has " << pod_basis.cols() << " columns." << std::endl;
}

template <int dim>
std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> OnlinePOD<dim>::getPODBasis() {
    return basis;
}

template <int dim>
MatrixXd OnlinePOD<dim>::getEigenPODBasis() {
    return pod_basis;
}

template <int dim>
VectorXd OnlinePOD<dim>::getReferenceState() {
    return reference_state;
}

template class OnlinePOD <PHILIP_DIM>;

}
}

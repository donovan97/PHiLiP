#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include "pod_basis_online.h"
#include <deal.II/base/index_set.h>

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template<int dim>
OnlinePOD<dim>::OnlinePOD(std::shared_ptr<DGBase<dim,double>> &dg_input)
        : basis(std::make_shared<dealii::TrilinosWrappers::SparseMatrix>())
        , mass_matrix_sparsity(dg_input->global_mass_matrix.trilinos_sparsity_pattern())
        , dg(dg_input)
        , snapshotMatrix(0,0)
{
    const bool compute_dRdW = true;
    dg_input->assemble_residual(compute_dRdW);
    massMatrix.reinit(dg_input->global_mass_matrix.m(), dg_input->global_mass_matrix.n());
    massMatrix.copy_from(dg_input->global_mass_matrix);
}

template <int dim>
void OnlinePOD<dim>::addSnapshot(dealii::LinearAlgebra::distributed::Vector<double> snapshot) {
    std::cout << "Adding new snapshot to snapshot matrix..." << std::endl;
    dealii::LinearAlgebra::ReadWriteVector<double> read_snapshot(snapshot.size());
    read_snapshot.import(snapshot, dealii::VectorOperation::values::insert);
    //snapshotVectors.push_back(read_snapshot);
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

template class OnlinePOD <PHILIP_DIM>;

}
}

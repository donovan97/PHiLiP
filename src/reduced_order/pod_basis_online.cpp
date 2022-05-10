#include "pod_basis_online.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template<int dim>
OnlinePOD<dim>::OnlinePOD(std::shared_ptr<DGBase<dim,double>> &dg_input)
        : basis(std::make_shared<dealii::TrilinosWrappers::SparseMatrix>())
        , dg(dg_input)
        , mpi_communicator(MPI_COMM_WORLD)
        , mpi_rank(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
        , pcout(std::cout, mpi_rank==0)
{
    const bool compute_dRdW = true;
    dg_input->assemble_residual(compute_dRdW);
}

template <int dim>
void OnlinePOD<dim>::addSnapshot(dealii::LinearAlgebra::distributed::Vector<double> snapshot) {
    pcout << "Adding new snapshot to snapshot matrix..." << std::endl;
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
    pcout << "Computing POD basis..." << std::endl;

    VectorXd reference_state = snapshotMatrix.rowwise().mean();

    //for(unsigned int i = 0 ; i < reference_state.size() ; i++){
    //    reference_state(i) = 1;
    //}

    referenceState.reinit(reference_state.size());
    for(unsigned int i = 0 ; i < reference_state.size() ; i++){
        referenceState(i) = reference_state(i);
    }

    snapshotMatrix = snapshotMatrix.colwise() - reference_state;

    Eigen::BDCSVD<MatrixXd> svd(snapshotMatrix, Eigen::DecompositionOptions::ComputeThinU);
    MatrixXd pod_basis = svd.matrixU();

    fullBasis.reinit(pod_basis.rows(), pod_basis.cols());

    for (unsigned int m = 0; m < pod_basis.rows(); m++) {
        for (unsigned int n = 0; n < pod_basis.cols(); n++) {
            fullBasis.set(m, n, pod_basis(m, n));
        }
    }

    std::ofstream out_file("POD_adaptation_basis.txt");
    unsigned int precision = 16;
    fullBasis.print_formatted(out_file, precision);


    Epetra_CrsMatrix *epetra_system_matrix  = const_cast<Epetra_CrsMatrix *>(&(this->dg->system_matrix.trilinos_matrix()));
    Epetra_Map system_matrix_map = epetra_system_matrix->RowMap();
    Epetra_CrsMatrix epetra_basis(Epetra_DataAccess::Copy, system_matrix_map, pod_basis.cols());

    const int numMyElements = system_matrix_map.NumMyElements(); //Number of elements on the calling processor

    for (int localRow = 0; localRow < numMyElements; ++localRow){
        const int globalRow = system_matrix_map.GID(localRow);
        //pcout << "global row: " << globalRow << std::endl;
        for(int n = 0 ; n < pod_basis.cols() ; n++){
            epetra_basis.InsertGlobalValues(globalRow, 1, &pod_basis(globalRow, n), &n);
        }
    }

    Epetra_MpiComm epetra_comm(MPI_COMM_WORLD);
    Epetra_Map domain_map((int)pod_basis.cols(), 0, epetra_comm);

    epetra_basis.FillComplete(domain_map, system_matrix_map);
    //epetra_basis.FillComplete();

    /*
    for (int localRow = 0; localRow < numMyElements; ++localRow){
        const int globalRow = system_matrix_map.GID(localRow);
        int numentries;
        double* values;
        epetra_basis.ExtractGlobalRowView(globalRow, numentries, values);

        for(int n = 0 ; n < pod_basis.cols() ; n++){
            pcout << "global row: " << globalRow << "entry: " << values[n] << std::endl;
        }
    }
    */
    basis->reinit(epetra_basis);

    pcout << "Done computing POD basis. Basis now has " << basis->n() << " columns." << std::endl;
}

template <int dim>
std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> OnlinePOD<dim>::getPODBasis() {
    return basis;
}

template <int dim>
dealii::LinearAlgebra::ReadWriteVector<double> OnlinePOD<dim>::getReferenceState() {
    return referenceState;
}

template class OnlinePOD <PHILIP_DIM>;

}
}
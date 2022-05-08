#include "pod_basis_online.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

template<int dim>
OnlinePOD<dim>::OnlinePOD(std::shared_ptr<DGBase<dim,double>> &dg_input)
        : basis(std::make_shared<dealii::TrilinosWrappers::SparseMatrix>())
        , dg(dg_input)
{
    const bool compute_dRdW = true;
    dg_input->assemble_residual(compute_dRdW);
}

template <int dim>
void OnlinePOD<dim>::addSnapshot(dealii::LinearAlgebra::distributed::Vector<double> snapshot) {
    std::cout << "Adding new snapshot to POD basis..." << std::endl;
    dealii::LinearAlgebra::ReadWriteVector<double> snapshot_vector(snapshot.size());
    snapshot_vector.import(snapshot, dealii::VectorOperation::values::insert);
    snapshotVectors.push_back(snapshot_vector);
}

template <int dim>
void OnlinePOD<dim>::computeBasis() {
    std::cout << "Computing POD basis..." << std::endl;

    Teuchos::RCP<const Teuchos::Comm<int> > comm (new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
    const int myRank = comm->getRank ();

    std::cout << "rank: " << myRank << std::endl;

    dealii::LAPACKFullMatrix<double> snapshot_matrix(snapshotVectors[0].size(), snapshotVectors.size());

    for (unsigned int n = 0; n < snapshotVectors.size(); n++) {
        for (unsigned int m = 0; m < snapshotVectors[0].size(); m++) {
            snapshot_matrix(m, n) = snapshotVectors[n][m];
        }
    }

    std::ofstream out_file_snap("snapshot_matrix.txt");
    unsigned int precision0 = 16;
    snapshot_matrix.print_formatted(out_file_snap, precision0);


    snapshot_matrix.compute_svd();
    dealii::LAPACKFullMatrix<double> svd_u = snapshot_matrix.get_svd_u();


    fullBasis.reinit(snapshot_matrix.m(), snapshot_matrix.n());

    for (unsigned int m = 0; m < snapshotVectors[0].size(); m++) {
        for (unsigned int n = 0; n < snapshotVectors.size(); n++) {
            fullBasis.set(m, n, svd_u(m, n));
        }
    }

    std::ofstream out_file("POD_adaptation_basis.txt");
    unsigned int precision = 16;
    fullBasis.print_formatted(out_file, precision);


    Epetra_CrsMatrix *epetra_system_matrix  = const_cast<Epetra_CrsMatrix *>(&(this->dg->system_matrix.trilinos_matrix()));
    Epetra_Map system_matrix_map = epetra_system_matrix->RowMap();
    Epetra_CrsMatrix epetra_basis(Epetra_DataAccess::Copy, system_matrix_map, snapshotVectors.size());

    const int numMyElements = system_matrix_map.NumMyElements(); //Number of elements on the calling processor

    for (int localRow = 0; localRow < numMyElements; ++localRow){
        const int globalRow = system_matrix_map.GID(localRow);
        std::cout << "global row: " << globalRow << std::endl;
        for(int n = 0 ; n < (int)snapshotVectors.size() ; n++){
            epetra_basis.InsertGlobalValues(globalRow, 1, &fullBasis(globalRow, n), &n);
        }
    }

    Epetra_MpiComm epetra_comm(MPI_COMM_WORLD);
    Epetra_Map domain_map((int)snapshotVectors.size(), 0, epetra_comm);

    epetra_basis.FillComplete(domain_map, system_matrix_map);
    //epetra_basis.FillComplete();

    for (int localRow = 0; localRow < numMyElements; ++localRow){
        const int globalRow = system_matrix_map.GID(localRow);
        int numentries;
        double* values;
        epetra_basis.ExtractGlobalRowView(globalRow, numentries, values);

        for(int n = 0 ; n < (int)snapshotVectors.size() ; n++){
            std::cout << "global row: " << globalRow << "entry: " << values[n] << std::endl;
        }
    }

    basis->reinit(epetra_basis);

    //Epetra_CrsMatrix output(Epetra_DataAccess::Copy, system_matrix_map, snapshotVectors.size());

    //EpetraExt::MatrixMatrix::Multiply(*epetra_system_matrix, false, epetra_basis, false, output);

    //dealii::TrilinosWrappers::SparseMatrix import_output;
    //import_output.reinit(output);

    std::cout << "Done computing POD basis. Basis now has " << basis->n() << " columns." << std::endl;
}

template <int dim>
std::shared_ptr<dealii::TrilinosWrappers::SparseMatrix> OnlinePOD<dim>::getPODBasis() {
    return basis;
}

template class OnlinePOD <PHILIP_DIM>;

}
}
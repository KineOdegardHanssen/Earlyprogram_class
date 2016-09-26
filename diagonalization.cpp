#include "diagonalization.h"

Diagonalization::Diagonalization()
{
}

Diagonalization::Diagonalization(Systems given)
{
    this->given = given;
}

void Diagonalization::using_armadillo()
{
    arma_n = 0;
    if(given.sectorbool==true)          arma_n = given.number_of_hits;
    else                                arma_n = given.no_of_states;
    eigenvalues_armadillo = arma::vec(arma_n);
    eigenvectors_armadillo = arma::mat(arma_n,arma_n);
    //double start_arma = clock();        // If I want to take the time
    //string method = "std";
    arma::eig_sym(eigenvalues_armadillo,eigenvectors_armadillo, given.armaH);
    //double end_arma = clock();
}

void Diagonalization::print_using_armadillo()
{
    for(unsigned long i= 0; i<arma_n; i++)
    {
        cout << "Eigenvalue = " << eigenvalues_armadillo(i) << endl;
        cout << "Its corresponding eigenvector:" << endl;
        for(unsigned long j=0; j<arma_n; j++)       cout << eigenvectors_armadillo(j,i) << " ";
        cout << endl;
    }   // End for-loop over i
}

void Diagonalization::print_sparse_Eigen()
{
    // See if i can find something a bit more clever...
    cout << given.sparseH << endl;
}


/*
void Diagonalization::using_dense_eigen()
{
    arma_n = 0;
    if(given.sectorbool==true)          arma_n = given.number_of_hits;
    else                                arma_n = given.no_of_states;

    Eigen::EigenSolver<Eigen::MatrixXd> es(given.eigenH.cast<double>);

    Eigen::VectorXd eigenvalues_H = es.eigenvalues();
    Eigen::MatrixXd eigenmatrix_H = es.eigenvectors();

}


void Diagonalization::print_dense_using_eigen()
{
    Eigen::EigenSolver<Eigen::MatrixXd> es(given.eigenH.cast<double>());
    cout << "The eigenvalues of H are:" << endl << es.eigenvalues() << endl;
    cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
}
*/




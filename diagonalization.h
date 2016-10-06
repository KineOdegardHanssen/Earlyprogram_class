#ifndef DIAGONALIZATION_H
#define DIAGONALIZATION_H
#include <iostream>
#include <systems.h>
#include <armadillo>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include<vector>
#include<complex>


using namespace std;
//using namespace arma;

class Diagonalization
{
public:

    Systems given;

    // For linking directly to LAPACK
    char JOBS;
    char UPLO;  // upper triangular matrix is stored after diag.
    int N;      // size of matrix
    int LDA;


    // For armadillo solvers
    unsigned long arma_n, eigen_n;
    double armatime, lapacktime;
    arma::vec eigenvalues_armadillo; // Should I use arma::vec? Probably a good idea. Eigen probably has something similar
    arma::mat eigenvectors_armadillo;

    // For Eigen using dense matrices
    Eigen::VectorXd eigenvalues_H;      // Not using these (yet...)
    Eigen::MatrixXd eigenmatrix_H;


    //Initializers
    Diagonalization();
    Diagonalization(Systems given);



    // Functions for using LAPACK
    void lapack_directly();

    // Functions for armadillo matrices (dense)
    void using_armadillo();
    void print_using_armadillo();

    void using_dense_eigen();
    void print_dense_using_eigen();

    // Functions for Eigen::SparseMatrix
    void print_sparse_Eigen();
};

#endif // DIAGONALIZATION_H

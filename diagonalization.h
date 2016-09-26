#ifndef DIAGONALIZATION_H
#define DIAGONALIZATION_H
#include <iostream>
#include <systems.h>
#include <armadillo>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>

using namespace std;
//using namespace arma;

class Diagonalization
{
public:

    Systems given;

    // For armadillo solvers
    unsigned long arma_n, eigen_n;
    arma::vec eigenvalues_armadillo; // Should I use arma::vec? Probably a good idea. Eigen probably has something similar
    arma::mat eigenvectors_armadillo;

    // For Eigen using dense matrices
    Eigen::VectorXd eigenvalues_H;
    Eigen::MatrixXd eigenmatrix_H;


    //Initializers
    Diagonalization();
    Diagonalization(Systems given);


    // Functions for armadillo matrices (dense)
    void using_armadillo();
    void print_using_armadillo();

    void using_dense_eigen();
    void print_dense_using_eigen();

    // Functions for Eigen::SparseMatrix
    void print_sparse_Eigen();
};

#endif // DIAGONALIZATION_H

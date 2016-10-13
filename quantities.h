#ifndef QUANTITIES_H
#define QUANTITIES_H
#include <iostream>
#include <systems.h>
#include <diagonalization.h>
#include <armadillo>
#include <Eigen/Dense>
//#include <Eigen/Sparse>
#include <Eigen/QR>            // I don't really use this...
#include <Eigen/Eigenvalues>   // I guess I don't need this here..
#include<vector>
#include<complex>


class Quantities
{
public:

    // One class on top of this, like runnit for 2p?


    Systems system;            // Should I do something like this and nest them?  // Worry about this later
    Diagonalization eig_all;

    int N, maxit, systemsize, li, lh;
    double Z, beta, min_ev, tolerance; // Only change smallest_ev this for every new instance of quantities...
    bool armadillobool;

    Eigen::VectorXd eigvals;
    Eigen::MatrixXd eigmat;

    arma::vec eigvals_a;  // This will probably be a problem...
    arma::mat eigmat_a;   // Could rename it at each step. Troublesome...


    // Initializer
    Quantities(int maxit, double tolerance, bool armadillobool, Systems system, Diagonalization eig_all);


    //Functions

    // Basic functions
    void sort_energies();
    arma::mat initialize_matrix_arma(int size);
    Eigen::MatrixXd initialize_matrix_Eigen(int size);
    int signcompare(double fa, double fc);
    void calculateZ();
    void calculateZ_arma();

    // Finding beta
    // To be run for each eigenstate/eigenvalue
    // I may be giving up some speed or accuracy here: can divide by e^(-beta E_lowest), but not set e^(-beta(E_lowest-E_lowest)) to 1 automatically
    void newtonsmethod(double eigenvalue);
    void bisectionmethod(double eigenvalue);
    double self_consistency_beta(double eigenvalue, double betatest);
    double self_consistency_beta_derivative(double eigenvalue, double betatest);
    // armadillo
    void newtonsmethod_arma(double eigenvalue);
    void bisectionmethod_arma(double eigenvalue);
    double self_consistency_beta_a(double eigenvalue, double betatest);
    double self_consistency_beta_derivative_a(double eigenvalue, double betatest);

    // Eigenstate Thermalization Hypothesis
    double ETH(int i);       // Should this really be a void?
    double ETH_arma(int i);
    double ETH_arma_sector(int i);
    double ETH_arma_maziero(int i);
    double ETH_Eigen(int i);
    double ETH_Eigen_sector(int i);
    double ETH_Eigen_maziero(int i);

    // Eigen
    Eigen::MatrixXd trace_Eigen(Eigen::MatrixXd A);
    Eigen::MatrixXd trace_Eigen_maziero(Eigen::MatrixXd A);
    Eigen::MatrixXd trace_Eigen_sector(Eigen::MatrixXd A);
    Eigen::MatrixXd thermalmat_Eigen();             // See if I change this a bit.
    Eigen::MatrixXd eigenstatemat_Eigen(int i);

    // Armadillo
    arma::mat trace_arma(arma::mat A);
    arma::mat trace_arma_maziero(arma::mat A);
    arma::mat trace_arma_sector(arma::mat A);
    arma::mat thermalmat_arma();                    // See if I change this a bit.
    arma::mat eigenstatemat_arma(int i);
};

#endif // QUANTITIES_H

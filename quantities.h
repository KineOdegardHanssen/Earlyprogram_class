#ifndef QUANTITIES_H
#define QUANTITIES_H
#include <iostream>
#include <systems.h>
#include <diagonalization.h>
#include <armadillo>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include<vector>
#include<complex>


class Quantities
{
public:

    // One class on top of this, like runnit for 2p?


    //Systems system;            // Should I do something like this and nest them?  // Worry about this later
    Diagonalization eig_all;

    int N;
    double Z, beta, min_ev; // Only change smallest_ev this for every new instance of quantities...

    Eigen::VectorXd eigvals;
    Eigen::MatrixXd eigmat;

    // Initializer
    Quantities(Systems system);


    //Functions

    // Basic functions
    void min_eigenvalue_of_sector();

    // Finding beta
    // To be run for each eigenstate/eigenvalue
    // I may be giving up some speed or accuracy here: can divide by e^(-beta E_lowest), but not set e^(-beta(E_lowest-E_lowest)) to 1 automatically
    void findbeta(double eigenvalue);
    double self_consistency_beta(double eigenvalue, double betatest);
    double self_consistency_beta_derivative(double eigenvalue, double betatest);
};

#endif // QUANTITIES_H

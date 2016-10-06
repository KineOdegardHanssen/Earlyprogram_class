#ifndef SYSTEMS_H
#define SYSTEMS_H
#include <iostream>
#include <vector>
#include <cmath>
#include <armadillo>
#include <random>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/QR> 
#include <Eigen/Eigenvalues>

using namespace std;

class Systems
{
public:

    int systemsize, number_of_hits, no_of_states, mysector;
    double h, J;

    bool dense, palhuse, sectorbool, ctbool ,armadillobool, testupip1downi, testdownip1upi;

    vector<double> hs;
    vector<int> sectorlist;
    typedef Eigen::Triplet<double> T;

    // For armadillo solvers        // Okay, the program flow certainly is awkward...
    arma::mat armaH;
    arma::cx_mat armaHcompl;

    // For Eigen solvers
    Eigen::Triplet<double> currentTriplet;
    vector< Eigen::Triplet<double> > tripletList; // Fetched this from the Eigen website
    Eigen::SparseMatrix<double> sparseH;          // Or get this from systems?

    Eigen::MatrixXd eigenH;
    Eigen::Matrix2Xcd complexH;

    // Initializer
    Systems();
    Systems(int systemsize, double J, double h, bool armadillobool);
    Systems(int systemsize, double J, double h, bool armadillobool, bool dense);

    void set_mysector(int mysector);
    void change_system(double h);
    void change_system(double h, double J);
    void create_armadillo_matrix();
    void create_armadillo_matrix(int size);     // This is intended if we consider sectors
    void create_dense_Eigen_matrix();
    void create_dense_Eigen_matrix(int size);
    void create_complex_Eigen_matrix();
    void create_complex_Eigen_matrix(int size);

    // Basic spin operations
    int give_spin(int i, int a);
    int set_spin_up(int i, int a);
    int set_spin_down(int i, int a);
    int flip_spin(int i, int a);

    // Spin operators
    double szi(int i, int a);
    double szip1szi(int i, int a);
    int upip1downi(int i, int a);  // A simpler version of spip1smi (if-test outside of function)
    int downip1upi(int i, int a);  // A simpler version of smip1spi (if-test outside of function)

    // State-specific functions
    int number_of_up(int a);
    int number_of_down(int a);
    void checktestupdown(int j, int i);
    void sector0();
    void sector1_2();          // Odd name, perhaps...
    void find_sector_sparse();
    void find_sector_dense();
    void trim_sectorlist();

    //The systems
    void randomize();
    void set_hs(vector<double> hs_in);
    void set_hs_hom();
    void set_hs_alt();


    //Hamiltonians: Different kinds of systems
    void set_elements(int i, int b);           // Should I use this instead of having different functions for each type of system
    void palhuse_set_elements(int i, int b);
    // Sector Hamiltonians
    void palhuse_interacting_sectorHamiltonian_dense();   // Merge these?
    void palhuse_interacting_sectorHamiltonian_sparse();  //
    void palhuse_random_sectorHamiltonian_dense();        // Don't need these?
    void palhuse_random_sectorHamiltonian_sparse();       //
    void palhuse_diagonal_sectorHamiltonian();

    // Total Hamiltonians
    void palhuse_interacting_totalHamiltonian();
    void palhuse_diagonal_totalHamiltonian();


};

#endif // SYSTEMS_H

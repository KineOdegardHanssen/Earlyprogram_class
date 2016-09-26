#ifndef SYSTEMS_H
#define SYSTEMS_H
#include <iostream>
#include <vector>
#include <sector.h>
#include <cmath>
#include <armadillo>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>

using namespace std;

class Systems
{
public:

    unsigned long systemsize, no_of_states, mysector, number_of_hits;

    double h, J;

    bool dense, palhuse, sectorbool, armadillobool, testupip1downi, testdownip1upi;

    vector<double> hs;
    vector<unsigned long int> sectorlist;
    typedef Eigen::Triplet<double> T;

    // For armadillo solvers        // Okay, the program flow certainly is awkward...
    arma::mat armaH;

    // For Eigen solvers
    Eigen::Triplet<double> currentTriplet;
    vector< Eigen::Triplet<double> > tripletList; // Fetched this from the Eigen website
    Eigen::SparseMatrix<double> sparseH;          // Or get this from systems?

    //Eigen::MatrixXd eigenH;

    // Initializer
    Systems();
    Systems(unsigned long systemsize, double J, double h, bool armadillobool);
    Systems(unsigned long systemsize, double J, double h, bool armadillobool, bool dense);

    void set_mysector(unsigned long mysector);
    void change_system(double h);
    void change_system(double h, double J);
    void create_armadillo_matrix();
    void create_armadillo_matrix(unsigned long size);     // This is intended if we consider sectors
    //void create_dense_Eigen_matrix();
    //void create_dense_Eigen_matrix(unsigned long size);

    // Basic spin operations
    unsigned long give_spin(unsigned long i, unsigned long a);
    unsigned long set_spin_up(unsigned long i, unsigned long a);
    unsigned long set_spin_down(unsigned long i, unsigned long a);
    unsigned long flip_spin(unsigned long i, unsigned long a);

    // Spin operators
    double szi(unsigned long i, unsigned long a);
    double szip1szi(unsigned long i, unsigned long a);
    unsigned long upip1downi(unsigned long i, unsigned long a);  // A simpler version of spip1smi (if-test outside of function)
    unsigned long downip1upi(unsigned long i, unsigned long a);  // A simpler version of smip1spi (if-test outside of function)

    // State-specific functions
    unsigned long number_of_up(unsigned long a);
    unsigned long number_of_down(unsigned long a);
    void checktestupdown(unsigned long j, unsigned long i);
    void sector0();
    void find_sector_sparse();
    void find_sector_dense();
    void trim_sectorlist();

    //The systems
    void randomize();
    void set_hs_hom();
    void set_hs_alt();


    //Hamiltonians: Different kinds of systems
    void set_elements(unsigned long i, unsigned long b);           // Should I use this instead of having different functions for each type of system
    void palhuse_set_elements(unsigned long i, unsigned long b);
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

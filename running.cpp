#include "running.h"

Running::Running()    // Just in case I forget.
{
    armadillobool = true;
    cout << "Forgot to give armadillobool. Running for armadillobool=true." << endl;
}

Running::Running(bool armadillobool)    // Constructor for big runs
{
    this->armadillobool = armadillobool;
}

Running::Running(int systemsize, int maxit, double h, double J, double tolerance, bool armadillobool)   // Constructor best suited for one realization
{
    this->systemsize = systemsize;
    this->h = h;
    this->J = J;
    this->armadillobool = armadillobool;
}

void Running::homogenous_field(int systemsize, double h, double J)
{
    // Only works for the total Hamiltonian so far...
    system = Systems(systemsize, J, h, armadillobool);

    system.set_hs_hom();
    system.palhuse_interacting_totalHamiltonian();
    system.palhuse_diagonal_totalHamiltonian();

    Diagonalization eig_all = Diagonalization(system);

    if(armadillobool)
    {
        eig_all.using_armadillo();
        //diag.print_using_armadillo();
    }
    else
    {
        eig_all.lapack_directly();
    }

    Quantities quants = Quantities(maxit, tolerance, armadillobool, system, eig_all);
    for(int i=quants.li; i<=quants.lh; i++)   // Something is wrong here
    {
        cout << "Eigenvalue no." << i << "; Abs. diff. between first matrix elements: " << quants.ETH(i) << endl;
    }
}


void Running::alternating_field(int systemsize, double h, double J)
{
    // Similar implementation as for the homogenous field.

}

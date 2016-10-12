#include "running.h"

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


}


void Running::alternating_field(int systemsize, double h, double J)
{

}

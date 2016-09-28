#include <iostream>
#include <vector>
#include <systems.h>
#include <diagonalization.h>

using namespace std;

void system_total_hom(unsigned long systemsize, double h, double J, bool armabool, bool densebool);
void system_sector0_hom(unsigned long systemsize, double h, double J, bool armabool, bool densebool);

int main()
{       
    unsigned long systemsize = 2;
    double J = 1;
    double h = 1;
    bool armabool = true;
    bool densebool = true;

    //system_total_hom(systemsize, h, J, armabool, densebool);
    system_sector0_hom(systemsize, h, J, armabool, densebool);

    /*
    unsigned long number_of_hits = 0;
    unsigned long no_of_up = 0;
    unsigned long nstates = pow(2, systemsize);
    unsigned long mysector = nstates/2;
    for(unsigned long state=0; state<nstates; state++)
    {
        no_of_up = 0;
        for(unsigned long i=0; i<systemsize; i++)    no_of_up += ((state&(1<<i))>>i);
        cout << "no_of_up: " << no_of_up << "; mysector: " << mysector << "; state: " << state << endl;
        if(no_of_up==mysector)
        {
            cout << "state = " << state << endl;
            number_of_hits ++;
        } // End if
    } // End for
    cout << "number of hits = " << number_of_hits << endl;


    unsigned long a = 7;
    unsigned long no_of_up = 0;
    for(unsigned long i=0; i<systemsize; i++)    no_of_up += ((a&(1<<i))>>i);
    cout << no_of_up << endl;
    */

    //cout << testsystem.sparseH << endl;

    /*
    typedef Eigen::SparseMatrix<double> SpM;
    double timestart = clock();
    Eigen::EigenSolver<SpM> es(testsystem.sparseH);
    double timeend = clock();
    double comptime = timeend - timestart;

    cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
    cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;

    double printendtime = clock();

    double printtime = printendtime - timeend;

    cout << "Computation time: " << comptime << endl;
    cout << "Printing time: "    << printtime << endl;
    */

}

void system_total_hom(unsigned long systemsize, double h, double J, bool armabool, bool densebool)
{
    Systems testsystem = Systems(systemsize, J, h, armabool, densebool);

    testsystem.set_hs_hom();
    testsystem.palhuse_interacting_totalHamiltonian();
    testsystem.palhuse_diagonal_totalHamiltonian();

    for(unsigned long i=0; i<testsystem.no_of_states; i++)
    {
        for(unsigned long j=0; j<testsystem.no_of_states; j++)    cout << testsystem.armaH(i,j) << " ";
        cout << endl;
    }

    Diagonalization giveit = Diagonalization(testsystem);
    giveit.using_armadillo();
    giveit.print_using_armadillo();
}


void system_sector0_hom(unsigned long systemsize, double h, double J, bool armabool, bool densebool)
{
    Systems testsystem = Systems(systemsize, J, h, armabool, densebool);


    testsystem.set_hs_hom();
    testsystem.sector0();
    if(densebool==true)     testsystem.palhuse_interacting_sectorHamiltonian_dense();
    // Maybe I should just drop palhuse_interacting_sectorHamiltonian_sparse()
    testsystem.palhuse_diagonal_sectorHamiltonian();


    for(unsigned long i=0; i<testsystem.number_of_hits; i++)
    {
        for(unsigned long j=0; j<testsystem.number_of_hits; j++)    cout << testsystem.armaH(i,j) << " ";
        cout << endl;
    }

    Diagonalization giveit = Diagonalization(testsystem);
    giveit.using_armadillo();
    giveit.print_using_armadillo();

}

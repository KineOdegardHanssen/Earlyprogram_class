#include <iostream>
#include <vector>
#include <systems.h>
#include <diagonalization.h>

using namespace std;

int main()
{

    // Okay, so the Hamiltonian is flipped...
    unsigned long systemsize = 2;
    double J = 1;
    double h = 1;
    bool armabool = true;
    bool densebool = true;

    Systems testsystem = Systems(systemsize, J, h, armabool, densebool);

    /*
    testsystem.set_hs_hom();
    testsystem.sector0();
    testsystem.palhuse_interacting_sectorHamiltonian_sparse();
    testsystem.palhuse_diagonal_sectorHamiltonian();
    */

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


    /*
    for(unsigned long i=0; i<testsystem.no_of_states; i++)
    {
        for(unsigned long j=0; j<testsystem.no_of_states; j++)    cout << testsystem.eigenH(i,j) << " ";
        cout << endl;
    }


    Diagonalization giveit = Diagonalization(testsystem);
    giveit.print_dense_using_eigen();
    */

    cout << testsystem.sparseH << endl;

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


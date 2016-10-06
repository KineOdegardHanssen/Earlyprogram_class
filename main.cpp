#include <iostream>
#include <vector>
#include <systems.h>
#include <diagonalization.h>

using namespace std;

void system_total_hom(int systemsize, double h, double J, bool armabool, bool densebool);
void system_sector0_hom(int systemsize, double h, double J, bool armabool, bool densebool);
void system_sector1_2_hom(int systemsize, double h, double J, bool armabool, bool densebool);

int main()
{       
    const bool TRACE = true;
    int systemsize = 2;
    double J = 1;
    double h = 1;
    bool armabool = false;
    bool densebool = true;

    if(TRACE)    cout << "Well, at least I made it to main..." << endl;

    system_total_hom(systemsize, h, J, armabool, densebool);
    //system_sector0_hom(systemsize, h, J, armabool, densebool);
    //system_sector1_2_hom(systemsize, h, J, armabool, densebool);

}

void system_total_hom(int systemsize, double h, double J, bool armabool, bool densebool)
{
    Systems testsystem = Systems(systemsize, J, h, armabool, densebool);

    double set_system_start = clock();
    testsystem.set_hs_hom();
    testsystem.palhuse_interacting_totalHamiltonian();
    testsystem.palhuse_diagonal_totalHamiltonian();   
    double set_system_end = clock();
    double set_system_time = (set_system_end - set_system_start)/CLOCKS_PER_SEC;

    // /*
    for(int i=0; i<testsystem.number_of_hits; i++)
    {
        for(int j=0; j<testsystem.number_of_hits; j++)    cout << testsystem.eigenH(i,j) << " ";
        cout << endl;
    }

    cout << "I'm going to run diagonalization now" << endl;
    Diagonalization giveit = Diagonalization(testsystem);
    //giveit.using_armadillo();
    //giveit.print_using_armadillo();
    // */
    cout << "Diagonalization initiated" << endl;
    giveit.lapack_directly();

    cout << "Time needed to create the matrix: " << set_system_time << endl;
    cout << "Computation time using lapack directly: " << giveit.lapacktime << endl;

}


void system_sector0_hom(int systemsize, double h, double J, bool armabool, bool densebool)
{

    Systems testsystem = Systems(systemsize, J, h, armabool, densebool);


    /*
    vector<double> hs_in = vector<double>(systemsize);
    hs_in[0] = 0.1;
    hs_in[1] = 4;
    hs_in[2] = 7;
    hs_in[3] = 0.35;
    testsystem.set_hs(hs_in);
    */

    testsystem.set_hs_hom();
    testsystem.sector0();
    if(densebool==true)     testsystem.palhuse_interacting_sectorHamiltonian_dense();
    // Maybe I should just drop palhuse_interacting_sectorHamiltonian_sparse()
    testsystem.palhuse_diagonal_sectorHamiltonian();


    /*
    for(int i=0; i<testsystem.number_of_hits; i++)
    {
        for(int j=0; j<testsystem.number_of_hits; j++)    cout << testsystem.armaH(i,j) << " ";
        cout << endl;
    }
    */
    Diagonalization giveit = Diagonalization(testsystem);
    giveit.print_dense_using_eigen();
    //giveit.print_using_armadillo();


    //cout << "Computation time using armadillo: " << giveit.armatime << endl;

    //Diagonalization giveit2 = Diagonalization(testsystem);
    //giveit2.lapack_directly();

    //cout << "Computation time using lapack directly: " << giveit2.lapacktime << endl;
}

void system_sector1_2_hom(int systemsize, double h, double J, bool armabool, bool densebool)
{
    Systems testsystem = Systems(systemsize, J, h, armabool, densebool);

    double set_system_start = clock();
    testsystem.set_hs_hom();
    testsystem.sector1_2();
    if(densebool==true)     testsystem.palhuse_interacting_sectorHamiltonian_dense();
    // Maybe I should just drop palhuse_interacting_sectorHamiltonian_sparse()
    testsystem.palhuse_diagonal_sectorHamiltonian();
    double set_system_end = clock();
    double set_system_time = (set_system_end - set_system_start)/CLOCKS_PER_SEC;

    // /*
    for(int i=0; i<testsystem.number_of_hits; i++)
    {
        for(int j=0; j<testsystem.number_of_hits; j++)    cout << testsystem.eigenH(i,j) << " ";
        cout << endl;
    }

    cout << "I'm going to run diagonalization now" << endl;
    Diagonalization giveit = Diagonalization(testsystem);
    //giveit.using_armadillo();
    //giveit.print_using_armadillo();
    // */
    cout << "Diagonalization initiated" << endl;
    giveit.lapack_directly();

    cout << "Time needed to create the matrix: " << set_system_time << endl;
    cout << "Computation time using lapack directly: " << giveit.lapacktime << endl;
}

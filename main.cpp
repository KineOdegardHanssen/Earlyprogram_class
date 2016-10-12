#include <iostream>
#include <vector>
#include <systems.h>
#include <diagonalization.h>
#include <quantities.h>      // Do I really need this?
#include <running.h>

using namespace std;

void system_total_hom(int systemsize, double h, double J);
void system_sector0_hom(int systemsize, double h, double J);
void system_sector1_2_hom(int systemsize, double h, double J);

int main()
{       
    const bool TRACE = false;
    int systemsize = 5;
    int maxit = 1e7;
    double J = 1;
    double h = 1;
    double tolerance = 1e-10;
    bool armabool = true;

    if(TRACE)    cout << "Well, at least I made it to main..." << endl;

    //system_total_hom(systemsize, h, J);
    //system_sector0_hom(systemsize, h, J);
    system_sector1_2_hom(systemsize, h, J);

    //Running runit = Running(systemsize, maxit, h, J, tolerance, armabool);
    //runit.homogenous_field(systemsize, h, J);  // Should I change this call?
}

void system_total_hom(int systemsize, double h, double J)
{
    bool armabool = true;
    Systems testsystem = Systems(systemsize, J, h, armabool);

    double set_system_start1 = clock();
    testsystem.set_hs_hom();
    testsystem.palhuse_interacting_totalHamiltonian();
    testsystem.palhuse_diagonal_totalHamiltonian();   
    double set_system_end1 = clock();
    double set_system_time1 = (set_system_end1 - set_system_start1)/CLOCKS_PER_SEC;

    /*
    for(int i=0; i<testsystem.number_of_hits; i++)
    {
        for(int j=0; j<testsystem.number_of_hits; j++)    cout << testsystem.armaH(i,j) << " ";
        cout << endl;
    }
    cout << endl << endl;
    */

    cout << "Suspiciously small non-zero armadillo elements:" << endl;
    double ay;
    for(int i=0; i<testsystem.number_of_hits; i++)
    {
        for(int j=0; j<testsystem.number_of_hits; j++)
        {
            ay = testsystem.armaH(i,j);
            if(ay!=0 && abs(ay)<0.25)            cout << "Element non-zero: " << ay  << endl;
        }

    }
   cout << endl;
    /*
    Diagonalization giveit = Diagonalization(testsystem);
    giveit.using_armadillo();
    giveit.print_using_armadillo();
    */

    armabool = false;
    Systems testsystem2 = Systems(systemsize, J, h, armabool);

    double set_system_start2 = clock();
    testsystem2.set_hs_hom();
    testsystem2.palhuse_interacting_totalHamiltonian();
    testsystem2.palhuse_diagonal_totalHamiltonian();
    double set_system_end2 = clock();
    double set_system_time2 = (set_system_end2 - set_system_start2)/CLOCKS_PER_SEC;

    /*
    for(int i=0; i<testsystem.number_of_hits; i++)
    {
        for(int j=0; j<testsystem.number_of_hits; j++)    cout << testsystem2.eigenH(i,j) << " ";
        cout << endl;
    }
    */

    cout << "Suspiciously small non-zero Eigen elements:" << endl;
    for(int i=0; i<testsystem.number_of_hits; i++)
    {
        for(int j=0; j<testsystem.number_of_hits; j++)
        {
            ay = testsystem2.eigenH(i,j);
            if(ay!=0 && abs(ay)<0.25)            cout << "Element non-zero: " << ay  << endl;
        }

    }
   cout << endl;

    cout << "Time needed to create the matrix with armadillo: " << set_system_time1 << endl;
    cout << "Time needed to create the matrix with Eigen: " << set_system_time2 << endl;

    /*
    Diagonalization giveit2 = Diagonalization(testsystem2);
    giveit2.print_dense_using_eigen();

    giveit2.lapack_directly();

    cout << "Time needed to create the matrix with armadillo: " << set_system_time1 << endl;
    cout << "Time needed to create the matrix with Eigen: " << set_system_time2 << endl;
    cout << "Computation time using armadillo: " << giveit.armatime << endl;
    cout << "Computation time using Eigen: " << giveit2.eigen_time << endl;
    cout << "Computation time using lapack directly: " << giveit2.lapacktime << endl;
    */
}


void system_sector0_hom(int systemsize, double h, double J)
{

    bool armabool = true;
    Systems testsystem = Systems(systemsize, J, h, armabool);

    double set_system_start1 = clock();
    testsystem.set_hs_hom();
    testsystem.sector0();
    testsystem.palhuse_interacting_sectorHamiltonian_dense();
    testsystem.palhuse_diagonal_sectorHamiltonian();
    double set_system_end1 = clock();
    double set_system_time1 = (set_system_end1 - set_system_start1)/CLOCKS_PER_SEC;


    Diagonalization giveit = Diagonalization(testsystem);
    giveit.using_armadillo();
    giveit.print_using_armadillo();

    armabool = false;
    Systems testsystem2 = Systems(systemsize, J, h, armabool);

    double set_system_start2 = clock();
    testsystem2.set_hs_hom();
    testsystem2.sector0();
    testsystem2.palhuse_interacting_sectorHamiltonian_dense();
    testsystem2.palhuse_diagonal_sectorHamiltonian();
    double set_system_end2 = clock();
    double set_system_time2 = (set_system_end2 - set_system_start2)/CLOCKS_PER_SEC;

    Diagonalization giveit2 = Diagonalization(testsystem2);
    giveit2.print_dense_using_eigen();

    giveit2.lapack_directly();

    cout << "Time needed to create the matrix with armadillo: " << set_system_time1 << endl;
    cout << "Time needed to create the matrix with Eigen: " << set_system_time2 << endl;
    cout << "Computation time using armadillo: " << giveit.armatime << endl;
    cout << "Computation time using Eigen: " << giveit2.eigen_time << endl;
    cout << "Computation time using lapack directly: " << giveit2.lapacktime << endl;

    cout << "sectorlist using arma: " << endl;
    for(int i=0; i<testsystem.number_of_hits; i++)        cout << testsystem.sectorlist[i] << "   ";
    cout << endl;
    cout << "sectorlist using Eigen: ";
    for(int i=0; i<testsystem.number_of_hits; i++)        cout << testsystem2.sectorlist[i] << "   ";
    cout << endl;
}

void system_sector1_2_hom(int systemsize, double h, double J)
{
    bool armabool = true;
    Systems testsystem = Systems(systemsize, J, h, armabool);

    double set_system_start1 = clock();
    testsystem.set_hs_hom();
    testsystem.sector1_2();
    //testsystem.set_mysector(1); // Looking at sectot S^tot_z = -1/2
    testsystem.palhuse_interacting_sectorHamiltonian_dense();
    testsystem.palhuse_diagonal_sectorHamiltonian();
    double set_system_end1 = clock();
    double set_system_time1 = (set_system_end1 - set_system_start1)/CLOCKS_PER_SEC;

    Diagonalization giveit = Diagonalization(testsystem);
    giveit.using_armadillo();
    giveit.print_using_armadillo();
    /*
    for(int i=0; i<testsystem.number_of_hits; i++)
    {
        for(int j=0; j<testsystem.number_of_hits; j++)    cout << testsystem.armaH(i,j) << " ";
        cout << endl;
    }
    */

    /*
    cout << "Suspiciously small non-zero armadillo elements:" << endl;
    double ay;
    for(int i=0; i<testsystem.number_of_hits; i++)
    {
        for(int j=0; j<testsystem.number_of_hits; j++)
        {
            ay = testsystem.armaH(i,j);
            if(ay!=0 && abs(ay)<0.25)            cout << "Element non-zero: " << ay  << endl;
        }

    }
   cout << endl;
   */

    armabool = false;
    Systems testsystem2 = Systems(systemsize, J, h, armabool);

    double set_system_start2 = clock();
    testsystem2.set_hs_hom();
    testsystem2.sector1_2();
    //testsystem2.set_mysector(1); // Looking at sectot S^tot_z = -1/2
    testsystem2.palhuse_interacting_sectorHamiltonian_dense();
    testsystem2.palhuse_diagonal_sectorHamiltonian();
    double set_system_end2 = clock();
    double set_system_time2 = (set_system_end2 - set_system_start2)/CLOCKS_PER_SEC;

    /*
    for(int i=0; i<testsystem2.number_of_hits; i++)
    {
        for(int j=0; j<testsystem2.number_of_hits; j++)    cout << testsystem2.eigenH(i,j) << " ";
        cout << endl;
    }


    cout << "Suspiciously small non-zero Eigen elements:" << endl;
    for(int i=0; i<testsystem2.number_of_hits; i++)
    {
        for(int j=0; j<testsystem2.number_of_hits; j++)
        {
            ay = testsystem2.eigenH(i,j);
            if(ay!=0 && abs(ay)<0.25)            cout << "Element non-zero: " << ay  << endl;
        }

    }
   cout << endl;
   */

    Diagonalization giveit2 = Diagonalization(testsystem2);
    //giveit2.print_dense_using_eigen();

    giveit2.lapack_directly();

    cout << "Time needed to create the matrix with armadillo: " << set_system_time1 << endl;
    cout << "Time needed to create the matrix with Eigen: " << set_system_time2 << endl;
    //cout << "Computation time using armadillo: " << giveit.armatime << endl;
    //cout << "Computation time using Eigen: " << giveit2.eigen_time << endl;
    //cout << "Computation time using lapack directly: " << giveit2.lapacktime << endl;

    cout << "sectorlist using arma: " << endl;
    for(int i=0; i<testsystem.number_of_hits; i++)        cout << testsystem.sectorlist[i] << "   ";
    cout << endl;
    cout << "sectorlist using Eigen: ";
    for(int i=0; i<testsystem.number_of_hits; i++)        cout << testsystem2.sectorlist[i] << "   ";
    cout << endl;

    /*
    //Eigen::VectorXd ev_diff(testsystem.number_of_hits);
    double average = 0;
    double thediff = 0;
    double largest_elem = 0;
    double no_not_zero = 0;
    for(int i=0; i<testsystem.number_of_hits; i++)
    {
        thediff = abs(giveit.eigenvalues_armadillo(i) - giveit2.eigenvalues_H(i));
        //ev_diff(i) = thediff;
        average += thediff;
        if(largest_elem<thediff)    largest_elem = thediff;
        if(thediff!=0)              no_not_zero++;
    }

    average = average/testsystem.number_of_hits;
    cout << "Max difference between eigenvalues: " << largest_elem << endl;
    cout << "Average absolute difference between the eigenvalues found by armadillo and LAPACK: " << average <<endl;
    cout << "Number of eigenvalues where armadillo and LAPACK differed:" << no_not_zero << endl;
    cout << "Total number of eigenvalues: " << testsystem.number_of_hits << endl;
    cout << "Percent of eigenvalues not exactly equal: " << (100*no_not_zero/testsystem.number_of_hits) << endl;
    */
}

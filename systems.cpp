#include "systems.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <armadillo>
//#include <random>             // Just leaving this out for now...
#include <Eigen/Dense>
//#include <Eigen/Sparse>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>


using namespace std;

Systems::Systems()
{
}

Systems::Systems(int systemsize, double J, double h, bool armadillobool)
{
    this->systemsize = systemsize;                // Is this notation really neccessary? Look it up
    this->no_of_states = pow(2, systemsize);
    this->J=J;
    this->h=h;
    hs = vector<double>(no_of_states);
    this->armadillobool = armadillobool;
    palhuse = true;                               // Default setting
    //randomize();
}


void Systems::set_mysector(int mysector)
{
    this->mysector = mysector;
    find_sector_dense();
}

void Systems::change_system(double h)
{
    hs = vector<double> (no_of_states);
    number_of_hits = 0;
    this->h = h;
    randomize();
}

void Systems::change_system(double h, double J)
{
    hs = vector<double> (no_of_states);
    number_of_hits = 0;
    this->h = h;
    this->J= J;
}

void Systems::create_armadillo_matrix()
{   // Creates an empty matrix, ready to use, in those cases where it is not too big.
    armaH = arma::mat(no_of_states, no_of_states);
    for(int i=0; i<no_of_states; i++)
    {
        for(int j=0; j<no_of_states; j++)    armaH(i,j)= 0.0;
    }
    //cout << "Arma matrix set. Max index no = " << no_of_states-1 << endl;
}

void Systems::create_armadillo_matrix(int size)
{   // Creates an empty matrix, ready to use, in those cases where it is not too big.
    armaH = arma::mat(size, size);
    for(int i=0; i<size; i++)        // This does not seem to be neccessary all the time, but it is safer.
    {
        for(int j=0; j<size; j++)    armaH(i,j)= 0.0;
    }
    //cout << "Arma matrix set. Max index no = " << size-1 << endl;
}


void Systems::create_dense_Eigen_matrix()
{   // Creates an empty matrix, ready to use, in those cases where it is not too big.
    //const bool TRACE = false;
    //if(TRACE)    cout << "I am going to set the size of an eigenmatrix" << endl;
    eigenH = Eigen::MatrixXd(no_of_states, no_of_states);
    for(int i=0; i<no_of_states; i++)
    {
        for(int j=0; j<no_of_states; j++)    eigenH(i,j)= 0.0;
    }
    //cout << "Eigen matrix set. Max index no. = " << no_of_states-1 << endl;
}

void Systems::create_dense_Eigen_matrix(int size)
{   // Creates an empty matrix, ready to use, in those cases where it is not too big.
    //const bool TRACE = false;
    //if(TRACE)    cout << "I am going to set the size of an eigenmatrix" << endl;
    eigenH = Eigen::MatrixXd(size, size);
    for(int i=0; i<size; i++)
    {
        for(int j=0; j<size; j++)    eigenH(i,j)= 0.0;
    }
    //cout << "Eigen matrix set. Max index no. = " << size-1 << endl;
}





//---------------------------------------BASIC SPIN OPERATIONS--------------------------------------------//

int Systems::give_spin(int i, int a)
{
    return ((a&(1<<i))>>i);
}

int Systems::set_spin_up(int i, int a)
{
    return (a |(1<<i));
}

int Systems::set_spin_down(int i, int a)
{
    return ~((~a)|(1<<i));
}

int Systems::flip_spin(int i, int a)      // This is not in use
{
    return (a ^(1<<i));
}


//--------------------------------------------SPIN OPERATORS----------------------------------------------//
//-------------------------------------------------/Sz/---------------------------------------------------//


double Systems::szi(int i, int a)
{
    if(give_spin(i,a)==0)    return -0.5;
    else                     return 0.5;
}

double Systems::szip1szi(int i, int a)
{
    double firstspin = szi(i, a);
    double secondspin = 0;
    if(i==(systemsize-1))    secondspin = szi(0,a);
    else                     secondspin = szi((i+1),a);
    return firstspin*secondspin;
}


//-------------------------------------------/S+ and S-/--------------------------------------------------//


int Systems::upip1downi(int i, int a)
{
    int a2 = 0;
    if(i==(systemsize-1))            a2 = set_spin_up(0, a);         // Periodic boundary conditions
    else if(i<(systemsize-1))        a2 = set_spin_up((i+1), a);     //
    else                             cout << "Check your indices, woman!!" << endl;
    return set_spin_down(i, a2);
}

int Systems::downip1upi(int i, int a)
{
    int a2 = 0;
    if(i==(systemsize-1))            a2 = set_spin_down(0, a);       // Periodic boundary conditions
    else if(i<(systemsize-1))        a2 = set_spin_down((i+1), a);
    else                             cout << "Check your indices, woman!!" << endl;
    return set_spin_up(i, a2);
}


//-----------------------------------------STATE-SPECIFIC FUNCTIONS---------------------------------------//
// Both of these functions work. It might be excessive to have number_of_down(a,systemsize) as a separate
// funtction.
int Systems::number_of_up(int a)
{
    int no_of_up = 0;
    for(int i=0; i<systemsize; i++)    no_of_up += give_spin(i, a);
    return no_of_up;
}

int Systems::number_of_down(int a)
{
    return (systemsize - number_of_up(a));
}

// function for s0 sector

void Systems::checktestupdown(int j,int i)
{
    int helpnumber = 0;


    if(j!=(systemsize-1))
    {
        testupip1downi = (give_spin(j,i)==1) && ( (give_spin((j+1),i))==0);
        testdownip1upi = (give_spin(j,i)==0) && ( (give_spin((j+1),i))==1);
    }
    else
    {
        testupip1downi = (give_spin(j,i)==1) && ( (give_spin(helpnumber,i))==0);
        testdownip1upi = (give_spin(j,i)==0) && ( (give_spin(helpnumber,i))==1);
    }
}


void Systems::sector0()
{
    if(systemsize%2==0)
    {
        mysector = systemsize/2;
        find_sector_dense();  // Change these commands.
    }
    else    cout << "A system of an odd no. of particles have no S_tot = 0 states. Try another sector." << endl;
}

void Systems::sector1_2()
{   // System with one more spin up than spin down. This is really sector1_2.
    if(systemsize%2==1)
    {
        mysector = (systemsize+1)/2;
        find_sector_dense();  // Change these commands.
    }
    else    cout << "A system of an even no. of particles have no S_tot = 1/2 states. Try another sector." << endl;
}


void Systems::find_sector_dense()
{   // Should keep some information on what states are in the list. Well, that is all in sectorlist.
    // And that is not destroyed. Retrieve it in main, perhaps.
    sectorbool = true;
    number_of_hits = 0;
    sectorlist = vector<int>(no_of_states);
    for(int state=0; state<no_of_states; state++)
    {
        if(number_of_up(state)==mysector)
        {
            sectorlist[number_of_hits] = state;
            number_of_hits ++;
        } // End if
    } // End for
    if(number_of_hits > 0)                   trim_sectorlist();
    else
    {
        sectorlist = vector<int>(1);
        sectorlist[0] = 0;
        number_of_hits = 1;
    }
}

void Systems::trim_sectorlist()
{   // A function that cuts sectorlist off as it detects a zero
    vector<int> sectorlist_short = vector<int>(number_of_hits);
    for(int i=0; i<number_of_hits; i++)     sectorlist_short[i] = sectorlist[i];
        sectorlist = sectorlist_short;
}

//--------------------------------------------//THE SYSTEM//----------------------------------------------//
//-----------------------------------------------/GENERAL/------------------------------------------------//


void Systems::randomize()
{
    /*
    std::default_random_engine generator;        // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution(-h,h);
    for(int i=0; i<no_of_states; i++)
    {
        hs[i] = distribution(generator);   // This should do it
    } // End for-loop
    */
}  // End function randomize

void Systems::set_hs(vector<double> hs_in)
{   // This function may be useful for testing. NB: Must make sure len(hs_in) = no_of_states.
    for(unsigned int i=0; i<no_of_states; i++)        hs[i] = hs_in[i];
}

void Systems::set_hs_hom()
{
    for(int i=0; i<no_of_states; i++)
    {
        hs[i] = h;   // This should do it
    } // End for-loop
}  // End function randomize

void Systems::set_hs_alt()
{
    for(int i=0; i<no_of_states; i++)
    {
        hs[i] = h*pow(-1,i);   // This should do it
    } // End for-loop
}  // End function randomize


//-------------------------------------------/HAMILTONIANS/--------------------------------------------//

void Systems::set_elements(int i, int b)
{
    if(palhuse==true)    palhuse_set_elements(i, b);
}

void Systems::palhuse_set_elements(int i, int b)
{
    int index1 = 0;
    int index2 = 0;
    if(sectorbool==false)
    {
        index1 = no_of_states - (b+1);
        index2 = no_of_states - (i+1);
    }
    else
    {   // Must correct this later. Do not want an upside-down matrix. Or? Not a big problem if I retrieve
        // sectorlist. That holds for all systems with the same no. of. particles.
        index1 = i;
        index2 = b;
        // Blimey. This redistribution of states is trickier than I thought.
    }


    double element = 0;
    if(systemsize == 2)        element = J;  // Special case. BCs and small
    else                       element = 0.5*J;

    if(armadillobool == true)
    {
        armaH(index1,index2) = element;
        armaH(index2,index1) = element;
    }
    else
    {
        eigenH(index1, index2) = element;
        eigenH(index2, index1) = element;

    }
}
//----------------------------------------SECTOR HAMILTONIANS---------------------------------------------//

void Systems::palhuse_interacting_sectorHamiltonian_dense()
{
    const bool TRACE = false;
    if(TRACE)    cout << "At least I am in palhuse_interacting_sectorHamiltonian" << endl;

    if(armadillobool)        create_armadillo_matrix(number_of_hits);

    if(armadillobool==false)        create_dense_Eigen_matrix(number_of_hits);

    int a = 0;
    int b = 0;
    for(int i=0; i<number_of_hits; i++)
    {
        a = sectorlist[i];
        for(int j=0; j<systemsize; j++)

        {   // Should I include a while loop, or something?

            checktestupdown(j,a);
            if(testupip1downi==true)
            {
                b = upip1downi(j, a);
                for(int k = 0; k<number_of_hits; k++)
                {
                    if(b==sectorlist[k])  set_elements(i, k);
                }
            }
            if(testdownip1upi==true)
            {
                b = downip1upi(j, a);
                for(int k = 0; k<number_of_hits; k++)
                {
                    if(b==sectorlist[k])  set_elements(i, k);
                }
            }
        } // End for j
    } // End for i
    if(TRACE)    cout << "Exiting palhuse_interacting_sectorHamiltonian" << endl;
}

void Systems::palhuse_diagonal_sectorHamiltonian()
{
    const bool TRACE = false;
    // Do something like index = no_of_states - i that works for this one.
    // For now, the matrix is upside down, but we have gotten a list of its entries: sectorlist.
    double element = 0;
    for(int i=0; i<number_of_hits; i++)
    {
        element = 0;
        int a = sectorlist[i];
        for(int j=0; j<systemsize; j++)  element += hs[j]*szi(j, a) + J*szip1szi(j,a);
        if(armadillobool == true)                            armaH(i,i) = element;
        else if(armadillobool == false)                      eigenH(i,i) = element;
    } // End for-loop over i
    if(TRACE)    cout << "Exiting palhuse_diagonal_sectorHamiltonian" << endl;
} // End function palhuse_random_sectorHamiltonian_dense

//-------------------------------------------TOTAL HAMILTONIAN--------------------------------------------//

void Systems::palhuse_interacting_totalHamiltonian()
{
    if(armadillobool)                                               create_armadillo_matrix();
    if(armadillobool==false)                                        create_dense_Eigen_matrix();
    sectorbool=false;

    number_of_hits = no_of_states; // Have a test if(sectorbool) and set this, maybe along with some other statements, and drop separate functions for total and sector Hamiltonians?

    int b = 0;
    for(int i=0; i<no_of_states; i++)
    {   // i is our state
        for(int j=0; j<systemsize; j++)
        {
            checktestupdown(j,i);
            if(testupip1downi==true)
            {
                b = upip1downi(j, i);
                set_elements(i,b);
            }
            if(testdownip1upi==true)
            {
                b = downip1upi(j, i);
                set_elements(i,b);
            }
        } // End for j
    } // End for i
} // End function



void Systems::palhuse_diagonal_totalHamiltonian()
{
    int index = 0;
    double element = 0;
    for(int i=0; i<no_of_states; i++)
    {
        element = 0;
        // The 2-particle case is special again...
        for(int j=0; j<systemsize; j++)  element += hs[j]*szi(j, i) + J*szip1szi(j,i);
        index = no_of_states - (i+1);
        if(armadillobool==true)                            armaH(index,index) = element;
        else if(armadillobool==false)                      eigenH(index,index) = element;
    } // End for-loop over i
    if(armadillobool==false)                                                  eigenH.cast<double>();
} // End function palhuse_random_sectorHamiltonian_dense

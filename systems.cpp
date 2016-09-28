#include "systems.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <armadillo>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>

/* Use

SparseMatrix<double> sm1(1000,1000);
SparseMatrix<std::complex<double>,RowMajor> sm2;


std::vector< Eigen::Triplet<double> > tripletList;
tripletList.reserve(estimation_of_entries);
// -- Fill tripletList with nonzero elements...
sm1.setFromTriplets(TripletList.begin(), TripletList.end());

...Bruke push_back ...?

make_compress ...
*/

// Find some way to reset tripleList or removing the last values.

using namespace std;

Systems::Systems()
{
}

Systems::Systems(unsigned long systemsize, double J, double h, bool armadillobool)
{
    this->systemsize = systemsize;                // Is this notation really neccessary? Look it up
    this->no_of_states = pow(2, systemsize);
    this->J=J;
    this->h=h;
    hs = vector<double>(no_of_states);
    this->dense = false;
    this->armadillobool = armadillobool;
    palhuse = true;                               // Default setting
    randomize();
}

Systems::Systems(unsigned long systemsize, double J, double h, bool armadillobool, bool dense)
{
    this->systemsize = systemsize;
    this->no_of_states = pow(2, systemsize);
    this->J = J;
    this->h = h;
    hs = vector<double>(no_of_states);
    this->dense = dense;
    this->armadillobool = armadillobool;
    palhuse = true;                               // Default setting
    randomize();
}


void Systems::set_mysector(unsigned long mysector)
{
    this->mysector = mysector;
    if(dense==true)    find_sector_dense();
    else               find_sector_sparse();
}

void Systems::change_system(double h)
{
    hs = vector<double> (no_of_states);
    number_of_hits = 0;
    this->h = h;
    randomize();
    // Removing the diagonal entries in tripleList, as these contain the dependence of the h's
    if(dense==false)    for(unsigned int i=0; i<no_of_states; i++)        tripletList.pop_back();
    // Because the diagonal part of the matrix is always set last.
    // Using dense armadillo matrices, these are just overwritten.
    // That may not be the case for dense matrices using Eigen.
}

void Systems::change_system(double h, double J)
{
    hs = vector<double> (no_of_states);
    number_of_hits = 0;
    this->h = h;
    this->J= J;
    if(dense==false)
    {
        while(tripletList.empty()==false)     tripletList.pop_back();   // Resetting tripletList
        // Maybe find a more effective way to do this? Just declare it again?
        // Don't think I will actually use this function.
    }
}

void Systems::create_armadillo_matrix()
{   // Creates an empty matrix, ready to use, in those cases where it is not too big.
    armaH = arma::mat(no_of_states, no_of_states);
    for(unsigned long i=0; i<no_of_states; i++)
    {
        for(unsigned long j=0; j<no_of_states; j++)    armaH(i,j)= 0.0;
    }
}

void Systems::create_armadillo_matrix(unsigned long size)
{   // Creates an empty matrix, ready to use, in those cases where it is not too big.
    armaH = arma::mat(size, size);
    for(unsigned long i=0; i<size; i++)
    {
        for(unsigned long j=0; j<size; j++)    armaH(i,j)= 0.0;
    }
}

/*
void Systems::create_dense_Eigen_matrix()
{   // Creates an empty matrix, ready to use, in those cases where it is not too big.
    eigenH = Eigen::MatrixXd(no_of_states, no_of_states);
    for(unsigned long i=0; i<no_of_states; i++)
    {
        for(unsigned long j=0; j<no_of_states; j++)    eigenH(i,j)= 0.0;
    }
}

void Systems::create_dense_Eigen_matrix(unsigned long size)
{   // Creates an empty matrix, ready to use, in those cases where it is not too big.
    eigenH = Eigen::MatrixXd(size, size);
    for(unsigned long i=0; i<size; i++)
    {
        for(unsigned long j=0; j<size; j++)    eigenH(i,j)= 0.0;
    }
}
*/

//---------------------------------------BASIC SPIN OPERATIONS--------------------------------------------//

unsigned long Systems::give_spin(unsigned long i, unsigned long a)
{
    return ((a&(1<<i))>>i);
}

unsigned long Systems::set_spin_up(unsigned long i, unsigned long a)
{
    return (a |(1<<i));
}

unsigned long Systems::set_spin_down(unsigned long i, unsigned long a)
{
    return ~((~a)|(1<<i));
}

unsigned long Systems::flip_spin(unsigned long i, unsigned long a)      // This is not in use
{
    return (a ^(1<<i));
}


//--------------------------------------------SPIN OPERATORS----------------------------------------------//
//-------------------------------------------------/Sz/---------------------------------------------------//


double Systems::szi(unsigned long i, unsigned long a)
{
    if(give_spin(i,a)==0)    return -0.5;
    else                     return 0.5;
}

double Systems::szip1szi(unsigned long i, unsigned long a)
{
    double firstspin = szi(i, a);
    double secondspin = 0;
    if(i==(systemsize-1))    secondspin = szi(0,a);
    else                     secondspin = szi((i+1),a);
    return firstspin*secondspin;
}


//-------------------------------------------/S+ and S-/--------------------------------------------------//


unsigned long Systems::upip1downi(unsigned long i, unsigned long a)
{
    unsigned long a2 = 0;
    if(i==(systemsize-1))            a2 = set_spin_up(0, a);         // Periodic boundary conditions
    if(i<(systemsize-1))             a2 = set_spin_up((i+1), a);     //
    if(i>(systemsize-1))             cout << "Check your indices, woman!!" << endl;
    return set_spin_down(i, a2);
}

unsigned long Systems::downip1upi(unsigned long i, unsigned long a)
{
    unsigned long a2 = 0;
    if(i==(systemsize-1))            a2 = set_spin_down(0, a);       // Periodic boundary conditions
    if(i<(systemsize-1))             a2 = set_spin_down((i+1), a);
    if(i>(systemsize-1))             cout << "Check your indices, woman!!" << endl;
    return set_spin_up(i, a2);
}


//-----------------------------------------STATE-SPECIFIC FUNCTIONS---------------------------------------//
// Both of these functions work. It might be excessive to have number_of_down(a,systemsize) as a separate
// funtction.
unsigned long Systems::number_of_up(unsigned long a)
{
    unsigned long no_of_up = 0;
    for(unsigned long i=0; i<systemsize; i++)    no_of_up += give_spin(i, a);
    return no_of_up;
}

unsigned long Systems::number_of_down(unsigned long a)
{
    return (systemsize - number_of_up(a));
}

// function for s0 sector

void Systems::checktestupdown(unsigned long j,unsigned long i)
{
    unsigned long helpnumber = 0;


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

        if(dense==true)        find_sector_dense();  // Change these commands.
        else                   find_sector_sparse();
    }
    else    cout << "A system of 3 particles have no S_tot = 0 states. Try another sector." << endl;
}

//Something is rotten in the state of sectors.
void Systems::find_sector_sparse()
{   // Should append to some list in some way... Armadillo? Eigen? Eigen.
    //unsigned long number = 0;
    sectorbool = true;
    number_of_hits = 0;
    sectorlist = vector<unsigned long int>(no_of_states);
    for(unsigned long state=0; state<no_of_states; state++)
    {
        if(number_of_up(state)==mysector) // Well, something... Append to some list
        {
            number_of_hits ++;
            sectorlist[state] = number_of_hits;
        } // End if
    } // End while
}

void Systems::find_sector_dense()
{   // Should keep some information on what states are in the list. Well, that is all in sectorlist.
    // And that is not destroyed. Retrieve it in main, perhaps.
    sectorbool = true;
    number_of_hits = 0;
    sectorlist = vector<unsigned long>(no_of_states);
    for(unsigned long state=0; state<no_of_states; state++)
    {
        if(number_of_up(state)==mysector)
        {
            sectorlist[number_of_hits] = state;
            cout << "sectorlist[number_of_hits] = " << sectorlist[number_of_hits] << endl;
            number_of_hits ++;
        } // End if
    } // End for
    if(number_of_hits > 0)                   trim_sectorlist();
    else
    {
        sectorlist = vector<unsigned long>(1);
        sectorlist[0] = 0;
        number_of_hits = 1;
    }
}

void Systems::trim_sectorlist()
{   // A function that cuts sectorlist off as it detects a zero
    vector<unsigned long> sectorlist_short = vector<unsigned long>(number_of_hits);
    for(unsigned long i=0; i<number_of_hits; i++)     sectorlist_short[i] = sectorlist[i];
        sectorlist = sectorlist_short;
}

//--------------------------------------------//THE SYSTEM//----------------------------------------------//
//-----------------------------------------------/GENERAL/------------------------------------------------//


void Systems::randomize()
{
    std::default_random_engine generator;        // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution(-h,h);
    for(unsigned long i=0; i<no_of_states; i++)
    {
        hs[i] = distribution(generator);   // This should do it
    } // End for-loop
}  // End function randomize

void Systems::set_hs_hom()
{
    for(unsigned long i=0; i<no_of_states; i++)
    {
        hs[i] = h;   // This should do it
    } // End for-loop
}  // End function randomize

void Systems::set_hs_alt()
{
    for(unsigned long i=0; i<no_of_states; i++)
    {
        hs[i] = h*pow(-1,i);   // This should do it
    } // End for-loop
}  // End function randomize


//-------------------------------------------/HAMILTONIANS/--------------------------------------------//

void Systems::set_elements(unsigned long i, unsigned long b)
{
    if(palhuse==true)    palhuse_set_elements(i, b);
}

void Systems::palhuse_set_elements(unsigned long i, unsigned long b)
{
    unsigned long index1 = 0;
    unsigned long index2 = 0;
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

    if(dense==false)
    {
        if(systemsize == 2)             currentTriplet = T(index1,index2,0.5*element); // ... ?
        else                            currentTriplet = T(index1,index2,element);
        tripletList.push_back(currentTriplet);
    }

    if(dense==true)
    {
        if(armadillobool == true)
        {
            armaH(index1,index2) = element;
            armaH(index2,index1) = element;
        }  // End if-test dense
        /*
        else   // Maybe set a Eigenbool sometime also?
        {  // Check if this returns the correct result
            int b_int = b & INT_MAX;
            int i_int = i & INT_MAX;
            eigenH(b_int,i_int) = element;
            eigenH(i_int,b_int) = element;
        }*/
    }
}
//----------------------------------------SECTOR HAMILTONIANS---------------------------------------------//

void Systems::palhuse_interacting_sectorHamiltonian_dense()
{
    if(armadillobool == true)    create_armadillo_matrix(number_of_hits);
    //else                         create_dense_Eigen_matrix(number_of_hits);

    unsigned long a = 0;
    unsigned long b = 0;
    for(unsigned long i=0; i<number_of_hits; i++)
    {
        a = sectorlist[i];
        for(unsigned long j=0; j<systemsize; j++)

        {   // Should I include a while loop, or something?

            checktestupdown(j,a);
            if(testupip1downi==true)
            {
                b = upip1downi(j, a);
                for(unsigned long k = 0; k<number_of_hits; k++)
                {
                    if(b==sectorlist[k])  set_elements(i, k);
                }
            }
            if(testdownip1upi==true)
            {
                b = downip1upi(j, a);
                for(unsigned long k = 0; k<number_of_hits; k++)
                {
                    if(b==sectorlist[k])  set_elements(i, k);
                }
            }
        } // End s+_{i+1}s-_{i} part
    } // End for j
}  // End for i



void Systems::palhuse_interacting_sectorHamiltonian_sparse()
{
    // The Hamiltonian is hermittian, i.e. symmetric about the diagonal, so we need not change the indices.
    unsigned long length = number_of_hits << 3;

    currentTriplet = T(0,0,0);
    tripletList.reserve(length);     // Is this too big?

    unsigned long a = 0;
    unsigned long b = 0;

    for(unsigned long i=0; i<no_of_states; i++)  // This is excessive, but whatever...
    {   // Should I include a while loop, or something?
        a = sectorlist[i];  // Get no contribution from the zero elements, but they do steal computation time
        // spi1psmi, upip1downi
        for(unsigned long j=0; j<systemsize; j++)
        {           
            checktestupdown(j,a);
            if(testupip1downi==true)
            {
                b = upip1downi(j, a);
                set_elements(a,b);
            }
            if(testdownip1upi==true)
            {
                b = downip1upi(j, a);
                set_elements(a,b);
            }
        } // End s+_{i+1}s-_{i} part
    } // End for j
}  // End for i


void Systems::palhuse_diagonal_sectorHamiltonian()
{
    // Do something like index = no_of_states - i that works for this one.
    // For now, the matrix is upside down, but we have gotten a list of its entries: sectorlist.
    double element = 0;
    for(unsigned long i=0; i<number_of_hits; i++)
    {
        element = 0;
        unsigned long a = sectorlist[i];
        for(unsigned long j=0; j<systemsize; j++)  element += hs[j]*szi(j, a) + J*szip1szi(j,a);
        if(dense==false)
        {
            currentTriplet = T(i,i,element);
            tripletList.push_back(currentTriplet);
        }
        if((dense==true) && (armadillobool == true))    armaH(i,i) = element;
        //if((dense==true) && (armadillobool == false))   eigenH(i,i) = element;
    } // End for-loop over i
} // End function palhuse_random_sectorHamiltonian_dense

//-------------------------------------------TOTAL HAMILTONIAN--------------------------------------------//

void Systems::palhuse_interacting_totalHamiltonian()
{
    unsigned long length = no_of_states << 3;

    currentTriplet = T(0,0,0);
    tripletList.reserve(length);     // Is this too big?

    if(dense==true)    create_armadillo_matrix();
    //else               create_dense_Eigen_matrix();
    sectorbool=false;

    unsigned long b = 0;
    for(unsigned long i=0; i<no_of_states; i++)
    {   // i is our state
        for(unsigned long j=0; j<systemsize; j++)
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
    unsigned long index = 0;
    double element = 0;
    for(unsigned long i=0; i<no_of_states; i++)
    {
        element = 0;
        // The 2-particle case is special again...
        for(unsigned long j=0; j<systemsize; j++)  element += hs[j]*szi(j, i) + J*szip1szi(j,i);
        index = no_of_states - (i+1);
        currentTriplet = T(index,index,element);
        tripletList.push_back(currentTriplet);

        if((dense==true) && (armadillobool==true))    armaH(index,index) = element;
        //if((dense==true) && (armadillobool==false))   eigenH(usual,usual) = element;
        if(dense==false)
        {
            sparseH = Eigen::SparseMatrix<double>(no_of_states, no_of_states);
        }
    } // End for-loop over i
    if(dense==false)             sparseH.setFromTriplets(tripletList.begin(), tripletList.end());
    //if((dense==true) && (armadillobool==false))    eigenH.cast<double>();
} // End function palhuse_random_sectorHamiltonian_dense

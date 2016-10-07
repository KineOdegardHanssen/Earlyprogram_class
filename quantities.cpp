#include "quantities.h"

// Okay, so this whole thing is based on only wanting to use lapack directly for the diagonalization

Quantities::Quantities()
{
    eig_all = Diagonalization();
    N = eig_all.N;
    eigvals = eig_all.eigenvalues_H;
    eigmat = eig_all.eigenmatrix_H;
    min_ev = eigvals.minCoeff();

}

//----------------------------------------BASIC FUNCTIONS-------------------------------------------------//

void Quantities::min_eigenvalue_of_sector()   // Don't need this? Can take eigvals.minCoeff() instead?
{
    // Don't need this?
    //for(int i=0; i<N; i++)
    //{
    //
    //}
}




//----------------------------------------FINDING BETA, ETC-----------------------------------------------//
void Quantities::findbeta(double eigenvalue)  // Uses Newtons method. Should I call it newtonsmethod, in case I want to test the bisection method too? The bisection method is safer, but must specify an interval.
{ // Ooops, loops?
    double betatest = 0.5; // Arbitrary small value, hopefully close to our zero point.
    double diff = 1000.0;
    int i = 0;
    double fbetan = 0.0;
    while((diff > tolerance) && (i < maxit))
    {
        fbetan = self_consistency_beta(eigenvalue, betatest);
        betaattempt -= fbetan/self_consistency_beta_derivative(eigenenergy, betatest);
        diff = abs(fbetan);
        i++;
        //if(i == maxit) cout<< "NB! Max no. of iterations exceeded in newtonsmethod. Diff="<< diff << endl;
        //prevprevdiff = prevdiff;
        //prevdiff = diff;
        //if(abs(prevprevdiff - diff) < 1e-15)  break; // To prevent the program from bouncing back and forth
    }
}

void Quantities::self_consistency_beta(double eigenvalue, double betatest)
{ // Ooops, loops?
    double Z_test = 0;               // Partition function
    double en_sum_test = 0;          // Weighted sum over energies.
    for(int i=0; i<N; i++)
    {
        Z_test += exp(betatest*(min_ev-eigvals[i])); // Or, Z_test/exp(min_ev), really. Will cancel it out of the eq.
        en_sum_test += eigvals[i]*exp(betatest*(min_ev-eigvals[i]));  // set exp(beta(min_em-eigvals[i])) instead of exp(-beta(eigvals[i]-min_em)) to save a really small amount of time... This is to be run a lot
    }
    return en_sum_test - Z_test*eigenvalue;
}

void Quantities::self_consistency_beta_derivative(double eigenvalue, double betatest)
{ // Ooops, loops?
    double Z_der_test = 0;               // Partition function
    double en_sum_der_test = 0;          // Weighted sum over energies.
    for(int i=0; i<N; i++)
    {
        Z_der_test -= eigvals[i]*exp(betatest*(min_ev-eigvals[i]));                  // Do I really need to take -= ? Yes, I think so.
        en_sum_der_test -= eigvals[i]*eigvals[i]*exp(betatest*(min_ev-eigvals[i]));
    }
    return en_sum_der_test - Z_der_test*eigenvalue;
}

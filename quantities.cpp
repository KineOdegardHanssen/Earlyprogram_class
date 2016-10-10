#include "quantities.h"

// Okay, so this whole thing is based on only wanting to use lapack directly for the diagonalization
// Armadillo is faster for some reason, though...
Quantities::Quantities(int maxit,double tolerance)
{
    this->maxit = maxit;
    this->tolerance = tolerance;
    eig_all = Diagonalization();
    N = eig_all.N;
    eigvals = eig_all.eigenvalues_H;
    eigmat = eig_all.eigenmatrix_H;
    min_ev = eigvals.minCoeff();

}

//----------------------------------------BASIC FUNCTIONS-------------------------------------------------//

void Quantities::calculateZ()
{   //Have divided by exp(-beta*min_ev). Do the same everywhere else it is needed, so that it cancels.
    Z = 0;
    for(int i=0; i<N; i++)
    {
        Z += exp(betatest*(min_ev-eigvals[i]));  // Is it wiser to point in introduce a temporary vealue for Z?
    }
}


//----------------------------------------FINDING BETA, ETC-----------------------------------------------//
void Quantities::newtonsmethod(double eigenvalue)  // Uses Newtons method. Should I call it newtonsmethod, in case I want to test the bisection method too? The bisection method is safer, but must specify an interval.
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
    beta = betaattempt;    // Is setting beta this way really the wisest?
}

void Quantities::bisectionmethod(double eigenvalue)
{
    double a = -1;
    double b = 1;
    double c = 0;
    double counter = 0;
    double fa = self_consistency_beta(eigenevalue,a);
    double fb = self_consistency_beta(eigenvalue, b);
    double diff = abs(fa);
    while(diff > tolerance && counter < maxit)
    {
        c = (a+b)/2;
        fc = self_consistency_beta(eigenvelue, c);
        if(sign(fa)==sign(fc))
        {
            a = c;
            fa = fc;
        }
        else
        {
            b = c;
            fb = fc;
        }
        counter++;
    }
    beta = c;
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

//----------------------------------------FUNCTIONS FOR THE TRACE-----------------------------------------//

 void Quantities::thermaltrace()
 {
     // Setting up the thermal matrix.
     calculateZ();
     double Zf = 1/Z;         // Since dividing is computationally expensive. But can we do something like A = A/Z ? seems a bit high-level.
     Eigen::MatrixXd A(N,N);  // Is this the correct notation?
     for(int i=0; i<N; i++)
     {
         for(int j=0; j<N; j++)    A(i,j) += Zf*exp(beta*(min_ev-eigval_s(i))*eigmat(i,i)*eigmat(j,i);   // Double check that this is correct. Think so since the eigenvectors in eigmat are stored as column vectors.s
     }
 }

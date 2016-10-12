#include "quantities.h"

// Okay, so this whole thing is based on only wanting to use lapack directly for the diagonalization
// Armadillo is faster for some reason, though...
Quantities::Quantities(int maxit, double tolerance, bool armadillobool, Systems system, Diagonalization eig_all)
{
    this->maxit = maxit;
    this->tolerance = tolerance;
    this->system = system;                         // Double check the implementation here
    this->armadillobool = armadillobool;           // Should I do something with this

    systemsize = system.systemsize;


    this->eig_all = eig_all;
    N = eig_all.N;
    if(armadillobool)   // So this is a bit complicated...
    {
        eigvals_a = eig_all.eigenvalues_armadillo;
        eigmat_a = eig_all.eigenvectors_armadillo;
        min_ev = eigvals_a(0);            // The vector is sorted so that the first eigenvalue is the smallest. But the same goes for eigenvalues_H, I guess?
    }
    else
    {
        eigvals = eig_all.eigenvalues_H;
        eigmat = eig_all.eigenmatrix_H;
        min_ev = eigvals.minCoeff();
    }


}


// I need a way to choose eigenstates to look at. Store information in a vector or matrix of some sort?


//----------------------------------------BASIC FUNCTIONS-------------------------------------------------//

/*
int Quantities::sign(double a)
{
    if(a>0.0)         return 1;
    else if(a<0.0)    return -1;
    else              return 0;  // Should have a better implementation for use in bisectionmethod, but this is a very rare case
}
*/

int Quantities::signcompare(double fa, double fc)
{
    if( (fa>0.0 && fc<0.0) || (fa<0.0 && fc>0.0) )        return    1;
    else if((fa=0.0) || (fc=0.0))                         return    0;
    else                                                  return   -1;
}


/*
int Quantities::sign_givethezero(double fa, double a, double c)  // We don't actually need this
{   // This function is only called if fa or fc is 0.
    if(fa=0.0)    return a;     // Returns a as betaattempt
    else          return c;     // Returns c as betaattempt
}
*/


void Quantities::calculateZ()
{   //Have divided by exp(-beta*min_ev). Do the same everywhere else it is needed, so that it cancels.
    Z = 0;
    for(int i=0; i<N; i++)
    {
        Z += exp(beta*(min_ev-eigvals[i]));  // Is it wiser to point in introduce a temporary vealue for Z?
    }
}


//----------------------------------------FINDING BETA, ETC-----------------------------------------------//

void Quantities::newtonsmethod(double eigenvalue)
{ // Ooops, loops?
    double betatest = 0.5; // Arbitrary small value, hopefully close to our zero point.
    double diff = 1000.0;
    int i = 0;
    double fbetan = 0.0;
    double betaattempt;
    while((diff > tolerance) && (i < maxit))
    {
        fbetan = self_consistency_beta(eigenvalue, betatest);
        betaattempt -= fbetan/self_consistency_beta_derivative(eigenvalue, betatest);
        diff = abs(fbetan);
        i++;
        //if(i == maxit) cout<< "NB! Max no. of iterations exceeded in newtonsmethod. Diff="<< diff << endl;
        //prevprevdiff = prevdiff;
        //prevdiff = diff;
        //if(abs(prevprevdiff - diff) < 1e-15)  break; // To prevent the program from bouncing back and forth
    }
    beta = betaattempt;    // Is setting beta this way really the wisest?
}

void Quantities::bisectionmethod(double eigenvalue)   // Is this a sufficiently large interval?
{
    int signcompfafb;
    double a = -1;
    double b = 1;
    double c = 0;
    double counter = 0;
    double fa = self_consistency_beta(eigenvalue,a);
    double fb = self_consistency_beta(eigenvalue, b);
    double fc;
    double diff = abs(fa);
    while(diff > tolerance && counter < maxit)
    {
        c = (a+b)/2;
        fc = self_consistency_beta(eigenvalue, c);
        signcompfafb = signcompare(fa,fc);
        if(signcompfafb==1)
        {
            a = c;
            fa = fc;
        }
        else if(signcompfafb==-1)
        {
            b = c;
            fb = fc;
        }
        else
        {
            break;  // Have encountered a zero. As we had a was in the last run of the loop and we did not encounter the zero then, c must be beta.
        }

        counter++;
    }
    beta = c;
}

double Quantities::self_consistency_beta(double eigenvalue, double betatest)
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

double Quantities::self_consistency_beta_derivative(double eigenvalue, double betatest)
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

//----------------------------------------FUNCTIONS FOR THE ETH-------------------------------------------//

void Quantities::ETH(int i)
{

}

//-------------------------------------------/USING EIGEN/------------------------------------------------//
Eigen::MatrixXd Quantities::trace_Eigen(Eigen::MatrixXd A)  // Should I just do armadillo instead?
{
    // Include an if-test
    int traceN = N >> 1; // This works for the whole matrix only. Not quite sure what will happen for sectors.
    Eigen::MatrixXd trace_matrix(traceN, traceN);  // Should I set all elements to zero?
    int index1, index2;
    for(int k=0; k<traceN; k++)
    {
        for(int l=0; l<traceN; l++)
        {
            for(int j=0; j<N; j++)
            {
                index1 = (k-1)*N+j;
                index2 = (l-1)*N+j;
                trace_matrix(k,l) += A(index1,index2);
            }
        }
    }
    return trace_matrix;
}

Eigen::MatrixXd Quantities::thermalmat_Eigen()
 {
     // Setting up the thermal matrix.
     calculateZ();
     double Zf = 1/Z;         // Since dividing is computationally expensive. But can we do something like A = A/Z ? seems a bit high-level.
     double eigmatii;
     Eigen::MatrixXd A(N,N);  // Is this the correct notation?
     for(int i=0; i<N; i++)
     {
         eigmatii = eigmat(i,i);
         for(int j=0; j<N; j++)    A(i,j) += Zf*exp(beta*(min_ev-eigvals(i)))*eigmatii*eigmat(j,i);   // Double check that this is correct. Think so since the eigenvectors in eigmat are stored as column vectors.
     }
     return A;
}

Eigen::MatrixXd Quantities::eigenstatemat_Eigen(int i)
{
    Eigen::MatrixXd B(N,N);
    for(int j=0; j<N; j++)
    {
        for(int k=0; k<N; k++)    B(i,j) += eigmat(k,i)*eigmat(j,i);
    }
    return B;
}

//----------------------------------------/USING ARMADILLO/-----------------------------------------------//
arma::mat Quantities::trace_arma(arma::mat A)
{
    // Include an if-test
    int traceN = N >> 1; // This works for the whole matrix only. Not quite sure what will happen for sectors.
    arma::mat trace_matrix(traceN, traceN);  // Should I set all elements to zero?
    int index1, index2;
    for(int k=0; k<traceN; k++)
    {
        for(int l=0; l<traceN; l++)
        {
            for(int j=0; j<N; j++)
            {
                index1 = (k-1)*N+j;
                index2 = (l-1)*N+j;
                trace_matrix(k,l) += A(index1,index2);
            }
        }
    }
    return trace_matrix;
}

arma::mat Quantities::thermalmat_arma()
{
    // Setting up the thermal matrix.
    calculateZ();
    double Zf = 1/Z;         // Since dividing is computationally expensive. But can we do something like A = A/Z ? seems a bit high-level.
    double eigmatii;
    arma::mat A(N,N);  // Is this the correct notation?
    for(int i=0; i<N; i++)
    {
        eigmatii = eigmat(i,i);  // 'Cause, getting the element each time is a bit tiring? Not important FLOPs, though...
        for(int j=0; j<N; j++)    A(i,j) += Zf*exp(beta*(min_ev-eigvals(i)))*eigmatii*eigmat(j,i);   // Double check that this is correct. Think so since the eigenvectors in eigmat are stored as column vectors.
    }
    return A;
}

arma::mat Quantities::eigenstatemat_arma(int i)
{
    arma::mat B(N,N);
    for(int j=0; j<N; j++)
    {
        for(int k=0; k<N; k++)    B(i,j) += eigmat(k,i)*eigmat(j,i);
    }
    return B;
}

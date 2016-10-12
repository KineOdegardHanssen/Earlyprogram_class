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
        cout << "Eigenvectors and eigenvalues in: " << endl;
        cout << eigvals_a << endl;
        cout << eigmat_a << endl;
    }
    else
    {
        eigvals = eig_all.eigenvalues_H;
        eigmat = eig_all.eigenmatrix_H;
        min_ev = eigvals.minCoeff();
    }
    sort_energies();
}


// I need a way to choose eigenstates to look at. Store information in a vector or matrix of some sort?


//----------------------------------------BASIC FUNCTIONS-------------------------------------------------//

void Quantities::sort_energies()
{   // Only looking at the middle energies
    li = floor(0.25*N); // Flooring just in case
    lh = ceil(0.75*N); // Is there such a function
    // Must do something about these
}

int Quantities::signcompare(double fa, double fc)
{
    if( (fa>0.0 && fc<0.0) || (fa<0.0 && fc>0.0) )        return    1;
    else if((fa=0.0) || (fc=0.0))                         return    0;
    else                                                  return   -1;
}

void Quantities::calculateZ()
{   //Have divided by exp(-beta*min_ev). Do the same everywhere else it is needed, so that it cancels.
    Z = 0;
    for(int i=0; i<N; i++)
    {
        Z += exp(beta*(min_ev-eigvals[i]));  // Is it wiser to point in introduce a temporary vealue for Z?
    }
}

void Quantities::calculateZ_arma()
{   //Have divided by exp(-beta*min_ev). Do the same everywhere else it is needed, so that it cancels.
    Z = 0;
    for(int i=0; i<N; i++)
    {
        Z += exp(beta*(min_ev-eigvals_a[i]));  // Is it wiser to point in introduce a temporary vealue for Z?
    }
}


//----------------------------------------FINDING BETA, ETC-----------------------------------------------//
//---------------------------------------------/EIGEN/----------------------------------------------------//
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
        Z_test += exp(betatest*(min_ev-eigvals(i)));                  // Or, Z_test/exp(min_ev), really. Will cancel it out of the eq.
        en_sum_test += eigvals(i)*exp(betatest*(min_ev-eigvals(i)));  // set exp(beta(min_em-eigvals[i])) instead of exp(-beta(eigvals[i]-min_em)) to save a really small amount of time... This is to be run a lot
    }
    return en_sum_test - Z_test*eigenvalue;
}

double Quantities::self_consistency_beta_derivative(double eigenvalue, double betatest)
{ // Ooops, loops?
    double Z_der_test = 0;               // Partition function
    double en_sum_der_test = 0;          // Weighted sum over energies.
    for(int i=0; i<N; i++)
    {
        Z_der_test -= eigvals(i)*exp(betatest*(min_ev-eigvals(i)));                  // Do I really need to take -= ? Yes, I think so.
        en_sum_der_test -= eigvals(i)*eigvals(i)*exp(betatest*(min_ev-eigvals(i)));
    }
    return en_sum_der_test - Z_der_test*eigenvalue;
}

//--------------------------------------------/ARMADILLO/-------------------------------------------------//

void Quantities::newtonsmethod_arma(double eigenvalue)
{ // Ooops, loops?
    double betatest = 0.5; // Arbitrary small value, hopefully close to our zero point.
    double diff = 1000.0;
    int i = 0;
    double fbetan = 0.0;
    double betaattempt;
    while((diff > tolerance) && (i < maxit))
    {
        fbetan = self_consistency_beta_a(eigenvalue, betatest);
        betaattempt -= fbetan/self_consistency_beta_derivative_a(eigenvalue, betatest);
        diff = abs(fbetan);
        i++;
        //if(i == maxit) cout<< "NB! Max no. of iterations exceeded in newtonsmethod. Diff="<< diff << endl;
        //prevprevdiff = prevdiff;
        //prevdiff = diff;
        //if(abs(prevprevdiff - diff) < 1e-15)  break; // To prevent the program from bouncing back and forth
    }
    beta = betaattempt;    // Is setting beta this way really the wisest?
}

void Quantities::bisectionmethod_arma(double eigenvalue)   // Is this a sufficiently large interval?
{
    int signcompfafb;
    double a = -1;
    double b = 1;
    double c = 0;
    double counter = 0;
    double fa = self_consistency_beta_a(eigenvalue,a);
    double fb = self_consistency_beta_a(eigenvalue, b);
    double fc;
    double diff = abs(fa);
    while(diff > tolerance && counter < maxit)
    {
        c = (a+b)/2;
        fc = self_consistency_beta_a(eigenvalue, c);
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

double Quantities::self_consistency_beta_a(double eigenvalue, double betatest)
{ // Ooops, loops?
    double Z_test = 0;               // Partition function
    double en_sum_test = 0;          // Weighted sum over energies.
    for(int i=0; i<N; i++)
    {
        Z_test += exp(betatest*(min_ev-eigvals_a(i)));                  // Or, Z_test/exp(min_ev), really. Will cancel it out of the eq.
        en_sum_test += eigvals_a(i)*exp(betatest*(min_ev-eigvals_a(i)));  // set exp(beta(min_em-eigvals[i])) instead of exp(-beta(eigvals[i]-min_em)) to save a really small amount of time... This is to be run a lot.
    }
    return en_sum_test - Z_test*eigenvalue;
}

double Quantities::self_consistency_beta_derivative_a(double eigenvalue, double betatest)
{ // Ooops, loops?
    double Z_der_test = 0;               // Partition function
    double en_sum_der_test = 0;          // Weighted sum over energies.
    for(int i=0; i<N; i++)
    {
        Z_der_test -= eigvals_a(i)*exp(betatest*(min_ev-eigvals_a(i)));                  // Do I really need to take -= ? Yes, I think so.
        en_sum_der_test -= eigvals_a(i)*eigvals_a(i)*exp(betatest*(min_ev-eigvals_a(i)));
    }
    return en_sum_der_test - Z_der_test*eigenvalue;
}

//----------------------------------------FUNCTIONS FOR THE ETH-------------------------------------------//
// I should really test these for small systems...

double Quantities::ETH(int i)
{
    if(armadillobool)    ETH_arma(i);
    else                 ETH_Eigen(i);
}

// So these are only suitable for the total system (yet). Should generalize them
double Quantities::ETH_arma(int i)         // Or should I return a list of doubles, so I have more flexibility regarding which quantity to use
{
    int j = 0;
    //newtonsmethod_arma(eigvals_a(i));   // Or should I call the bisection method? // newtonsmethod works, but is incredibly slow...
    //cout << "newtonsmethod_arma run" << endl;
    beta = 0; // Incorrect, but want to focus on removing any errors for now.
    arma::mat thm = thermalmat_arma();
    arma::mat esm = eigenstatemat_arma(i);
    cout << "The density matrix is:" << endl;
    cout << esm << endl << endl;
    cout << "The thermal matrix is: " << endl;
    cout << thm << endl << endl;

    // Trace procedures: Should we trace over all spins except our state?
    while(j<(systemsize-1)) // Want to trace over all particles except one
    {
        thm = trace_arma(thm);  // Should I declare it again, or just let it stand like this?
        esm = trace_arma(esm);
        j++;
    }
    cout << "Reduced thermal matrix: " << endl;
    cout << thm << endl << endl;
    cout << "Reduced density matrix: " << endl;
    cout << esm << endl << endl;
    arma::mat diff_mat = thm - esm;
    cout << "matrix diff_mat created" << endl;
    double diff_elem = diff_mat(0,0);   // Difference between the first elements
    //double diff_norm_frob = norm(diff_mat, "fro"); // Frobenius norm
    //double diff_norm_maxsumrow = norm(diff_mat, "inf");
    //double diff_norm_maxsumcol = norm(diff_mat, 1);

    // Do something with some of these
    return diff_elem;
}

double Quantities::ETH_Eigen(int i)
{
    int j=0;
    newtonsmethod(eigvals(i));
    Eigen::MatrixXd thm = thermalmat_Eigen();
    Eigen::MatrixXd esm = eigenstatemat_Eigen(i);
    // Trace procedures: Should we trace over all spins except our state?
    while(j<(systemsize-1)) // Want to trace over all particles except one
    {
        thm = trace_Eigen(thm);  // Should I declare it again, or just let it stand like this?
        esm = trace_Eigen(esm);
        j++;
    }
    Eigen::MatrixXd diff_mat = thm - esm;
    double diff_elem = diff_mat(0,0);   // Difference between the first elements
    //double diff_norm_frob = Eigen::MatrixBase<Eigen::MatrixXd>::norm(diff_mat);   // The Frobenius norm

    return diff_elem;
}

//-------------------------------------------/USING EIGEN/------------------------------------------------//
/*
Eigen::MatrixXd Quantities::trace_Eigen(Eigen::MatrixXd A)  // Should I just do armadillo instead?
{
    // Include an if-test
    int traceN = N >> 1; // This works for the whole matrix only. Not quite sure what will happen for sectors.
    Eigen::MatrixXd trace_matrix(traceN, traceN);  // Should I set all elements to zero?
    for(int k=0; k<traceN; k++)
    {
        for(int l=0; l<traceN; l++)    trace_matrix(k,l) = 0.0;
    }
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
*/

Eigen::MatrixXd Quantities::trace_Eigen(Eigen::MatrixXd A)  // Should I just do armadillo instead?
{
    // Include an if-test
    int n = N >> 1; // This works for the whole matrix only. Not quite sure what will happen for sectors.
    Eigen::MatrixXd trace_matrix(n, n);  // Should I set all elements to zero?
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)        trace_matrix(i,j) = A(2*i,2*j) + A(2*i+1,2*j+1);    // This should work, though, I haven't looked it up.
    }
    return trace_matrix;
}



Eigen::MatrixXd Quantities::thermalmat_Eigen()
 {
     // Setting up the thermal matrix.
     calculateZ();
     double Zf = 1/Z;         // Since dividing is computationally expensive. But can we do something like A = A/Z ? seems a bit high-level.
     double eigmatik;
     Eigen::MatrixXd A(N,N);  // Is this the correct notation?
     for(int k=0; k<N; k++)
     {
         for(int i=0; i<N; i++)
         {
             eigmatik = eigmat_a(i,k);  // 'Cause, getting the element each time is a bit tiring? Not important FLOPs, though...
             cout << "i: " << i <<"eigmatii: " << eigmatik << endl;
             for(int j=0; j<N; j++)    A(i,j) += Zf*exp(beta*(min_ev-eigvals_a(k)))*eigmatik*eigmat_a(j,k);   // Double check that this is correct. Think so since the eigenvectors in eigmat are stored as column vectors.
         }
     }
     return A;
}

Eigen::MatrixXd Quantities::eigenstatemat_Eigen(int i)
{

    Eigen::MatrixXd B(N,N);   
    for(int j=0; j<N; j++)
    {
        for(int k=0; k<N; k++)    B(k,j) = eigmat(k,i)*eigmat(j,i);
    }
    return B;
}

//----------------------------------------/USING ARMADILLO/-----------------------------------------------//
/*
arma::mat Quantities::trace_arma(arma::mat A)
{
    // Include an if-test
    int traceN = N >> 1; // This works for the whole matrix only. Not quite sure what will happen for sectors.
    int n = N-1; // because the counting in Maziero was different.
    arma::mat trace_matrix(traceN, traceN);  // Should I set all elements to zero?
    int index1, index2;
    for(int k=0; k<traceN; k++)   // If neccessary
    {
        for(int l=0; l<traceN; l++)    trace_matrix(k,l) = 0.0;
    }

    for(int k=0; k<traceN; k++)
    {
        for(int l=0; l<traceN; l++)
        {
            for(int j=0; j<N; j++)
            {
                index1 = abs((k-1)*n+j);  // Well, so maybe this isn't a good way to take the trace...
                index2 = abs((l-1)*n+j);
                //cout << "Index1: " << index1 << "; Index2: " << index2 << endl;
                trace_matrix(k,l) += A(index1,index2);
            }
        }
    }
    return trace_matrix;
}
*/

arma::mat Quantities::trace_arma(arma::mat A)  // Should I just do armadillo instead?
{
    // Include an if-test
    int n = N >> 1; // This works for the whole matrix only. Not quite sure what will happen for sectors.
    arma::mat trace_matrix(n, n);  // Should I set all elements to zero?
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)        trace_matrix(i,j) = A(2*i,2*j) + A(2*i+1,2*j+1);    // This should work, though, I haven't looked it up.
    }
    return trace_matrix;
}

arma::mat Quantities::thermalmat_arma()
{
    // Setting up the thermal matrix.
    calculateZ_arma();
    double Zf = 1/Z;         // Since dividing is computationally expensive. But can we do something like A = A/Z ? seems a bit high-level.
    double eigmatki;
    arma::mat A(N,N);  // Is this the correct notation?
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)    A(i,j) = 0.0;   // See if this helps
    }
    cout << "And now, I have set A" << endl;
    for(int k=0; k<N; k++)
    {
        for(int i=0; i<N; i++)
        {
            eigmatki = eigmat_a(i,k);  // 'Cause, getting the element each time is a bit tiring? Not important FLOPs, though...
            cout << "i: " << i <<"eigmatii: " << eigmatki << endl;
            for(int j=0; j<N; j++)    A(i,j) += Zf*exp(beta*(min_ev-eigvals_a(k)))*eigmatki*eigmat_a(j,k);   // Double check that this is correct. Think so since the eigenvectors in eigmat are stored as column vectors.
        }
    }

    return A;
}

arma::mat Quantities::eigenstatemat_arma(int i)
{
    arma::mat B(N,N);
    for(int j=0; j<N; j++)
    {
        for(int k=0; k<N; k++)    B(k,j) = eigmat_a(k,i)*eigmat_a(j,i);
    }
    return B;
}

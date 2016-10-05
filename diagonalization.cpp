#include "diagonalization.h"
#include<lapacke.h>

Diagonalization::Diagonalization()
{
}

Diagonalization::Diagonalization(Systems given)
{
    this->given = given;
}

void Diagonalization::lapack_directly()
{
    /*
    if(given.armadillobool)    arma::mat A = given.armaH;        // What kind of matrix should A be?
    else                       Eigen::MatrixXd A = given.eigenH;
    if(given.sectorbool)       N = given.number_of_hits;
    else                       N = given.no_of_states;
    */
    double start_time_lapack_directly = clock();

    // Tried to go complex. Didn't work.
    /*
    arma::mat complexpart = arma::zeros(N,N);       // The complex part is zero.
    arma::mat C  = arma::mat(N,N);
    C(0,1) = 1.1;
    C(1,1) = 1.4;
    arma::cx_mat A = arma::cx_mat(B, C);  // Not correct, just going to see if it will run
    */

    // I just use lapacke.h
    /*
    #define LAPACKE

    #ifdef LAPACKE   // Can I put this somewhere that is a bit more sensible?
    {
    */

    #define LAPACK_DISABLE_NAN_CHECK
    #define LAPACK_COMPLEX_CUSTOM
    #define lapack_complex_float std::complex<float>
    #define lapack_complex_double std::complex<double>
    /*
    #include<lapacke.h>

    }
    #else
    {
    //namespace mylapacke{
    #define LAPACK_DISABLE_NAN_CHECK
    #define LAPACK_COMPLEX_CUSTOM
    #define lapack_complex_float std::complex<float>
    #define lapack_complex_double std::complex<double>
    #include "mkl_lapacke.h"
    }
    #endif
    */

    // Maybe convert elements to what they should be? Use a cast function
    JOBS = 'V';
    UPLO = 'U';
    LDA  = 'N';
    int LWORK = 2*N-1;
    //int INFO = 0;
    std::vector<double> WORK  = std::vector<double> (LWORK); // What filetype is this?
    std::vector<double> RWORK = std::vector<double>(3*N-2);
    std::vector<double> W     = std::vector<double>(N);

    /*
    #ifdef LAPACKE
        lapack::LAPACKE_zheev_work(LAPACK_COL_MAJOR, JOBS, UPLO,
                                      N, &A[0],LDA, &W[0],&WORK[0],LWORK,&RWORK[0]);
    #else
    LAPACKE_zheev(LAPACK_COL_MAJOR, JOBS, UPLO,N, reinterpret_cast<lapack_complex_double*>(&A[0]),LDA, &W[0]);
    #endif
    */

    //complex<double> *A1 = given.complexH.data();

    double *A1 = given.eigenH.data();

    // dsyev (JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)

    //LAPACKE_zheev_work(LAPACK_COL_MAJOR, JOBS, UPLO, N, A1, LDA, &W[0],&WORK[0],LWORK,&RWORK[0]);

    LAPACKE_dsyev_work(LAPACK_COL_MAJOR, JOBS, UPLO, N, A1, LDA, &W[0],&WORK[0],LWORK);

    double end_time_lapack_directly = clock();

    lapacktime = (end_time_lapack_directly - start_time_lapack_directly)/CLOCKS_PER_SEC;


    //cout << &A << endl;

    /*
    for(int i=0; i<N; i++)
    {
        cout << "Eigenvalue " << (i+1) << " is: " << W[i] << endl;
    }
    */


}


void Diagonalization::using_armadillo()
{
    double start_time_arma = clock();        // If I want to take the time
    arma_n = 0;
    if(given.sectorbool==true)          arma_n = given.number_of_hits;
    else                                arma_n = given.no_of_states;
    eigenvalues_armadillo = arma::vec(arma_n);
    eigenvectors_armadillo = arma::mat(arma_n,arma_n);

    //string method = "std";
    arma::eig_sym(eigenvalues_armadillo,eigenvectors_armadillo, given.armaH);
    double end_time_arma = clock();
    armatime = (end_time_arma - start_time_arma)/CLOCKS_PER_SEC;
}

void Diagonalization::print_using_armadillo()
{
    for(unsigned long i= 0; i<arma_n; i++)
    {
        cout << "Eigenvalue = " << eigenvalues_armadillo(i) << endl;
        cout << "Its corresponding eigenvector:" << endl;
        for(unsigned long j=0; j<arma_n; j++)       cout << eigenvectors_armadillo(j,i) << " ";
        cout << endl;
    }   // End for-loop over i
}

void Diagonalization::print_sparse_Eigen()
{
    // See if i can find something a bit more clever...
    cout << given.sparseH << endl;
}


// Here is the danger zone!

void Diagonalization::using_dense_eigen()
{
    Eigen::EigenSolver<Eigen::MatrixXd> es(given.eigenH);  // Why won't this work?

    //Eigen::VectorXd eigenvalues_H = es.eigenvalues();
    //Eigen::MatrixXd eigenmatrix_H = es.eigenvectors();

}


void Diagonalization::print_dense_using_eigen()
{
    Eigen::EigenSolver<Eigen::MatrixXd> es(given.eigenH.cast<double>());
    cout << "The eigenvalues of H are:" << endl << es.eigenvalues() << endl;
    cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
}




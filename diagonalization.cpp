#include "diagonalization.h"
#include<lapacke.h>
//#include <stdlib.h>
//#include <stdio.h>

Diagonalization::Diagonalization()
{
}

Diagonalization::Diagonalization(Systems given)
{
    this->given = given;
}

void Diagonalization::lapack_directly()
{
    cout << "In diagonalization" << endl;
    const bool TRACE = true;
    double start_time_lapack_directly = clock();

    cout << "Did I assign N?" << endl;
    if(given.sectorbool==true)          N = given.number_of_hits;
    else                                N = given.no_of_states;
    cout << "Yes, I assigned N. Assigning the other input parameters:" << endl;

    // note, to understand this part take a look in the MAN pages, at section of parameters.
    //double wkopt; // Some optimalization parameter...

    // My thing

    JOBS = 'V';
    if(TRACE)          cout << "JOBS OK," << endl;
    UPLO = 'L';
    if(TRACE)          cout << "UPLO OK," << endl;
    LDA  = N;
    if(TRACE)          cout << "LDA OK," << endl;
    int LWORK = 3*N-1;
    if(TRACE)          cout << "LWORK OK," << endl;
    int INFO;
    if(TRACE)          cout << "INFO declared," << endl;

    std::vector<double> WORK  = std::vector<double> (LWORK); // What filetype is this?
    if(TRACE)          cout << "WORK OK," << endl;
    //std::vector<double> RWORK = std::vector<double>(3*N-2);
    std::vector<double> W     = std::vector<double>(N);
    if(TRACE)          cout << "W OK," << endl;


    // Intel's thing:
    //int n = N, lda = LDA, info, lwork;
    //double wkopt;
    /* Local arrays */
    //double w[N];

    // My thing:
    if(TRACE)    cout << "Before declaration double *A1" << endl;
    double *A = given.eigenH.data(); // Must I rewrite this?
    if(TRACE)    cout << "After declaration double *A1" << endl;

    // end of declarations


    // Intel's thing:
    /*
    lwork = -1;
    LAPACK_dsyev( "Vectors".c_str(), "Upper".c_str(), &n, A, &lda, w, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    std::vector<double> work(lwork);
    /* Solve eigenproblem */
    //LAPACK_dsyev( "Vectors".c_str(), "Upper".c_str(), &n, A, &lda, w, &work[0], &lwork, &info );
    /* Check for convergence
    if( info > 0 )
    {
        printf( "The algorithm failed to compute eigenvalues.\n" );
        exit( 1 );
    }
    double end_time_lapack_directly = clock();
    lapacktime = (end_time_lapack_directly - start_time_lapack_directly)/CLOCKS_PER_SEC;

    cout << "Eigenvectors: " << endl;
    for(int i = 0; i<N; i++ )
    {
        for(int j = 0; j < n; j++ )     cout << A[i+j*lda] << " ";
        cout << endl;
    }

    cout << "Eigenvalues: " << endl;
    for(int i = 0; i<N; i++ )        cout << w[i+lda] << " ";
    cout << endl;
    */

    // My thing:
 

    cout << "compute the LU factorization..." << endl << endl;

    for(int i=0; i<N*N; i++)
    {
        cout << A[i] << " ";
    }
    cout << endl;

    cout << "W" << endl;
    for(int i=0; i<N; i++)
    {
        cout << W[i] << " ";
    }
    cout << endl;

    cout << "JOBS: " << JOBS << "; UPLO: " << UPLO << endl;

    LAPACK_dsyev(&JOBS, &UPLO, &N, A, &LDA, &W[0],&WORK[0],&LWORK, &INFO);
    //LWORK = (int)wkopt;
    //WORK = (double*)malloc( LWORK*sizeof(double) );
    //LAPACK_dsyev(&JOBS, &UPLO, &N, A, &LDA, &W[0],&WORK[0],&LWORK, &INFO);
    double end_time_lapack_directly = clock();
    if(TRACE)    cout << "Have run LAPACK_dsyev" << endl;

    cout << "JOBS: " << JOBS << "; UPLO: " << UPLO << endl;

    for(int i=0; i<N*N; i++)
    {
        cout << A[i] << " ";
    }
    cout << endl;

    cout << "W" << endl;
    for(int i=0; i<N; i++)
    {
        cout << W[i] << " ";
    }
    cout << endl;


    cout << "Info has the value: " << INFO << endl;

    cout << "Before using Eigen::Map: " << endl;

    // checks INFO, if INFO != 0 something goes wrong, for more information see the MAN page of dgetrf.
    if(INFO!=0)
    {
        cout << "an error occured : "<< INFO << endl << endl;
    }else{
        cout << "system solved..."<< endl << endl;
        lapacktime = (end_time_lapack_directly - start_time_lapack_directly)/CLOCKS_PER_SEC;
        if(INFO!=0)
        {
            // checks INFO, if INFO != 0 something goes wrong, for more information see the MAN page of dgetrs.
            cout << "an error occured : "<< INFO << endl << endl;
            lapacktime = 1e16;
        }else{
            cout << "Eigenvectors: ";
            for (int i=0;i<N;i++)
            {
                cout << "[" << endl;
                for(int j=0; j<N; j++)
                {
                    cout << A[i+j*LDA] << " ";
                }
                cout << "]" << endl;
            }
            cout << endl;

            cout << "Corresponding eigenvalues: " << endl;
            for(int k=0; k<N; k++)    cout << W[k] << ", ";
            cout << endl;
        }
    }

    cout << "After using Eigen::Map: " << endl;
    Eigen::MatrixXd B = Eigen::Map<Eigen::Matrix<double, N, N> > (A);
    //Eigen::MatrixXd B = Eigen::Map<Eigen::Matrix<double, N, N> > (A);
    //Eigen::MatrixXd B = Eigen::Map<Eigen::MatrixXd(N,N) > (A); // The dimensions ...

    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)    cout << B(i,j) << " ";
        cout << endl;
    }




    if(TRACE)    cout << "Everything should be running just fine" << endl;

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
    for(int i= 0; i<arma_n; i++)
    {
        cout << "Eigenvalue = " << eigenvalues_armadillo(i) << endl;
        cout << "Its corresponding eigenvector:" << endl;
        for(int j=0; j<arma_n; j++)       cout << eigenvectors_armadillo(j,i) << " ";
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




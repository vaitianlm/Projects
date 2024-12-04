#include <cmath>
#include "utils.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <complex>
using namespace std::complex_literals;
using namespace std;

// Contructor
Double_slit::Double_slit(int points, double timestep, double spacestep, arma::cx_vec a, arma::cx_vec b)
{
    M = points;
    dt = timestep;
    h = spacestep;
    r = 1i*dt/(2*pow(h, 2));
    //A = CN_matA(r, a);
    //B = CN_matB(r, b);
}

// Translates two indices i,j of (M-2 x M-2) matrix to a single index k
int Double_slit::vec_ind(int i, int j, int M)
{
    return i + j*(M - 2);
}

// generates vector containing diagonal elements of matrix A
arma::cx_vec Double_slit::a_vec(const arma::sp_cx_mat& V, arma::cx_double r, arma::cx_double dt)
{
    int N = V.n_elem;
    int M = V.n_cols + 2;

    arma::cx_double a0 = 1.0 + 4.0 * r;
    arma::cx_vec a = arma::cx_vec(N, arma::fill::ones)*a0;

    arma::sp_cx_mat::const_iterator it = V.begin();
    arma::sp_cx_mat::const_iterator it_end = V.end();

    int i;
    int j;
    int k;

    for(; it != it_end; ++it)
    {
        i = it.row();
        j = it.col();
        k = vec_ind(i, j, N+2);
        
        a(k) += 1i*dt* (*it) /2.0;
    }

    return a;
}

// generates vector containing diagonal elements of matrix B
arma::cx_vec Double_slit::b_vec(const arma::sp_cx_mat& V, arma::cx_double r, arma::cx_double dt)
{
    int N = V.n_elem;
    int M = V.n_cols + 2;

    arma::cx_double a0 = 1.0 - 4.0 * r;
    arma::cx_vec a = arma::cx_vec(N, arma::fill::ones)*a0;

    arma::sp_cx_mat::const_iterator it = V.begin();
    arma::sp_cx_mat::const_iterator it_end = V.end();

    int i;
    int j;
    int k;

    for(; it != it_end; ++it)
    {
        i = it.row();
        j = it.col();
        k = vec_ind(i, j, N+2);
        
        a(k) -= 1i*dt* (*it) /2.0;
    }

    return a;
}

// Helper function for CN_matA and CN_matB
arma::sp_cx_mat Double_slit::tridiag(arma::cx_double r, const arma::cx_vec& di)
{
    int N = di.n_elem;

    // Creating diagonal matrix filled with the vector di
    arma::sp_cx_mat Mi(N, N);
    Mi.diag(0) = di;

    // Vector for filling the sub- and superdiagonal
    arma::cx_vec off_diag(N - 1, arma::fill::value(r));

    // Subdiagonal
    Mi.diag(-1) = off_diag;

    // Superdiagonal
    Mi.diag(1) = off_diag;
  
  return Mi;
}

// Constructs matrices A and B for CN scheme (OPPGAVE 2)
arma::sp_cx_mat Double_slit::CN_matAB(arma::cx_double r, arma::cx_vec d)
{
    // Creating matrix A/B and filling the diagonal
    int N = d.n_elem;
    arma::sp_cx_mat Mat(N, N);
    Mat.diag(0) = d;

    // Vector for filling the sub- and superdiagonal
    arma::cx_vec off_diag(N - 1, arma::fill::value(r));
    arma::uvec indices_to_zero = arma::regspace<arma::uvec>(M - 3, M - 2, N - 2);
    off_diag.elem(indices_to_zero).zeros();

    // Subdiagonal
    Mat.diag(-1) = off_diag;

    // Superdiagonal
    Mat.diag(1) = off_diag;

    // Filling the last diagonals
    arma::cx_vec other_diag(N - (M - 2), arma::fill::value(r));
    Mat.diag(M - 2) = other_diag;
    Mat.diag(-(M - 2)) = other_diag;

    return Mat;
}

#include <armadillo>
#include <vector>
#include <string>

// A function that prints the structure of a sparse matrix to screen.
void print_sp_matrix_structure(const arma::sp_cx_mat& A)
{
    using namespace std;
    using namespace arma;

    // Declare a C-style 2D array of strings.
    string S[A.n_rows][A.n_cols];  

    // Initialise all the strings to " ".
    for (int i =0; i < A.n_rows; i++)
    {
        for (int j = 0; j < A.n_cols; j++)
        {
            S[i][j] = " ";
        }
    }

    // Next, we want to set the string to a dot at each non-zero element.
    // To do this we use the special loop iterator from the sp_cx_mat class
    // to help us loop over only the non-zero matrix elements.
    sp_cx_mat::const_iterator it     = A.begin();
    sp_cx_mat::const_iterator it_end = A.end();

    int nnz = 0;
    for(it; it != it_end; ++it)
    {
        S[it.row()][it.col()] = "•";
        nnz++;
    }

    // Finally, print the matrix to screen.
    cout << endl;
    for (int i =0; i < A.n_rows; i++)
    {
        cout << "| ";
        for (int j = 0; j < A.n_cols; j++)
        {
            cout << S[i][j] << " ";
        }
        cout <<  "|\n";
    }

    cout << endl;
    cout << "matrix size: " << A.n_rows << "x" << A.n_cols << endl;
    cout << "non-zero elements: " << nnz << endl ;
    cout << endl;
}




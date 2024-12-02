#include <cmath>
#include "utils.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

// Contructor
Double_slit::Double_slit(double timestep, double spacestep, arma::cx_vec a, arma::cx_vec b)
{
    dt = timestep;
    h = spacestep;
    r = 1i*dt/(2*pow(h, 2));
    //A = CN_matA(r, a);
    //B = CN_matB(r, b);
}

// Translates two indices i,j to a single index k
int Double_slit::vec_ind(int i, int j, int M)
{
    return i + j*(M - 2);
}


// // Generates vector containing diagonal elements of matrix A
// arma::cx_vec Double_slit::a_vec(arma::cx_mat V, arma::cx_double r, arma::cx_double dt)
// {
//     1;
// }

// // Generates vector containing diagonal elements of matrix B
// arma::cx_vec Double_slit::b_vec(arma::cx_mat V, arma::cx_double r, arma::cx_double dt)
// {
//     1;
// }

// Helper function for CN_matA and CN_matB
arma::cx_mat Double_slit::tridiag(arma::cx_double r, arma::cx_vec di)
{
    int N = di.n_elem;
    arma::cx_mat Mi = arma::cx_mat(N, N, arma::fill::zeros);

    // Filling the first row
    Mi(0,0) = di(0);
    Mi(0,1) = r;

    // Loop that fills rows 2 to n-1 (row indices 1 to n-2)
    for (int i = 1; i <= N-2; i++)
    {
        Mi(i, i-1) = r;
        Mi(i, i) = di(i);
        Mi(i, i+1) = r;
    }

    // Fill last row (row index n-1)
    Mi(N-1, N-2) = r;
    Mi(N-1, N-1) = di(N-1);
  return Mi;
}

// // Constructs matrix A for CN scheme (OPPGAVE 2)
// arma::cx_mat Double_slit::CN_matA(int M, arma::cx_double r, arma::cx_vec a)
// {
//     1;
// }

// // Constructs matrix B for CN scheme (OPPGAVE 2)
// arma::cx_mat Double_slit::CN_matB(int M, arma::cx_double r, arma::cx_vec b)
// {
//     1;
// }





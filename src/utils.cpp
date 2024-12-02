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
    A = CN_matA(r, a);
    B = CN_matB(r, b);
}

// Translates two indices i,j to single indice k  (OPPGAVE 2)
int Double_slit::vec_ind(int i, int j)
{
    1;
}

// generates vector containing diagonal elements of matrix A
arma::cx_vec Double_slit::a_vec(arma::cx_mat V, arma::cx_double r, arma::cx_double dt)
{
    1;
}

// generates vector containing diagonal elements of matrix B
arma::cx_vec Double_slit::b_vec(arma::cx_mat V, arma::cx_double r, arma::cx_double dt)
{
    1;
}

// Helper function for CN_matA and CN_matB
arma::cx_mat Double_slit::tridiag(arma::cx_double r, arma::cx_vec diag)
{
    1;
}

// Constructs matrix A for CN scheme (OPPGAVE 2)
arma::cx_mat Double_slit::CN_matA(arma::cx_double r, arma::cx_vec a)
{
    1;
}

// Constructs matrix B for CN scheme (OPPGAVE 2)
arma::cx_mat Double_slit::CN_matB(arma::cx_double r, arma::cx_vec b)
{
    1;
}

// COPY PASTE FRA PROSJEKT 2
// This function creates a matrix tridiag(a,d,e) of size n*n, 
// from scalar input a, d, and e. Taken from the code snippets provided in the problem text.
arma::mat create_tridiagonal(int N, double a, double d, double e)
{
    // Start from d* identity matrix
    arma::mat A = arma::mat(N, N, arma::fill::eye) * d;

    // Fill the first row (row index 0), e.g.
    A(0,1) = e;

    // Loop that fills rows 2 to n-1 (row indices 1 to n-2)
    for (int i = 1; i <= N-2; i++)
    {
        A(i, i-1) = a;
        A(i, i+1) = e;
    }

    // Fill last row (row index n-1)
    A(N-1, N-2) = a;
  return A;
}




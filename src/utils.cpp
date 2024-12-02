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




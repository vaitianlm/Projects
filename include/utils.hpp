// Include guard
#ifndef __utils_hpp__  
#define __utils_hpp__

#include <armadillo>;

class Double_slit
{
public:
    
    double dt;          // Stepsize time
    double h;           // Stepsize space
    arma::cx_double r;  // i*dt/2h^2
    arma::cx_mat A;     // CN scheme matrix A
    arma::cx_mat B;     // CN scheme matrix B

    // Contructor
    Double_slit(double dt, double h, arma::cx_vec a, arma::cx_vec b);

    // Translates two indices i,j to single indice k  (O)
    int vec_ind(int i, int j);

    // Constructs matrix A for CN scheme
    arma::cx_mat CN_matA(arma::cx_double r, arma::cx_vec a);

    // Constructs matrix B for CN scheme
    arma::cx_mat CN_matB(arma::cx_double r, arma::cx_vec b);
}




#endif  // end of include guard __utils_hpp__
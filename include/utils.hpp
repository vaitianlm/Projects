// Include guard
#ifndef __utils_hpp__  
#define __utils_hpp__

#include <armadillo>

class Double_slit
{
public:
    int M;              // Points along the spacial axes
    double dt;          // Stepsize time
    double h;           // Stepsize space
    arma::cx_double r;  // i*dt/2h^2
    // arma::sp_cx_mat A;     // CN scheme matrix A
    // arma::sp_cx_mat B;     // CN scheme matrix B

    // Contructor
    Double_slit(int points, double dt, double h, arma::cx_vec a, arma::cx_vec b);

    // Translates two indices i,j to single indice k  (O)
    int vec_ind(int i, int j, int M);

    // Generates vector containing diagonal elements of matrix A
    arma::cx_vec a_vec(const arma::sp_cx_mat& V, arma::cx_double r, arma::cx_double dt);

    // generates vector containing diagonal elements of matrix B
    arma::cx_vec b_vec(const arma::sp_cx_mat& V, arma::cx_double r, arma::cx_double dt);

    // Helper function for CN_matA and CN_matB
    arma::sp_cx_mat tridiag(arma::cx_double r, const arma::cx_vec& di);

    // Constructs matrices A and B for CN scheme
    arma::sp_cx_mat CN_matAB( arma::cx_double r, arma::cx_vec d);

};

void print_sp_matrix_structure(const arma::sp_cx_mat& A);

#endif  // end of include guard __utils_hpp__
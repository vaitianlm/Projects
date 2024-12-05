// Include guard
#ifndef __double_slit_hpp__  
#define __double_slit_hpp__

#include <armadillo>

class Double_slit
{
public:
    double dt;          // Stepsize time
    double h;           // Stepsize space
    int M;              // Points along the spacial axes
    arma::sp_cx_mat V;  // Potential matrix
    arma::cx_vec u;     // Wave function
    arma::cx_double r;  // i*dt/2h^2
    arma::sp_cx_mat A;     // CN scheme matrix A
    arma::sp_cx_mat B;     // CN scheme matrix B

    // Contructor
    Double_slit(int points, double dt, double h, double xc, double yc, double px, double py, double sigmax, double sigmay);

    // Translates two indices i,j to single indice k  
    int vec_ind(int i, int j, int M);

    // Generates vector containing diagonal elements of matrix A or B
    // Input 1 for A and -1 for B
    arma::cx_vec ab_vec(double pm);

    // Constructs matrices A or B for CN scheme
    // Input -r, ab_vec(1) for A and r, ab_vec(-1) for B
    arma::sp_cx_mat CN_matAB(arma::cx_double r, arma::cx_vec d);

    // Evolves wavefunction one timestep using Crank-Nicolson scheme
    void CN_step();

    // Creates the wavefunction u(x,y,t=0)
    arma::cx_vec u_init(double xc,double yc, double px, double py, double sigmax, double sigmay);

    // Creates potential matrix V containing potantial at every point x,y
    void double_slit_init();
};

void print_sp_matrix_structure(const arma::sp_cx_mat& A);

#endif  // end of include guard __double_slit_hpp__
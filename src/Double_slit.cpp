#include <cmath>
#include "Double_slit.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <complex>
using namespace std::complex_literals;
using namespace std;

// Contructor
Double_slit::Double_slit(int points, double timestep, double spacestep, double x_c, double y_c, double p_x, double p_y, double sigma_x, double sigma_y)
{
    dt = timestep;
    h = spacestep;
    M = round(1/h);
    r = 1i*dt/(2*pow(h, 2));
    A = CN_matAB(-r, ab_vec(1));
    B = CN_matAB(r, ab_vec(-1));

    V = arma::sp_cx_mat(M-2, M-2);
    double_slit_init();

    u = u_init(x_c, y_c, p_x, p_y, sigma_x, sigma_y);
}

// Translates two indices i,j of (M-2 x M-2) matrix to a single index k
int Double_slit::vec_ind(int i, int j, int M)
{
    return i + j*(M - 2);
}

// Generates vector containing diagonal elements of matrix A.
// Takes 1 as argument for A and -1 as argument for B.
arma::cx_vec Double_slit::ab_vec(double pm)
{
    int N = V.n_elem;

    arma::cx_double a0 = 1.0 + pm * 4.0 * r;
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
        
        a(k) += pm * 1i*dt* (*it) /2.0;
    }

    return a;
}



// Outputs matrix A or B for CN scheme.
// Input -r, ab_vec(1) for A and r, ab_vec(-1) for B
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

// Evolves wavefunction one timestep using Crank-Nicolson scheme
void Double_slit::CN_step()
{
    arma::cx_vec b = B*u;
    u = arma::spsolve(A, b);
}

arma::cx_vec Double_slit::u_init(double xc, double yc, double px, double py, double sigmax, double sigmay)
{
    int N_u = (M-2) * (M-2);
    arma::cx_vec u(N_u);

    double x, y, dx, dy;
    for (int i = 0; i <= M-3; i++)
    {
        x = i*h;
        dx = x - xc;
        for (int j = 0; j <= M-3; j++)
        {
            y = j*h;
            dy = y - yc;
            u(vec_ind(i, j, M)) = exp( -dx*dx/(2.0*sigmax*sigmax) -dy*dy/(2.0*sigmay*sigmay) + 1i * px * x + 1i * py * y);
        }
    }

    u = arma::normalise(u);
    return u;
}

// Creates potential matrix V containing potential at every point x,y
void Double_slit::double_slit_init()
{
    // Potential strength
    double v0 = 1e10;

    // Barrier measurements
    double x_pos_barrier = 0.5;
    double barrier_half_width = 0.01;
    double y_middle = 0.05;
    double midbar_half_width = 0.025;
    double aperture = 0.05;
    
    // First and last x-indice of barrier
    int barrier_xf = floor((x_pos_barrier - barrier_half_width)/h);
    int barrier_xl = ceil((x_pos_barrier + barrier_half_width)/h);

    // y-indices of barrier ends
    int y_top_end = floor((y_middle-midbar_half_width-aperture)/h);
    int y_midbar_top = floor((y_middle-midbar_half_width)/h);
    int y_midbar_bot = ceil((y_middle+midbar_half_width)/h);
    int y_bot_end = ceil((y_middle+midbar_half_width+aperture)/h);

    // Filling V to make potential wall with slits
    for (int j = 0; j <= y_top_end; j++)
    {
        for (int i = barrier_xf; i <= barrier_xl; i++)
        {
            V(i, j) = v0;
        }
    }

    for (int j = y_midbar_top; j <= y_midbar_bot; j++)
    {
        for (int i = barrier_xf; i <= barrier_xl; i++)
        {
            V(i, j) = v0;
        }
    }

    int N = M - 2;
    for (int j = y_bot_end; j <= N-1; j++)
    {
        for (int i = barrier_xf; i <= barrier_xl; i++)
        {
            V(i, j) = v0;
        }
    }
}






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




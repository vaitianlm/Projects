#include <cmath>
#include "Quantum_box.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <complex>
using namespace std::complex_literals;
using namespace std;

// Contructor
Quantum_box::Quantum_box(double T_in, double dt_in, double h_in, int slits, double xc_in, 
                            double yc_in, double px_in, double py_in, double sigmax_in, double sigmay_in)
{
    // Setting parameters
    T = T_in;
    dt = dt_in;
    h = h_in;
    M = round(1/h)+1;

    // Initialising wavefunction
    u = u_init(xc_in, yc_in, px_in, py_in, sigmax_in, sigmay_in);

    // Initialising potential
    slits_init(slits);
    
    // Initialising A and B matrices
    r = 1i*dt/(2*h*h);
    A = CN_matAB(-r, ab_vec(1));
    B = CN_matAB(r, ab_vec(-1));

    // 3D-matrix to store wavefunction at every timestep
    int t_steps = ceil(T/dt)+1;
    S.set_size(M-2, M-2, t_steps);
    save_wf(0);

    norm_dev.set_size(t_steps);
    deviation(0);
}

// Translates two indices i,j of (M-2 x M-2) matrix to a single index k
int Quantum_box::vec_ind(int i, int j)
{
    return i + j*(M - 2);
}

// Generates vector containing diagonal elements of matrix A.
// Takes 1 as argument for A and -1 as argument for B.
arma::cx_vec Quantum_box::ab_vec(double pm)
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
        k = vec_ind(i, j);
        
        a(k) += pm * 1i*dt* (*it) /2.0;
    }

    return a;
}

// Outputs matrix A or B for CN scheme.
// Input -r, ab_vec(1) for A and r, ab_vec(-1) for B
arma::sp_cx_mat Quantum_box::CN_matAB(arma::cx_double r, arma::cx_vec d)
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

arma::cx_vec Quantum_box::u_init(double xc, double yc, double px, double py, double sigmax, double sigmay)
{
    int N_u = (M-2) * (M-2);
    arma::cx_vec u(N_u);

    double x, y, dx, dy;
    for (int i = 0; i <= M-3; i++)
    {
        y = i*h;
        dy = y - yc;
        for (int j = 0; j <= M-3; j++)
        {
            x = j*h;
            dx = x - xc;
            u(vec_ind(i, j)) = exp( -dx*dx/(2.0*sigmax*sigmax) -dy*dy/(2.0*sigmay*sigmay) + 1i * px * x + 1i * py * y);
        }
    }

    u = arma::normalise(u);
    return u;
}

// Creates potential matrix V containing potential at every point x,y
void Quantum_box::slits_init(int slits)
{
    // Setting correct size of V
    int N = M - 2;
    V.set_size(N,N);

    // Barrier measurements
    double x_pos_barrier = 0.5;
    double barrier_half_width = 0.01;
    double y_middle = 0.5;
    double midbar_half_width = 0.025;
    double aperture = 0.05;
    double v0 = 1e10;

    // First and last x-indice of barrier
    int barrier_xf = floor((x_pos_barrier - barrier_half_width)/h);
    int barrier_xl = ceil((x_pos_barrier + barrier_half_width)/h);
    
    if (slits == 0)
    {
        // This will be the v0 = 0 case, no wall at all
    }
    else if (slits == 1)
    {
        // y-indices of barrier ends
        int y_top_end = floor((y_middle - aperture/2)/h);
        int y_bot_end = ceil((y_middle + aperture/2)/h);

        // Filling V to make potential wall with slits
        for (int i = 0; i <= y_top_end; i++)
        {
            for (int j = barrier_xf; j <= barrier_xl; j++)
            {
                V(i, j) = v0;
            }
        }

        for (int i = y_bot_end; i <= N-1; i++)
        {
            for (int j = barrier_xf; j <= barrier_xl; j++)
            {
                V(i, j) = v0;
            }
        }
    }

    else if(slits == 2)
    {
        // y-indices of barrier ends
        int y_top_end = floor((y_middle-midbar_half_width-aperture)/h);
        int y_midbar_top = floor((y_middle-midbar_half_width)/h);
        int y_midbar_bot = ceil((y_middle+midbar_half_width)/h);
        int y_bot_end = ceil((y_middle+midbar_half_width+aperture)/h);

        // Filling V to make potential wall with slits
        for (int i = 0; i <= y_top_end; i++)
        {
            for (int j = barrier_xf; j <= barrier_xl; j++)
            {
                V(i, j) = v0;
            }
        }

        for (int i = y_midbar_top; i <= y_midbar_bot; i++)
        {
            for (int j = barrier_xf; j <= barrier_xl; j++)
            {
                V(i, j) = v0;
            }
        }

        for (int i = y_bot_end; i <= N-1; i++)
        {
            for (int j = barrier_xf; j <= barrier_xl; j++)
            {
                V(i, j) = v0;
            }
        }
    }

    else if (slits == 3)
    {
        // y-indices of barrier ends
        int y_top_end = floor((y_middle - 2*midbar_half_width - aperture*1.5)/h);
        int y_midbar1_top = floor((y_middle - 2*midbar_half_width - aperture/2)/h);
        int y_midbar1_bot = floor((y_middle - midbar_half_width)/h);
        int y_midbar2_top = ceil((y_middle + midbar_half_width)/h);
        int y_midbar2_bot = ceil((y_middle + 2*midbar_half_width + aperture/2)/h);
        int y_bot_end = ceil((y_middle + 2*midbar_half_width + aperture*1.5)/h);

        // Filling V to make potential wall with slits
        for (int i = 0; i <= y_top_end; i++)
        {
            for (int j = barrier_xf; j <= barrier_xl; j++)
            {
                V(i, j) = v0;
            }
        }

        for (int i = y_midbar1_top; i <= y_midbar1_bot; i++)
        {
            for (int j = barrier_xf; j <= barrier_xl; j++)
            {
                V(i, j) = v0;
            }
        }

            for (int i = y_midbar2_top; i <= y_midbar2_bot; i++)
        {
            for (int j = barrier_xf; j <= barrier_xl; j++)
            {
                V(i, j) = v0;
            }
        }

        for (int i = y_bot_end; i <= N-1; i++)
        {
            for (int j = barrier_xf; j <= barrier_xl; j++)
            {
                V(i, j) = v0;
            }
        }
    }
    else
    {
        cout << "wrong input for slits_init(). Shoud be 0, 1, 2 or 3." << endl;
        exit(1);
    }
}

// Evolves wavefunction one timestep using Crank-Nicolson scheme
void Quantum_box::CN_step()
{
    arma::cx_vec b = B*u;
    u = arma::spsolve(A, b);
}

// Stores wavefuntion at every timestep
void Quantum_box::save_wf(int t_ind)
{
    S.slice(t_ind) = arma::reshape(u, M-2, M-2); // Might need a transpose if results are weird
}

void Quantum_box::deviation(int t_ind)
{
    float d = 1 - arma::norm(u);
    norm_dev(t_ind) = d;
}

void Quantum_box::run_simulation()
{
    int t_steps = ceil(T/dt)+1;
    for (int t_ind = 0; t_ind <= t_steps-1; t_ind++)
    {
        CN_step();
        save_wf(t_ind);
        deviation(t_ind);
    }
}

// A function that prints the structure of a sparse matrix to terminal.
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




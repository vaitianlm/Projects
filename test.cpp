// Execute with:
// g++ test.cpp src/utils.cpp -I include -larmadillo -o test.exe; ./test.exe

#include <cmath>
#include "utils.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    arma::cx_vec salami = arma::cx_vec({1, 2, 3, 4, 5, 6, 7, 8 , 9, 10, 11, 12, 13, 14, 15, 16});
    arma::cx_vec pølse = arma::cx_vec({9, 8, 7, 6, 5, 4, 3, 2, 1});
    Double_slit slit = Double_slit(6, 1, 1, salami, pølse);
    arma::sp_cx_mat A = slit.CN_matAB(slit.r, salami);

    print_sp_matrix_structure(A);
}


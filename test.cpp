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

    arma::cx_vec salami = arma::cx_vec({1, 2, 3, 4});
    arma::cx_vec pølse = arma::cx_vec({4, 3, 2, 1});
    Double_slit slit = Double_slit(1, 1, salami, pølse);

    cout << slit.tridiag(-slit.r, pølse) << endl;
}


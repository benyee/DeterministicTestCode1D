//
//  main.cpp
//  
//
//  Created by Ben Yee on 9/20/13.
//
//
//  Unit tests for various features.

#include <iostream>
#include <string>
#include <sstream>

#include "InputDeck.h"
#include "Utilities.h"
#include "SourceIteration.h"


using namespace std;

int main ()
{
    cout << "Hello world!"<<endl;
    
    
    
    cout << "================ Utilities::solve_ndiag() ===============" << endl;
    unsigned int n = 3;
    unsigned int m = 25;
    vector<vector<double> > A( m, vector<double>(n,0.0) );
    vector<double> b( m , 0.0 );
    
    for(unsigned int i = 0; i < m; i++){
        for(unsigned int j = 0; j < n ; j++)
            A[i][j] = rand() % 100;
        b[i] = rand() % 100;
    }
    
//    cout << "Inputs: " << endl;
//    Utilities::print_dmatrix(A);
//    Utilities::print_dvector(b);
    cout << "Running solve_ndiag..." << endl;
    vector<double> x_n = Utilities::solve_ndiag( A , b , n );
    cout << "Running solve_tridiag..." << endl;
    vector<double> x_3 = Utilities::solve_tridiag( A , b );
//    cout << "Outputs: (should be the same)" << endl;
//    Utilities::print_dvector(x_n);
//    Utilities::print_dvector(x_3);
    
    if( Utilities::p_norm(x_n,x_3,2) > 1e-6 )
        cout << "UNIT TEST FAILED." << endl;
    
    cout << "=========================================================" << endl;



    
    cout << "Goodbye world!"<< endl;
    return 0;
}

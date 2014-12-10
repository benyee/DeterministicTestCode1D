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
    unsigned int n = 5;
    unsigned int m = 25;
    vector<vector<double> > A( m, vector<double>(n,0.0) );
    vector<double> b( m , 0.0 );
    
    for(unsigned int i = 0; i < m; i++){
        for(unsigned int j = 0; j < n ; j++)
            A[i][j] = rand();
        b[i] = rand();
    }
    
    cout << "=========================================================" << endl;



    
    cout << "Goodbye world!"<< endl;
    return 0;
}

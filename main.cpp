//
//  main.cpp
//  
//
//  Created by Ben Yee on 9/20/13.
//
//

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
    
    InputDeck *input = new InputDeck();
    input->setfileName("defaultinput8.txt");
    int debug = input->loadInputDeck();
//    input->readValues();
//    return 0;
    if(debug){
        return 1;
    }
    //input->readValues();
//    input->diffusionSolve();
    SourceIteration *input_run = new SourceIteration(input);
    input_run->iterate(true, false);
    input_run->printOutput(false,20,true);
//    Utilities::print_dmatrix(input_run->get_solution());
    
    /*
     //Testing tridiagonal matrix function:
     vector<vector<double> > A(5, vector<double>(3,0));
     vector<double> b(5,0);
     A[0][1] = -1; A[0][2] = 2.;
     A[1][0] = 3; A[1][1] = 4.; A[1][2] = -5;
     A[2][0] = 6; A[2][1] = 7; A[2][2] = 8;
     A[3][0] = -9; A[3][1] = 10.; A[3][2] = 11.;
     A[4][0] = 12; A[4][1] = -13.;
     b[0] = 5;
     b[1] = -4;
     b[2] = 3;
     b[3] = -2;
     b[4] = 1;
     b = Utilities::solve_tridiag(A,b);
     Utilities::print_dvector(b);
     
     
    //Read in input:
    InputDeck *input = new InputDeck();
    int debug =input->loadInputDeck();
    
    //Test Gaussian quadrature functions:
    vector<double> mu_n = Utilities::calc_mu_n(5);
    Utilities::print_dvector(mu_n);
    Utilities::print_dvector(Utilities::calc_w_n(5));
    Utilities::print_dvector(Utilities::calc_w_n(mu_n));
    */
    
    cout << "Goodbye world!"<<endl;
    return 0;
}

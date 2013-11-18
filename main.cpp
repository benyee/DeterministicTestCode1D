//
//  main.cpp
//  
//
//  Created by Ben Yee on 9/20/13.
//
//

#include <iostream>
#include <string>

#include "InputDeck.h"
#include "Utilities.h"
#include "SourceIteration.h"

using namespace std;

int main ()
{
    cout << "Hello world!"<<endl;
    
    //Generate input files:
    static const double X_arr[] = {5,20,40};
    vector<double> X(X_arr,X_arr+sizeof(X_arr)/sizeof(X_arr[0]));
    static const double SigS_arr[] = {0.5,0.9,0.99,0.999,1};
    vector<double> SigS(SigS_arr,SigS_arr+sizeof(SigS_arr)/sizeof(SigS_arr[0]));
    static const double dx_arr[] = {2,0.5,0.125,0.01};
    vector<double> dx(dx_arr,dx_arr+sizeof(dx_arr)/sizeof(dx_arr[0]));
    static const double N_arr[] = {4};
    vector<double> N(N_arr,N_arr+sizeof(N_arr)/sizeof(N_arr[0]));
    
    for(unsigned int i = 0; i<X.size();i++){
        for(unsigned int j = 0; j<dx.size();j++){
            for(unsigned int k = 0; k<SigS.size();k++){
                
            }
        }
    }
    
    
    
    //Read in input:
    InputDeck *input = new InputDeck();
    int debug =input->loadInputDeck();
    
    //Check input:
    input->readValues();
    if(debug){
        cout << "There is an issue with the sizes of the input vectors!" << endl;
        return 1;
    }
    cout << "Vector sizes look good!"<<endl;
    
    SourceIteration *input_run = new SourceIteration(input);
    input_run->iterate();
    input_run->printOutput();
    /*
    //Test Gaussian quadrature functions:
    vector<double> mu_n = Utilities::calc_mu_n(5);
    Utilities::print_dvector(mu_n);
    Utilities::print_dvector(Utilities::calc_w_n(5));
    Utilities::print_dvector(Utilities::calc_w_n(mu_n));
    */
    
    cout << "Goodbye world!"<<endl;
    return 0;
}

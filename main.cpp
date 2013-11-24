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
    input->setfileName("defaultinput.txt");
    int debug = input->loadInputDeck();
    
    
    //Generate input files:
    static const double X_arr[] = {5,20,40};
    vector<double> X(X_arr,X_arr+sizeof(X_arr)/sizeof(X_arr[0]));
    static const double SigS_arr[] = {0.5,0.9,0.99,0.999,1};
    vector<double> SigS(SigS_arr,SigS_arr+sizeof(SigS_arr)/sizeof(SigS_arr[0]));
    static const double dx_arr[] = {2,0.5,0.125,0.01};
    vector<double> dx(dx_arr,dx_arr+sizeof(dx_arr)/sizeof(dx_arr[0]));
    static const double N_arr[] = {4};
    vector<double> N(N_arr,N_arr+sizeof(N_arr)/sizeof(N_arr[0]));
    /*
    vector<double> temp;
    temp.push_back(X_arr[1]);
    input->setX(temp);
    SourceIteration *input_run = new SourceIteration(input,"testtest.txt");
    input_run->iterate();
    input_run->printOutput();
    */
    
    input->setN(N[0]);
    for(unsigned int i = 0; i<X.size();i++){
        for(unsigned int j = 0; j<dx.size();j++){
            for(unsigned int k = 0; k<SigS.size();k++){
                vector<double> temp;
                temp.push_back(X_arr[i]);
                input->setX(temp);
                
                vector<unsigned int> temp2;
                temp2.push_back(X_arr[i]/dx[j]);
                input->setdiscret(temp2);
                
                vector<double> temp3(temp2[0],0);
                input->setphi_0_0(temp3);
                
                temp[0] = SigS_arr[k];
                input->setsigma_s0(temp);
                temp[0] = 1-SigS_arr[k];
                input->setsigma_a(temp);
                
                //Check input:
                //input->readValues();
                if(debug){
                    cout << "There is an issue with the sizes of the input vectors!" << endl;
                    cout<<"i = "<<i<<", j = "<<j<<", k = "<<k<<endl;
                    return 1;
                }
                //cout << "Vector sizes look good!"<<endl;
                
                ostringstream ss;
                ss<<"OutputFiles/output_I_"<<i<<"_"<<j<<"_"<<k<<"_.txt";
                SourceIteration *input_run = new SourceIteration(input,ss.str());
                input_run->iterate();
                input_run->printOutput();
                cout<<ss.str()<<endl;
                cout<<"-----------------------------------"<<endl;
                
            }
        }
    }
    
    
    /*
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

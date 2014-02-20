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
    input->setfileName("defaultinput1.txt");
    int debug = input->loadInputDeck();
    
    
    /* Problem I*/
    //Generate input files:
    static const double X = 20;
    static const double dx_arr[] = {4,2,1,0.5,0.25,0.125,0.0625, 0.03125, 0.01};
    vector<double> dx(dx_arr,dx_arr+sizeof(dx_arr)/sizeof(dx_arr[0]));
    static const double SigS_arr[] = {1,0.999,0.99,0.95,0.9,0.8};
    vector<double> SigS(SigS_arr,SigS_arr+sizeof(SigS_arr)/sizeof(SigS_arr[0]));
    static const double N_arr[] = {4};
    vector<double> N(N_arr,N_arr+sizeof(N_arr)/sizeof(N_arr[0]));
    static const double alpha_arr[] = {0,2,30};
    vector<double> alpha_mode(alpha_arr,alpha_arr+sizeof(alpha_arr)/sizeof(alpha_arr[0]));
    
    SourceIteration *input_run;
    
    input->setN(N[0]);
    
    vector<double> 
    
    for(unsigned int m = 0; m<alpha_mode.size();m++){
        for(unsigned int j = 0; j<dx.size();j++){
            for(unsigned int k = 0; k<SigS.size();k++){
                cout<<"Running m = "<<m<<", j = "<<j<<", k = "<<k<<endl;
                vector<double> temp;
                temp.push_back(X);
                input->setX(temp);
                
                input->setalpha_mode(alpha_mode[m]);
                
                vector<unsigned int> temp2;
                temp2.push_back(X/dx[j]);
                input->setdiscret(temp2);
                
                vector<double> temp3(temp2[0],0);
                input->setphi_0_0(temp3);
                input->setphi_1_0(temp3);
                input->setphi_0_0_lin(temp3);
                input->setphi_1_0_lin(temp3);
                
                temp[0] = SigS_arr[k];
                input->setsigma_s0(temp);
                temp[0] = 1.0-SigS_arr[k];
                input->setsigma_a(temp);
                
                //Check input:
                //input->readValues();
                if(debug){
                    cout << "There is an issue with the sizes of the input vectors!" << endl;
                    cout<<"m = "<<m<<", j = "<<j<<", k = "<<k<<endl;
                    return 1;
                }
                //cout << "Vector sizes look good!"<<endl;
                
                ostringstream ss;
                if(m == 0){
                    ss<<"OutputFiles/output_9_DD_"<<j<<"_"<<k<<"_.txt";
                }else if (alpha_mode[m] == 30){
                    ss<<"OutputFiles/output_9_new_"<<j<<"_"<<k<<"_.txt";
                }else if (alpha_mode[m] == 33){
                    ss<<"OutputFiles/output_9_new3_"<<j<<"_"<<k<<"_.txt";
                }else if (alpha_mode[m] == 2){
                    ss<<"OutputFiles/output_9_SC_"<<j<<"_"<<k<<"_.txt";
                }else{
                    ss<<"OutputFiles/output_9_"<<alpha_mode[m]<<"_"<<j<<"_"<<k<<"_.txt";
                }
                input_run = new SourceIteration(input,ss.str());
                input_run->iterate();
                input_run->printOutput(false);
                cout<<ss.str()<<endl;
                
                cout<<"Finished running m = "<<m<<", j = "<<j<<", k = "<<k<<endl;
                cout<<"-----------------------------------"<<endl;
                
                delete input_run;
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
    
    std:cout << "Goodbye world!"<<endl;
    return 0;
}

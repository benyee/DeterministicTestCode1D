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
    static const double alpha_arr[] = {0,2,30,33};
    vector<double> alpha_mode(alpha_arr,alpha_arr+sizeof(alpha_arr)/sizeof(alpha_arr[0]));
    
    static const double accel_arr[] = {0,1,2};
    vector<double> accel_mode(accel_arr,accel_arr+sizeof(accel_arr)/sizeof(accel_arr[0]));
    
    SourceIteration *input_run;
    
    input->setN(N[0]);
    double spec_rad[8][10][10];
    double p_esc[8][10][10];
    double it_num[8][10][10][4];
    
    vector<double> mu_n = Utilities::calc_mu_n(N[0]);
    vector<double> w_n = Utilities::calc_w_n(mu_n);
    
    for(unsigned int m = 0; m<alpha_mode.size();m++){
        for(unsigned int j = 0; j<dx.size();j++){
            for(unsigned int k = 0; k<SigS.size();k++){
                for(unsigned int i = 0; i<accel_mode.size();i++){
                    if(i == 0){
                        cout<<"Running i = "<<i<<", m = "<<m<<", j = "<<j<<", k = "<<k<<endl;
                        vector<double> temp;
                        temp.push_back(X);
                        input->setX(temp);
                        
                        input->setalpha_mode(alpha_mode[m]);
                        
                        vector<unsigned int> temp2;
                        temp2.push_back(X/dx[j]);
                        input->setdiscret(temp2);
                        input->setdiscret_CM(temp2);
                        
                        vector<double> temp3(temp2[0],0);
                        input->setphi_0_0(temp3);
                        input->setphi_1_0(temp3);
                        input->setphi_0_0_lin(temp3);
                        input->setphi_1_0_lin(temp3);
                        
                        temp[0] = SigS_arr[k];
                        input->setsigma_s0(temp);
                        temp[0] = 1.0-SigS_arr[k];
                        input->setsigma_a(temp);
                    }
                    input->setaccel_mode(accel_mode[i]);
                
                
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
                    input_run->iterate(0);
    //                input_run->printOutput(false);
                    if(i == 0){
                        spec_rad[m][j][k] = input_run->get_spec_rad();
                        
                        vector<vector<double> > psi_e = input_run->get_psi_e();
                        unsigned int last = psi_e.size()-1;
                        double numerator = (psi_e[last][0]*mu_n[0]*w_n[0] + psi_e[last][1]*mu_n[1]*w_n[1]);
                        double denominator = (psi_e[0][0]*mu_n[0]*w_n[0] + psi_e[0][1]*mu_n[1]*w_n[1]);
                        cout<<psi_e[1][0]<<endl;
                        p_esc[m][j][k] = numerator/denominator;
                    }
                    
                    //input_run->printOutput(false);
                    cout<<ss.str()<<endl;
                    
                    it_num[m][j][k][i] = input_run->get_it_num();
                    
                    cout<<"Finished running m = "<<m<<", j = "<<j<<", k = "<<k<<endl;
                    cout<<"-----------------------------------"<<endl;
                    
                    delete input_run;
                }
            }
        }
    }
    
    for(unsigned int m = 0; m<alpha_mode.size();m++){
        cout<<"------------------------------------------"<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"--SPECTRAL RADIUS FOR alpha_mode ="<<alpha_mode[m]<<" --"<<endl;
        cout<<"Sigma_s = ..."<<'\t';
        for(unsigned int k = 0; k<SigS.size();k++){
            cout<<setprecision(4)<<SigS[k]<<setw(20);
        }
        cout<<" "<<endl;
        for(unsigned int j = 0; j<dx.size();j++){
            cout<<"dx = "<<dx[j]<<": "<<'\t';
            for(unsigned int k = 0; k<SigS.size();k++){
                cout<<spec_rad[m][j][k]<<setw(20);
            }
            cout<<" "<<endl;
        }
        cout<<"------------------------------------------"<<endl;
    }
    for(unsigned int m = 0; m<alpha_mode.size();m++){
        cout<<"------------------------------------------"<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"--p_esc FOR alpha_mode ="<<alpha_mode[m]<<" --"<<endl;
        cout<<"Sigma_s = ..."<<'\t';
        for(unsigned int k = 0; k<SigS.size();k++){
            cout<<setprecision(4)<<SigS[k]<<setw(20);
        }
        cout<<" "<<endl;
        for(unsigned int j = 0; j<dx.size();j++){
            cout<<"dx = "<<dx[j]<<": "<<'\t';
            for(unsigned int k = 0; k<SigS.size();k++){
                cout<<p_esc[m][j][k]<<setw(20);
            }
            cout<<" "<<endl;
        }
        cout<<"------------------------------------------"<<endl;
    }
    for(unsigned int m = 0; m<alpha_mode.size();m++){
        cout<<"------------------------------------------"<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"--p_esc FOR alpha_mode ="<<alpha_mode[m]<<" --"<<endl;
        cout<<"Sigma_s = ..."<<'\t'<<'\t'<<'\t';
        for(unsigned int k = 0; k<SigS.size();k++){
            cout<<setprecision(4)<<SigS[k]<<'\t';
        }
        cout<<" "<<endl;
        for(unsigned int j = 0; j<dx.size();j++){
            for(unsigned int i = 0; i<accel_mode.size();i++){
                cout<<"dx = "<<dx[j]<<", accel_mode = "<<accel_mode[i]<<" : ";
                for(unsigned int k = 0; k<SigS.size();k++){
                    cout<<'\t'<<it_num[m][j][k][i];
                }
                cout<<" "<<endl;
            }
        }
        cout<<"------------------------------------------"<<endl;
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

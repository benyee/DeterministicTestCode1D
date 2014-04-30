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
    int runnum=7;
    cout << "Hello world!"<<endl;
    
    InputDeck *input = new InputDeck();
    input->setfileName("defaultinput7.txt");
    int debug = input->loadInputDeck();
    
    /* Epsilon Problem */
    static const double X = 20;
    static const double dx_arr[] = {4,2,1,0.5,0.25,0.01};
    vector<double> dx(dx_arr,dx_arr+sizeof(dx_arr)/sizeof(dx_arr[0]));
    static const double eps_arr[] = {0.01,0.03,0.1,0.25,0.5,1};
    vector<double> eps(eps_arr,eps_arr+sizeof(eps_arr)/sizeof(eps_arr[0]));
    static const double alpha_arr[] = {2,30};
    vector<double> alpha_mode(alpha_arr,alpha_arr+sizeof(alpha_arr)/sizeof(alpha_arr[0]));
    
    SourceIteration *input_run;
    
    static const int eps_size = sizeof(eps_arr)/sizeof(eps_arr[0]);
    static const int dx_size = sizeof(dx_arr)/sizeof(dx_arr[0]);
    static const int alpha_size =sizeof(alpha_arr)/sizeof(alpha_arr[0]);
    static const int accel_start = 1;
    static const int accel_size = 2;
    int it_num[eps_size][dx_size][alpha_size][accel_size];
    double error[eps_size][dx_size-1][alpha_size];
    
    //Store all the reference solutions here:
    vector<vector<vector<double> > > ref_soln;//ref_soln[soln_#][0 for x, 1 for phi][space #]
    
    bool writeToFile = 1;
    bool printToScreen = 0;
    
    bool solutionFound;
    
    input->setN(4);
    for(unsigned int k = 0; k<eps.size();k++){
        for(int j = dx.size()-1; j>=0;j--){
            for(unsigned int m = 0; m<alpha_mode.size();m++){
                for(int i = accel_size-1; i>=accel_start;i--){
                    if(i == 0 && solutionFound){
                        it_num[k][j][m][i] = 0;
                        continue; //If we already have a solution, don't bother with SI.
                    }
                    if(i==accel_size-1){
                        solutionFound = false;
                    }
                    cout<<"Running k = "<<k<<", j = "<<j<<", m = "<<m<<",i = "<<i<<endl;
                    cout<<"Running eps = "<<eps[k]<<", dx = "<<dx[j]<<", m = "<<m<<",i = "<<i<<endl;
                    input->setaccel_mode(i);
                    
                    if(i==accel_size-1){
                        input->setalpha_mode(alpha_mode[m]);
                        
                        if(m==0){
                            vector<double> temp;
                            temp.push_back(X);
                            input->setX(temp);
                            
                            vector<unsigned int> temp2;
                            temp2.push_back((int)(X/dx[j]));
                            input->setdiscret(temp2);
                            input->setdiscret_CM(temp2);
                            
                            vector<double> temp3(temp2[0],0);
                            input->setphi_0_0(temp3);
                            input->setphi_1_0(temp3);
                            input->setphi_0_0_lin(temp3);
                            input->setphi_1_0_lin(temp3);
                            
                            if(j==dx.size()-1){
                                temp[0] = eps[k];
                                input->setQ(temp);
                                temp[0] = 1.0/eps[k];
                                input->setsigma_s0(temp);
                            }
                        }
                    }
                    //if(i == 0){
                    //    input->settol(1e-10);
                    //}else{
                        //input->settol(1e-6);
                    //}
                    
                    //input->readValues();
                    ostringstream ss;
                    if(alpha_mode[m] == 0){
                        ss<<"OutputFiles/Runset"<<runnum<<"/noscale/output_"<<runnum<<"_DD_"<<k<<"_"<<j<<"_"<<i<<"_.txt";
                    }else if (alpha_mode[m] ==30){
                        ss<<"OutputFiles/Runset"<<runnum<<"/noscale/output_"<<runnum<<"_new_"<<k<<"_"<<j<<"_"<<i<<"_.txt";
                    }else if (alpha_mode[m] ==2){
                        ss<<"OutputFiles/Runset"<<runnum<<"/noscale/output_"<<runnum<<"_SC_"<<k<<"_"<<j<<"_"<<i<<"_.txt";
                    }else{
                        ss<<"OutputFiles/Runset"<<runnum<<"/noscale/output_"<<runnum<<"_"<<alpha_mode[m]<<"_"<<k<<"_"<<j<<"_"<<i<<"_.txt";
                    }
                    
                    //input->diffusionSolve();
//                  cout<<"Finished setting parameters..."<<endl;
                    input_run = new SourceIteration(input,ss.str());
                    input_run->iterate(printToScreen,writeToFile,false);
                    if(writeToFile){
                        input_run->printOutput(false,20);
                        if(input_run->get_it_num()>0 && !solutionFound){
                            cout<<"Found a solution with i = "<<i<<endl;
                            ostringstream new_ss;
                            if(alpha_mode[m] == 0){
                                new_ss<<"OutputFiles/Runset"<<runnum<<"/noscale/output_soln_"<<runnum<<"_DD_"<<k<<"_"<<j<<"_"<<0<<"_.txt";
                            }else if (alpha_mode[m] ==30){
                                new_ss<<"OutputFiles/Runset"<<runnum<<"/noscale/output_soln_"<<runnum<<"_new_"<<k<<"_"<<j<<"_"<<0<<"_.txt";
                            }else if (alpha_mode[m] ==2){
                                new_ss<<"OutputFiles/Runset"<<runnum<<"/noscale/output_soln_"<<runnum<<"_SC_"<<k<<"_"<<j<<"_"<<0<<"_.txt";
                            }else{
                                new_ss<<"OutputFiles/Runset"<<runnum<<"/noscale/output_soln_"<<runnum<<"_"<<alpha_mode[m]<<"_"<<k<<"_"<<j<<"_0_.txt";
                            }
                            input_run->set_outputfilename(new_ss.str());
                            input_run->printOutput(false,20,true);
                            solutionFound = true;
                        }
                    }
                    
                    if(eps_size == 1){
                        Utilities::print_dmatrix(input_run->get_solution());
                    }
                    it_num[k][j][m][i] = input_run->get_it_num();
                    
                    
                    //Check input:
                    if(debug){
                        cout << "There is an issue with the sizes of the input vectors!" << endl;
                        cout<<"m = "<<m<<", j = "<<j<<", k = "<<k<<",i = "<<i<<endl;
                        return 1;
                    }
                    cout<<"-----------------------------------"<<endl;
                    
                    delete input_run;
                }
            }
        }
    }
    
    
    for(unsigned int m = 0; m<alpha_mode.size();m++){
        if(alpha_mode[m] == 0){
            cout<<"Diamond Difference:"<<endl;
        }else if (alpha_mode[m] ==30){
            cout<<"Finite Volume:"<<endl;
        }else if (alpha_mode[m] ==2){
            cout<<"Step Characteristic:"<<endl;
        }else{
            cout<<"m = "<<m<<endl;
        }
        cout<<"dx"<<'\t'<<"Accel. Mode"<<'\t';
        for(unsigned int k  = 0;k<eps.size();k++){
            cout<<"eps="<<eps[k]<<'\t';
        }
        cout<<endl;
        for(unsigned int j = 0; j<dx.size();j++){
            for(unsigned int i = accel_start; i<accel_size;i++){
                cout<<dx[j]<<'\t';
                if(i == 0){
                    cout<<"No Accel."<<'\t';
                }else if(i == 1){
                    cout<<"CMFD"<<'\t';
                }else{
                    cout<<"pCMFD"<<'\t';
                }
                for(unsigned int k = 0; k<eps.size();k++){
                    cout<<it_num[k][j][m][i]<<'\t';
                }
                cout<<endl;
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

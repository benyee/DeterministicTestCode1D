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
    
    static const double c[] = {0.99,0.9,0.8};//{0.99,0.95,0.9,0.8};
//    static const double dx[] = {1.0,3.0,5.0,10.0,15.0};
//    static const double dx[] = {0.01,0.1,1.0};
    static const double dx[] = {0.01,0.05,0.1,0.5,1.0,2.0,3.0,5.0};
    static const unsigned int methods[] = {30,35,40,42};
    
    static const unsigned int c_size = sizeof(c)/sizeof(c[0]);
    static const unsigned int dx_size = sizeof(dx)/sizeof(dx[0]);
    static const unsigned int methods_size = sizeof(methods)/sizeof(methods[0]);
    
    static const double sigma_t = 1.0;
    static const double X = 30;
    
    static const double ref_dx = 1e-3;
    static const unsigned int ref_method = 30;
    
    InputDeck *input = new InputDeck();
    input->setfileName("20141006_input.txt");
    int debug = input->loadInputDeck();
    //    input->readValues();
    //    return 0;
    if(debug){
        return 1;
    }
    
    double last_Jpl_ref[c_size]; //Reference solution for each c value
    double last_Jpl[c_size][dx_size][methods_size];
    double error[c_size][dx_size][methods_size];
    double Qhat_norm[c_size][dx_size][methods_size];
    double Qhat_left[c_size][dx_size][methods_size];
    double Qhat_right[c_size][dx_size][methods_size];
    double rho[c_size][dx_size][methods_size];
    double kappa[c_size][dx_size][methods_size];
    double neg_ang_flux[c_size][dx_size][methods_size];
    double scal_flux[c_size][dx_size][methods_size];
    
    vector<double> w_n;
    
    for(unsigned int i = 0; i < c_size; i++){
        
        vector<double> sigma_s0;
        vector<double> sigma_a;
        sigma_s0.push_back( sigma_t * c[i] );
        sigma_a.push_back( sigma_t * (1 - c[i]) );
        input->setsigma_s0( sigma_s0 );
        input->setsigma_a( sigma_a );
        
        //Obtain the reference solution first:
        
        vector<unsigned int> discret0;
        discret0.push_back(X / ref_dx);
        input->setdiscret(discret0);
        input->setdiscret_CM(discret0);
        
        vector<double> temp0(discret0[0],0);
        input->setphi_0_0(temp0);
        input->setphi_1_0(temp0);
        input->setphi_0_0_lin(temp0);
        input->setphi_1_0_lin(temp0);
        temp0.push_back(0);
        input->setedgePhi0_0(temp0);
        input->setedgePhi1_0(temp0);
        
        input->setalpha_mode(ref_method);
        input->setaccel_mode(1);
        
        cout << "Solving ref soln for c = " << c[i] << endl;
        ostringstream ss_ref;
        ss_ref << "OutputFiles/20141006_Runset/c_" << c[i] << "_.txt";
        SourceIteration *input_run0 = new SourceIteration(input,ss_ref.str());
        input_run0->iterate(false, true, true);
        input_run0->printOutput(false,20);
        
        vector<vector<double> > ref_soln = input_run0 -> get_solution();
        last_Jpl_ref[i] = input_run0 -> get_exitcurr();
        
        input->setaccel_mode(0);
    

        for(unsigned int j = 0; j < dx_size; j++){
            
            vector<unsigned int> discret;
            discret.push_back(X / dx[j]);
            input->setdiscret(discret);
            input->setdiscret_CM(discret);
            
            vector<double> temp(discret[0],0);
            input->setphi_0_0(temp);
            input->setphi_1_0(temp);
            input->setphi_0_0_lin(temp);
            input->setphi_1_0_lin(temp);
            
            temp.push_back(0);
            input->setedgePhi0_0(temp);
            input->setedgePhi1_0(temp);
            
            for(unsigned int k = 0; k < methods_size; k++){
                input->setalpha_mode(methods[k]);
                if(methods[k] == 30){
                    input->setaccel_mode(0);
                }else{
                    input->setaccel_mode(0);
                }
                //input->readValues();
                //    input->diffusionSolve();
                
                ostringstream ss;
                ss << "OutputFiles/20141006_Runset/c_" << c[i] << "_dx_" << dx[j] << "_method_" << methods[k] << "_.txt";
                
                cout << "Running " << ss.str() << endl;
                
                SourceIteration *input_run = new SourceIteration(input,ss.str());
                input_run->iterate(false, true, true);
                input_run->printOutput(false,20);
                
                vector<vector<double> > soln = input_run -> get_solution();
                vector<double> Qhat_edge = input_run -> get_Qhat_edge();
                
                last_Jpl[i][j][k] = input_run -> get_exitcurr();
                
                error[i][j][k] = fabs( (last_Jpl[i][j][k] - last_Jpl_ref[i]) / last_Jpl_ref[i] );
                
                Qhat_norm[i][j][k] = input_run -> get_Qhatnorm();
                Qhat_left[i][j][k] = Qhat_edge[0];
                Qhat_right[i][j][k] = Qhat_edge[Qhat_edge.size()-1];
                rho[i][j][k] = input_run -> get_rho()[0];
                kappa[i][j][k] = input_run -> get_kappa();
                neg_ang_flux[i][j][k] = input_run -> checkNegativeAngularFlux();
                scal_flux[i][j][k] = input_run -> checkNegativeFlux();
            }
            
        }
    }
    
    cout << "======================================="<<endl;
    cout << "error: " << endl;
    for(unsigned int i = 0; i < c_size; i++){
        cout << " c = " << c[i] << endl;
        for(unsigned int j = 0; j < dx_size; j++){
            cout << "[ ";
            for(unsigned int k = 0; k < methods_size; k++){
                cout << " " << error[i][j][k] << " ";
            }
            cout << " ]"<<endl;
        }
    }
    cout << "======================================="<<endl;
    
    cout << "======================================="<<endl;
    cout << "====== last_Jpl ==========" << endl;
    for(unsigned int i = 0; i < c_size; i++){
        cout << " c = " << c[i] << ", last_Jpl_ref[i] = " << last_Jpl_ref[i] << endl;
        for(unsigned int j = 0; j < dx_size; j++){
            cout << "[ ";
            for(unsigned int k = 0; k < methods_size; k++){
                cout << " " << last_Jpl[i][j][k] << " ";
            }
            cout << " ]"<<endl;
        }
    }
    cout << "======================================="<<endl;
    
    cout << "======================================="<<endl;
    cout << "========= Qhat_norm ==========" << endl;
    for(unsigned int i = 0; i < c_size; i++){
        cout << " c = " << c[i] << ", last_Jpl_ref[i] = " << last_Jpl_ref[i] << endl;
        for(unsigned int j = 0; j < dx_size; j++){
            cout << "[ ";
            for(unsigned int k = 0; k < methods_size; k++){
                cout << " " << Qhat_norm[i][j][k] << " ";
            }
            cout << " ]"<<endl;
        }
    }
    
    cout << "======================================="<<endl;
    cout << "========= Qhat_L ==========" << endl;
    for(unsigned int i = 0; i < c_size; i++){
        cout << " c = " << c[i] << ", last_Jpl_ref[i] = " << last_Jpl_ref[i] << endl;
        for(unsigned int j = 0; j < dx_size; j++){
            cout << "[ ";
            for(unsigned int k = 0; k < methods_size; k++){
                cout << " " << Qhat_left[i][j][k] << " ";
            }
            cout << " ]"<<endl;
        }
    }
    
    cout << "======================================="<<endl;
    cout << "========= Qhat_R ==========" << endl;
    for(unsigned int i = 0; i < c_size; i++){
        cout << " c = " << c[i] << ", last_Jpl_ref[i] = " << last_Jpl_ref[i] << endl;
        for(unsigned int j = 0; j < dx_size; j++){
            cout << "[ ";
            for(unsigned int k = 0; k < methods_size; k++){
                cout << " " << Qhat_right[i][j][k] << " ";
            }
            cout << " ]"<<endl;
        }
    }
    
    cout << "======================================="<<endl;
    cout << "========= rho ==========" << endl;
    for(unsigned int i = 0; i < c_size; i++){
        cout << " c = " << c[i] << ", last_Jpl_ref[i] = " << last_Jpl_ref[i] << endl;
        for(unsigned int j = 0; j < dx_size; j++){
            cout << "[ ";
            for(unsigned int k = 0; k < methods_size; k++){
                cout << " " << rho[i][j][k] << " ";
            }
            cout << " ]"<<endl;
        }
    }
    
    cout << "======================================="<<endl;
    cout << "========= kappa ==========" << endl;
    for(unsigned int i = 0; i < c_size; i++){
        cout << " c = " << c[i] << ", last_Jpl_ref[i] = " << last_Jpl_ref[i] << endl;
        for(unsigned int j = 0; j < dx_size; j++){
            cout << "[ ";
            for(unsigned int k = 0; k < methods_size; k++){
                cout << " " << kappa[i][j][k] << " ";
            }
            cout << " ]"<<endl;
        }
    }
    
    cout << "======================================="<<endl;
    cout << "========= # of negative angular fluxes ==========" << endl;
    for(unsigned int i = 0; i < c_size; i++){
        cout << " c = " << c[i] << ", last_Jpl_ref[i] = " << last_Jpl_ref[i] << endl;
        for(unsigned int j = 0; j < dx_size; j++){
            cout << "[ ";
            for(unsigned int k = 0; k < methods_size; k++){
                cout << " " << neg_ang_flux[i][j][k] << " ";
            }
            cout << " ]"<<endl;
        }
    }
    
    cout << "======================================="<<endl;
    cout << "========= # of negative scalar fluxes ==========" << endl;
    for(unsigned int i = 0; i < c_size; i++){
        cout << " c = " << c[i] << ", last_Jpl_ref[i] = " << last_Jpl_ref[i] << endl;
        for(unsigned int j = 0; j < dx_size; j++){
            cout << "[ ";
            for(unsigned int k = 0; k < methods_size; k++){
                cout << " " << scal_flux[i][j][k] << " ";
            }
            cout << " ]"<<endl;
        }
    }
    
    cout << "Goodbye world!"<<endl;
    return 0;
}

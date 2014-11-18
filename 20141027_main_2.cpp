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
    
    static const double eps[] = {1, 0.9, 0.5, 0.4,0.3,0.2,0.1, 0.09,0.08,0.07,0.05,0.04,0.03, 0.01, 0.007, 0.005, 0.003, 0.001, 1e-4, 1e-5, 1e-6};
    static const unsigned int methods[] = {42};
    static const unsigned int accel_modes[] = {1};
    
    static const unsigned int eps_size = sizeof(eps)/sizeof(eps[0]);
    static const unsigned int methods_size = sizeof(methods)/sizeof(methods[0]);
    static const unsigned int accel_modes_size = sizeof(accel_modes)/sizeof(accel_modes);
    
    static const double sigma_t0 = 1.0;
    static const double sigma_a0 = 0.1;
    static const double Q0 = 0.1;
    static const double X = 10;
    static const double dx = 0.1;
    
    static const double root3 = pow( 3.0 , 0.5 );
    static const double exp_root3 = exp(root3);
    
    InputDeck *input = new InputDeck();
    input->setfileName("InputFiles/20141027_input.txt");
    int debug = input->loadInputDeck();
    if(debug){
        return 1;
    }
    
    cout << "lalal"<<endl;
    vector<double> X_vec(1,X);
    input->setX(X_vec);
    vector<unsigned int> discret(1,X/dx);
    input->setdiscret(discret);
    input->setdiscret_CM(discret);
    cout<<"LALALA" << endl;
    vector<double> zeros_vec(discret[0],0.0);
    input->setphi_0_0(zeros_vec);
    input->setphi_1_0(zeros_vec);
    input->setphi_0_0_lin(zeros_vec);
    input->setphi_1_0_lin(zeros_vec);
    
    double it_num[accel_modes_size][methods_size][eps_size];
    double error[accel_modes_size][methods_size][eps_size];
    double error2norm[accel_modes_size][methods_size][eps_size];
    unsigned int num_neg_fluxes[accel_modes_size][methods_size][eps_size];
    unsigned int num_neg_ang_fluxes[accel_modes_size][methods_size][eps_size];
    
    vector<double> w_n;
    
    for(unsigned int i = 0; i < eps_size; i++){
        vector<double> sigma_s0;
        vector<double> sigma_a;
        vector<double> Q;
        sigma_s0.push_back( sigma_t0 / eps[i]- sigma_a0 * eps[i] );
        sigma_a.push_back( sigma_a0 * eps[i] );
        Q.push_back( Q0 * eps[i] );
        input->setsigma_s0( sigma_s0 );
        input->setsigma_a( sigma_a );
        input->setQ( Q );
        
        double A[2][2];
        A[0][0] = 0;
        A[1][0] = 1;
        A[0][1] = 1;
        A[1][1] = eps[i];
        double b[] = {0,0};
        
        vector<double> diff_soln = Utilities::diffusion_solve( sigma_t0 , sigma_a0 , Q0 , X , dx / 2 , A, b);
        
        
        //Diffusion solution via Wolfram alpha:
        // - y'' + 3*y = 3 ; y'(0) = 0; y(10)+y'(10)/(3/eps) = 0
        //Note that the boundary condition depends on eps
//        vector<double> diff_soln((int)(2*X/dx)+1,0.0);
//        double x_ref = 0.0;
//        for(unsigned int n = 0; n < diff_soln.size() ; n++){
//            double expx = exp(root3*x_ref);
//            diff_soln[n] = (3.0/eps[i]-2*root3) * expx - 3.0/eps[i] * exp( 2 * root3 * (x_ref + 5) ) + (3.0/eps[i]+2*root3)*exp(root3 * (x_ref + 20)) - 3.0/eps[i] * exp(10 * root3);
//            diff_soln[n] /= expx;
//            diff_soln[n] /=  3.0/eps[i] - 2*root3 + (3.0/eps[i]+2*root3)*exp(20*root3);
//            x_ref += dx/2;
//        }
//        Utilities::print_dvector(diff_soln);
        

        for(unsigned int j = 0; j < methods_size; j++){
            
//            vector<unsigned int> discret;
//            discret.push_back(X / dx);
//            input->setdiscret(discret);
//            input->setdiscret_CM(discret);
//            
//            vector<double> temp(discret[0],0);
//            input->setphi_0_0(temp);
//            input->setphi_1_0(temp);
//            input->setphi_0_0_lin(temp);
//            input->setphi_1_0_lin(temp);
//            
//            temp.push_back(0);
//            input->setedgePhi0_0(temp);
//            input->setedgePhi1_0(temp);
            
            input->setalpha_mode(methods[j]);
            
            for(unsigned int k = 0; k < accel_modes_size; k++){
                input->setaccel_mode(accel_modes[k]);
                input->diffusionSolve();
                
                ostringstream ss;
                ss << "OutputFiles/20141027_Runset/eps_" << eps[i] << "_method_" << methods[j] << "_accelmodes_" << accel_modes[k]<< "_.txt";
                
                cout << "Running " << ss.str() << endl;
                
                SourceIteration *input_run = new SourceIteration(input,ss.str());
                input_run->iterate(false, true, true);
                input_run->printOutput(false,20);
                
                vector<vector<double> > soln = input_run -> get_solution();
//                Utilities::print_dmatrix(soln);
                
                it_num[k][j][i] = input_run->get_it_num();
                error[k][j][i] = Utilities::p_norm_of_rel_error(soln[1],diff_soln,2);
                error2norm[k][j][i] = Utilities::p_norm(soln[1],diff_soln,2);
                num_neg_fluxes[k][j][i] = input_run->checkNegativeFlux(true);
                num_neg_ang_fluxes[k][j][i] = input_run->checkNegativeAngularFlux();
                
            }
            
        }
    }

    for(unsigned int i = 0; i < eps_size; i++){
        cout << "eps = " << eps[i];
        cout<< "   , it_num =  " << it_num[0][0][i];
        cout<< "  , error = " << error[0][0][i];
        cout<< "  , error2norm = " << error2norm[0][0][i];
        cout<<"  , num_neg_fluxes= " << num_neg_fluxes[0][0][i];
        cout<<"  , num_neg_ang_fluxes= " << num_neg_ang_fluxes[0][0][i];
        cout<< endl;
    }
    
    cout << "Goodbye world!"<<endl;
    return 0;
}

//
//  main.cpp
//  
//
//  Created by Ben Yee on 9/20/13.
//
//
//  MB-3 problem with step function initial guess, Q = 0, solution is flat

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
    
    bool isPrintingToWindow = false;
    
    static const double dx[] = {0.01, 0.1,0.5,1,2,5};
    static const double c[] = {0.2,0.5,0.7,0.9,0.99};
    
    static const unsigned int dx_size = sizeof(dx)/sizeof(dx[0]);
    static const unsigned int c_size = sizeof(c)/sizeof(c[0]);
    
    static const double sigma_t0 = 1.0;
    static const double X = 50.0;
    static const double Q = 0.0;
    
    int bc[2] = {1,0};
    
    InputDeck *input = new InputDeck();
    input->setfileName("InputFiles/20141117_input.txt");
    int debug = input->loadInputDeck();
    if(debug){
        return 1;
    }
    
    double spec_rad[dx_size][c_size];
    double it_num[dx_size][c_size];
    
    vector<double> Q_vec(1,Q);
    input->setQ(Q_vec);
    vector<double> X_vec(1,X);
    input->setX(X_vec);
    
    input->setbc(bc);
    
    for(unsigned int i = 0; i < dx_size; i++){
        vector<unsigned int> discret(1,X/dx[i]);
        input->setdiscret(discret);
        input->setdiscret_CM(discret);
        
        vector<double> zeros_vec(discret[0],0.0);
        input->setphi_1_0(zeros_vec);
        input->setphi_0_0_lin(zeros_vec);
        input->setphi_1_0_lin(zeros_vec);
        zeros_vec.push_back(0.0);
        input->setedgePhi1_0(zeros_vec);
        
        if(Q == 0){
            vector<double> step_vec(discret[0],2.0);
            for(unsigned int n = 0; n <= discret[0]/2; n++)
                step_vec[n] = 1.0;
            input->setphi_0_0(step_vec);
            step_vec.push_back(2.0);
            input->setedgePhi0_0(step_vec);
        }
        else{
            vector<double> step_vec(discret[0],Q/2);
            for(unsigned int n = 0; n <= discret[0]/2; n++)
                step_vec[n] = 2*Q;
            input->setphi_0_0(step_vec);
            step_vec.push_back(Q/2);
            input->setedgePhi0_0(step_vec);
        }
        

        for(unsigned int j = 0; j < c_size; j++){
            vector<double> sigma_s0;
            vector<double> sigma_a;
            sigma_s0.push_back( sigma_t0 * c[j] );
            sigma_a.push_back( sigma_t0 - sigma_s0[0] );
            input->setsigma_s0( sigma_s0 );
            input->setsigma_a( sigma_a );
            
//            input -> diffusionSolve();
            
//            input -> readValues();
            
            ostringstream ss;
            ss << "OutputFiles/20141112_Runset/dx_" << dx[i] << "_c_" << c[j] <<"_.txt";
            
            cout << "Running " << ss.str() << endl;
            
            SourceIteration *input_run = new SourceIteration(input,ss.str());
            input_run->iterate(isPrintingToWindow);
            
            if(Q == -1){//spec_rad calculation for zero solution:
                vector<vector<double> > soln = input_run->get_solution( 0 );
                vector<vector<double> > old_soln = input_run->get_solution( 1 );
                spec_rad[i][j] = Utilities::p_norm(soln[1],2) / Utilities::p_norm(old_soln[1],2);
            }else{
                spec_rad[i][j] = input_run->get_spec_rad();
            }
            it_num[i][j] = input_run->get_it_num();
            
        }
    }

    cout << "================== spectral radius =====================" << endl;
    cout << "dx ";
    for(unsigned int j=0; j < c_size; j++){
        cout << "    " << c[j];
    }
    cout << endl;
    for(unsigned int i = 0; i < dx_size; i++){
        cout << dx[i];
        for(unsigned int j=0; j < c_size; j++){
            cout << "    " << spec_rad[i][j];
        }
        cout << endl;
    }
    cout << "================== spectral radius =====================" << endl;
    cout << "================== iteration number =====================" << endl;
    cout << "dx ";
    for(unsigned int j=0; j < c_size; j++){
        cout << "    " << c[j];
    }
    cout << endl;
    for(unsigned int i = 0; i < dx_size; i++){
        cout << dx[i];
        for(unsigned int j=0; j < c_size; j++){
            cout << "    " << it_num[i][j];
        }
        cout << endl;
    }
    cout << "================== iteration number =====================" << endl;
    
    cout << "Goodbye world!"<<endl;
    return 0;
}

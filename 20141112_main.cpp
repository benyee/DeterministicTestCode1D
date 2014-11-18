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
    
    static const double dx[] = {0.1,0.5,1,2,5,10};
    static const double c[] = {0.2,0.5,0.7,0.9,0.99};
    
    static const unsigned int dx_size = sizeof(dx)/sizeof(dx[0]);
    static const unsigned int c_size = sizeof(c)/sizeof(c[0]);
    
    static const double sigma_t0 = 1.0;
    static const double X = 20;
    
    InputDeck *input = new InputDeck();
    input->setfileName("InputFiles/20141112_input.txt");
    int debug = input->loadInputDeck();
    if(debug){
        return 1;
    }
    
    double spec_rad[dx_size][c_size];
    double it_num[dx_size][c_size];
    
    for(unsigned int i = 0; i < dx_size; i++){
        
        vector<double> X_vec(1,X);
        input->setX(X_vec);
        vector<unsigned int> discret(1,X/dx[i]);
        input->setdiscret(discret);
        input->setdiscret_CM(discret);
        vector<double> zeros_vec(discret[0],0.0);
        input->setphi_0_0(zeros_vec);
        input->setphi_1_0(zeros_vec);
        input->setphi_0_0_lin(zeros_vec);
        input->setphi_1_0_lin(zeros_vec);
        zeros_vec.push_back(0);
        input->setedgePhi0_0(zeros_vec);
        input->setedgePhi1_0(zeros_vec);
        

        for(unsigned int j = 0; j < c_size; j++){
            vector<double> sigma_s0;
            vector<double> sigma_a;
            sigma_s0.push_back( sigma_t0 * c[j] );
            sigma_a.push_back( sigma_t0 - sigma_s0[0] );
            input->setsigma_s0( sigma_s0 );
            input->setsigma_a( sigma_a );
            
            input -> diffusionSolve();
            
//            input -> readValues();
            
            ostringstream ss;
            ss << "OutputFiles/20141112_Runset/dx_" << dx[i] << "_c_" << c[j] <<"_.txt";
            
            cout << "Running " << ss.str() << endl;
            
            SourceIteration *input_run = new SourceIteration(input,ss.str());
            input_run->iterate(false);
            
//            Utilities::print_dmatrix(input_run->get_solution());

            spec_rad[i][j] = input_run->get_spec_rad();
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

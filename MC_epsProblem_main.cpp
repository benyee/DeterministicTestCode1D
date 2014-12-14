//
//  main.cpp
//  
//
//  Created by Ben Yee on 12/13/2014.
//
//
//  MB-3 problem with diffusion initial guess, epsilon --> 0, sigma_t0 = 1, c = 0.9
//		sigma_t = 1 / eps, sigma_s =  1 / eps - 0.1 * eps, sigma_a = 0.1 * eps
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
    
    bool isPrintingToWindow = false;
    
    static const double eps[] = { 1.0 , 0.1 , 0.01 , 1.e-3 , 1.e-4 , 1.e-5 , 1.e-6 };
    //Need to show this works for MB-0?
    static const unsigned int alpha[] = { 0 , 30 , 42};
    
    static const unsigned int N = 4;
    
    static const unsigned int eps_size = sizeof(eps)/sizeof(eps[0]);
    static const unsigned int alpha_size = sizeof(alpha)/sizeof(alpha[0]);
    
    cout << "eps = [";
    for(unsigned int j = 0; j < eps_size; j++)
    	cout << eps[j] << " , ";
    cout << "]"<< endl;
    
    static const double sigma_t0 = 1.0;
    static const double X = 10.0;
    static const double dx = 1.0;
    static const double Q = 0.0;
    static const double c = 0.95;
    double sigma_a0 = (1 - c ) * sigma_t0;
    
    int bc[2] = {2,0};
    
    vector<double> psi_bl;
    for(unsigned int i = 0; i < N/2-1 ; i++)
		psi_bl.push_back(0.0);
	psi_bl.push_back(1.0);
	       
    InputDeck *input = new InputDeck();
    input->setfileName("InputFiles/MC_epsproblem_input.txt");
    int debug = input->loadInputDeck();
    if(debug){
        return 1;
    }
    
    double spec_rad[eps_size][alpha_size];
    double it_num[eps_size][alpha_size];
    
    //----------------------------------------------------
    // Diffusion solver to get a solution to compare to:
    //---------------------------------------------------- 
       
    double A_diff_bc[2][2];
    double b_diff_bc[2];
    A_diff_bc[0][0] = 1.0;
    A_diff_bc[0][1] = -2. / 3 / sigma_t0;
    b_diff_bc[0] = 0;
    vector<double> mu_n = Utilities::calc_mu_n(N);
    vector<double> w_n = Utilities::calc_w_n(mu_n);
    for(unsigned int m = 0; m < N/2 ; m++)
    	b_diff_bc[0] += psi_bl[m] * w_n[m] * mu_n[m];
    b_diff_bc[0] *= 4.0;
    A_diff_bc[1][0] = 1.0;
    b_diff_bc[1] = 0;
    
	//----------------------------------------------------
	//----------------------------------------------------
	
    
    vector<double> Q_vec(1,Q);
    input->setQ(Q_vec);
    vector<double> X_vec(1,X);
    input->setX(X_vec);
    
    input->setbc(bc);
    input->setpsi_bl(psi_bl);
    
	vector<unsigned int> discret(1,X/dx);
	input->setdiscret(discret);
	input->setdiscret_CM(discret);
    
    
    for(unsigned int i = 0; i < eps_size; i++){  
        
		A_diff_bc[0][1] = - 2. / 3 / sigma_t0 * eps[i];
		A_diff_bc[1][1] = 2. / 3 / sigma_t0 * eps[i];
		vector<double> diff_soln = Utilities::diffusion_solve( sigma_t0 , sigma_a0 , Q , X , dx / 2. , A_diff_bc , b_diff_bc );
		cout << "diff_soln["<<i<<"] = " ; Utilities::print_dvector(diff_soln, ',');

        
        vector<double> zeros_vec(discret[0],0.0);
        input->setphi_0_0(zeros_vec);
        input->setphi_1_0(zeros_vec);
        input->setphi_0_0_lin(zeros_vec);
        input->setphi_1_0_lin(zeros_vec);
        zeros_vec.push_back(0.0);
        input->setedgePhi0_0(zeros_vec);
        input->setedgePhi1_0(zeros_vec);
        
		vector<double> sigma_s0;
		vector<double> sigma_a;
		sigma_a.push_back( sigma_a0 * eps[i] );
		sigma_s0.push_back( sigma_t0 / eps[i] - sigma_a[0]  );
		input->setsigma_s0( sigma_s0 );
		input->setsigma_a( sigma_a );
            
        for(unsigned int j = 0; j < alpha_size; j++){
		
			input -> diffusionSolve();
        	input -> setalpha_mode( alpha[j] );
            
//             input -> readValues();
            
            ostringstream ss;
            ss << "OutputFiles/20141213_Runset/eps_" << eps[i] << "_alpha_" << alpha[j] <<"_.txt";
            
            cout << "# Running " << ss.str() << endl;
            
            SourceIteration *input_run = new SourceIteration(input,ss.str());
            input_run->iterate(isPrintingToWindow);
            
            if(Q == -1){//spec_rad calculation for zero solution:
                vector<vector<double> > soln = input_run->get_solution( 0 );
                vector<vector<double> > old_soln = input_run->get_solution( 1 );
                spec_rad[i][j] = Utilities::p_norm(soln[1],2) / Utilities::p_norm(old_soln[1],2);
            }else{
                spec_rad[i][j] = input_run->get_spec_rad();
            }
            vector<vector<double> > soln = input_run->get_solution( 0 );
            if( i == 0 && j == 0 ){
                cout << "x = "; Utilities::print_dvector(soln[0],',');
			}
            cout << "phi["<<i<<"]["<<j<<"] = "; Utilities::print_dvector(soln[1],',');
            input_run->printOutput( isPrintingToWindow );
            
            it_num[i][j] = input_run->get_it_num();
            
        }
    }

    cout << "================== spectral radius =====================" << endl;
    cout << "eps ";
    for(unsigned int j=0; j < alpha_size; j++){
        cout << "    " << alpha[j];
    }
    cout << endl;
    for(unsigned int i = 0; i < eps_size; i++){
        cout << eps[i];
        for(unsigned int j=0; j < alpha_size; j++){
            cout << "    " << spec_rad[i][j];
        }
        cout << endl;
    }
    cout << "================== spectral radius =====================" << endl;
    cout << "================== iteration number =====================" << endl;
    cout << "eps ";
    for(unsigned int j=0; j < alpha_size; j++){
        cout << "    " << alpha[j];
    }
    cout << endl;
    for(unsigned int i = 0; i < eps_size; i++){
        cout << eps[i];
        for(unsigned int j=0; j < alpha_size; j++){
            cout << "    " << it_num[i][j];
        }
        cout << endl;
    }
    cout << "================== iteration number =====================" << endl;
    
    cout << "Goodbye world!"<<endl;
    return 0;
}

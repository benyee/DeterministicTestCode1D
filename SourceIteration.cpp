//
//  SourceIteration.cpp
//  
//
//  Created by Ben Yee on 9/21/13.
//
//

#include "SourceIteration.h"

SourceIteration::SourceIteration(InputDeck input){
    data = input;
    phi_0 = data.getphi_0_0();
    phi_1 = data.getphi_1_0();
    mu_n = Utilities::calc_mu_n(data.getN());
    w_n = Utilities::calc_w_n(mu_n);
}

int SourceIteration::iterate(){
    unsigned int N = data.getN();
    int *bc = data.getbc();
    
    vector< vector<double> > source;
    vector<int> discret = data.getdiscret();
    vector<double> sigma_s0 = data.getsigma_s0();
    vector<double> sigma_s1 = data.getsigma_s1();
    vector<double> sigma_a = data.getsigma_a();
    vector<double> Q = data.getQ();
    
    double tol = abs(data.gettol());
    double error;
    
    //Iterate until tolerance is achieved:
    do{
        //Go either right then left or left then right:
        if(bc[0] == 1 && bc[1] == 0){
            leftIteration();
            rightIteration();
        }else{
            rightIteration();
            leftIteration();
        }
        
        //Update phi, calculate source:
        vector<double> old_phi_0(phi_0);
        unsigned int region = 0;
        unsigned int within_region_counter = 0;
        for(unsigned int j = 0; j<phi_0.size();j++){
            phi_0[j] = 0;
            phi_1[j] = 0;
            //Integrate over angle:
            for(unsigned int m = 0; m<N;m++){
                phi_0[j]+= w_n[m]*psi_c[j][m];
                phi_1[j]+= mu_n[m]*w_n[m]*psi_c[j][m];
            }
            //Calculate and update source term:
            for(unsigned int m=0;m<N;m++){
                source[j][m] = (sigma_s0[region]*phi_0[j]+3*mu_n[m]*sigma_s1[region]*phi_1[j]+Q[region])/2;
            }
            within_region_counter++;
            if(within_region_counter==discret[region]){
                within_region_counter = 0;
                region++;
            }
        }
        error = Utilities::inf_norm(old_phi_0-phi_0);
        cout<<error<<endl;
        error = 0;
    }while(error>tol);
    return 0;
}

void SourceIteration::printOutput(){
    cout<<"Printing to a file named"<<Utilities::remove_ext(data.getfileName())<<endl;
    cout<<"Output goes here"<<endl;
}

void leftIteration(){
    return;
}
void rightIteration(){
    return;
}
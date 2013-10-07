//
//  SourceIteration.cpp
//  
//
//  Created by Ben Yee on 9/21/13.
//
//

#include "SourceIteration.h"

SourceIteration::SourceIteration(InputDeck *input){
    data = input;
    phi_0 = data->getphi_0_0();
    phi_1 = data->getphi_1_0();
    J = phi_0.size();
    bc = data->getbc();
    
    sigma_s0 = data->getsigma_s0();
    sigma_s1 = data->getsigma_s1();
    sigma_a = data->getsigma_a();
    
    N = data->getN();
    mu_n = Utilities::calc_mu_n(N);
    w_n = Utilities::calc_w_n(mu_n);
    
    //Initialize vectors:
    //Note that psi_e has one extra set of elements compared to psi_c and source
    for(unsigned int j = 0; j<=J;j++){
        if(j!=J){
            vector<double> temp;
            vector<double> temp2;
            psi_c.push_back(temp);
            source.push_back(temp2);
        }
        vector<double> temp3;
        psi_e.push_back(temp3);
        for(unsigned int m=0; m<N;m++){
            if(j!=J){
                psi_c[j].push_back(0);
                source[j].push_back(0);
            }
            psi_e[j].push_back(0);
        }
    }
    initializeGrid();
    initializeAlpha();
    
    updatePhi_calcSource();
}

int SourceIteration::iterate(){
    double tol = abs(data->gettol());
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
        
        finiteDifference();
        
        error = updatePhi_calcSource();
        cout<<"error is "<<error<<endl;
    }while(error>tol);
    return 0;
}

void SourceIteration::printOutput(string outfilename){
    ofstream outfile;
    outfile.open (outfilename.c_str());
    cout<<"Writing to file named "+outfilename<<endl;
    outfile<<"<--phi_0-->\n";
    cout<<"<--phi_0-->"<<endl;
    for(unsigned int j = 0; j<J; j++){
        outfile<<phi_0[j]<<'\n';
        cout<<phi_0[j]<<endl;
    }
    outfile<<"<--phi_1-->\n";
    cout<<"<--phi_1-->"<<endl;
    for(unsigned int j = 0; j<J; j++){
        outfile<<phi_1[j]<<'\n';
        cout<<phi_1[j]<<endl;
    }
    outfile<<"<--psi_c-->\n";
    cout<<"<--psi_c-->"<<endl;
    for(unsigned int j = 0; j<J; j++){
        for(unsigned int m = 0; m<N;m++){
            outfile<<psi_c[j][m]<<'\t';
            cout<<psi_c[j][m]<<'\t';
        }
        outfile<<'\n';
        cout<<endl;
    }
    outfile<<"<--psi_e-->\n";
    cout<<"<--psi_e-->"<<endl;
    for(unsigned int j = 0; j<=J; j++){
        for(unsigned int m = 0; m<N;m++){
            outfile<<psi_e[j][m]<<'\t';
            cout<<psi_e[j][m]<<'\t';
        }
        outfile<<'\n';
        cout<<endl;
    }
    outfile.close();
}

void SourceIteration::leftIteration(){
    if(bc[1]==1){
        for(unsigned int m=0; m<N/2;m++){
            psi_e[J][N-m-1] = psi_e[J][m];
        }
    }else if(bc[1]==2){
        vector<double> temp = data->getpsi_br();
        for(unsigned int m=0; m<N/2;m++){
            psi_e[J][m+N/2] = temp[m];
        }
    }
    int region;
    unsigned int within_region_counter;
    for(unsigned int m = N/2; m<N;m++){
        region = discret.size()-1;
        within_region_counter=0;
        for(int j = J-1; j>=0;j--){
            double numerator = (-mu_n[m]-(sigma_s0[region]+sigma_a[region])*h[j]/2.0*(1.0+alpha[j][m]))*psi_e[j+1][m]+source[j][m]*h[j]/2.0;
            double denominator = -mu_n[m]+(sigma_s0[region]+sigma_a[region])*h[j]/2.0*(1.0-alpha[j][m]);
            psi_e[j][m] = numerator/denominator;
            within_region_counter++;
            if(within_region_counter==discret[region]){
                within_region_counter = 0;
                region--;
            }
        }
    }
    return;
}
void SourceIteration::rightIteration(){
    //Reflective boundary conditions:
    if(bc[0]==1){
        for(unsigned int m=0; m<N/2;m++){
            psi_e[0][m] = psi_e[0][N-m-1];
        }
    }else if(bc[0]==2){
        vector<double> temp = data->getpsi_bl();
        for(unsigned int m=0; m<N/2;m++){
            psi_e[0][m] = temp[m];
        }
    }
    int region;
    unsigned int within_region_counter;
    for(unsigned int m = 0; m<N/2;m++){
        region = 0;
        within_region_counter=0;
        for(unsigned int j = 0; j<J;j++){
            double numerator = (mu_n[m]-(sigma_s0[region]+sigma_a[region])*h[j]/2.0*(1.0-alpha[j][m]))*psi_e[j][m]+source[j][m]*h[j]/2.0;
            double denominator = mu_n[m]+(sigma_s0[region]+sigma_a[region])*h[j]/2.0*(1.0+alpha[j][m]);
            psi_e[j+1][m] = numerator/denominator;
            within_region_counter++;
            if(within_region_counter==discret[region]){
                within_region_counter = 0;
                region++;
            }
        }
    }
    return;
}

void SourceIteration::finiteDifference(){
    for(unsigned int j = 0; j<J;j++){
        for(unsigned int m = 0; m<N;m++){
            psi_c[j][m] = ((1.0+alpha[j][m])*psi_e[j+1][m]+(1.0-alpha[j][m])*psi_e[j][m])/2.0;
        }
    }
    return;
}

void SourceIteration::initializeGrid(){
    discret = data->getdiscret();
    if(x.size()){return;}
    double tempL = 0;
    double tempR;
    vector<double> X = data->getX();
    for(unsigned int i = 0; i < X.size();i++){
        tempR = X[i];
        for(unsigned int j = 0; j<discret[i];j++){
            x_e.push_back(j*(tempR-tempL)/discret[i]);
            x.push_back((j+0.5)*(tempR-tempL)/discret[i]);
            h.push_back((tempR-tempL)/discret[i]);
        }
        tempL = tempR;
    }
    x_e.push_back(X[X.size()-1]);
}

void SourceIteration::initializeAlpha(){
    if(alpha.size()){return;}
    
    vector<double> sigma_s0 = data->getsigma_s0();
    vector<double> sigma_a = data->getsigma_a();
    vector<double> sigma_t = Utilities::vector_add(sigma_s0,sigma_a);
    unsigned int region = 0;
    unsigned int within_region_counter = 0;
    for(unsigned int j = 0; j<=J;j++){
        vector<double> temp;
        alpha.push_back(temp);
        if(data->getalpha_mode() == 1){
            for(unsigned int m=0; m<N/2;m++){
                alpha[j].push_back(1);
            }
            for(unsigned int m=N/2;m<N;m++){
                alpha[j].push_back(-1);
            }
        }else if(data->getalpha_mode()==0){
            for(unsigned int m=0; m<N;m++){
                alpha[j].push_back(0);
            }
        }else{
            for(unsigned int m=0; m<N;m++){
                double tau = sigma_t[region]*h[j]/2/mu_n[m];
                alpha[j].push_back(1/tanh(tau)-1/tau);
                within_region_counter++;
            }
        }
        if(within_region_counter == discret[region]){
            within_region_counter=0;
            region++;
        }
    }
}

//Update phi, calculate source:
double SourceIteration::updatePhi_calcSource(){
    vector<double> old_phi_0(phi_0);
    vector<double> Q = data->getQ();
    int region = 0;
    unsigned int within_region_counter = 0;
    for(unsigned int j = 0; j<J;j++){
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
    return (Utilities::p_norm(old_phi_0,phi_0,2));
    
}
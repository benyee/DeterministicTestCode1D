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
    mu_n = Utilities::calc_mu_n(data->getN());
    w_n = Utilities::calc_w_n(mu_n);
    
    //Initialize vectors:
    //Note that psi_e has one extra set of elements compared to psi_c and source
    for(unsigned int i = 0; i<=phi_0.size();i++){
        if(i!=phi_0.size()){
            vector<double> temp;
            vector<double> temp2;
            psi_c.push_back(temp);
            source.push_back(temp2);
        }
        vector<double> temp3;
        psi_e.push_back(temp3);
        for(unsigned int j=0; j<mu_n.size();j++){
            if(i!=phi_0.size()){
                psi_c[i].push_back(0);
                source[i].push_back(0);
            }
            psi_e[i].push_back(0);
        }
    }
    initializeGrid();
    initializeAlpha();
}

int SourceIteration::iterate(){
    unsigned int N = data->getN();
    int *bc = data->getbc();
    
    vector<unsigned int> discret = data->getdiscret();
    vector<double> sigma_s0 = data->getsigma_s0();
    vector<double> sigma_s1 = data->getsigma_s1();
    vector<double> sigma_a = data->getsigma_a();
    vector<double> Q = data->getQ();
    
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
        error = Utilities::inf_norm(old_phi_0,phi_0);
        cout<<error<<endl;
        error = 0;
    }while(error>tol);
    return 0;
}

void SourceIteration::printOutput(string outfilename){
    ofstream outfile;
    outfile.open (outfilename.c_str());
    cout<<"Writing to file named "+outfilename<<endl;
    outfile<<"<--phi_0-->\n";
    cout<<"<--phi_0-->"<<endl;
    for(unsigned int j = 0; j<phi_0.size(); j++){
        outfile<<phi_0[j]<<'\n';
        cout<<phi_0[j]<<endl;
    }
    outfile<<"<--phi_1-->\n";
    cout<<"<--phi_1-->"<<endl;
    for(unsigned int j = 0; j<phi_1.size(); j++){
        outfile<<phi_1[j]<<'\n';
        cout<<phi_1[j]<<endl;
    }
    outfile<<"<--psi_c-->\n";
    cout<<"<--psi_c-->"<<endl;
    for(unsigned int j = 0; j<psi_c.size(); j++){
        for(unsigned int m = 0; m<psi_c[j].size();m++){
            outfile<<psi_c[j][m]<<'\t';
            cout<<psi_c[j][m]<<'\t';
        }
        outfile<<'\n';
        cout<<endl;
    }
    outfile.close();
}

void SourceIteration::leftIteration(){
    return;
}
void SourceIteration::rightIteration(){
    return;
}

void SourceIteration::finiteDifference(){
    return;
}

void SourceIteration::initializeGrid(){
    if(x.size()){return;}
    double tempL = 0;
    double tempR;
    vector<double> X = data->getX();
    vector<unsigned int> discret = data->getdiscret();
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
    for(unsigned int j = 0; j<=phi_0.size();j++){
        vector<double> temp;
        alpha.push_back(temp);
        if(data->getalpha_mode() == 1){
            for(unsigned int m=0; m<(unsigned int)(mu_n.size()/2);m++){
                alpha[j].push_back(1);
            }
            for(unsigned int m=(unsigned int)(mu_n.size()/2);m<mu_n.size();m++){
                alpha[j].push_back(-1);
            }
        }else if(data->getalpha_mode()==0){
            for(unsigned int m=0; m<mu_n.size();m++){
                alpha[j].push_back(0);
            }
        }else{
            for(unsigned int m=0; m<mu_n.size();m++){
                double tau = sigma_t[j]*h[j]/2/mu_n[m];
                alpha[j].push_back(1/tanh(tau)-1/tau);
            }
        }
    }
}
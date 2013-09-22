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
    
    for(unsigned int i = 0; i<=phi_0.size();i++){
        vector<double> temp;
        vector<double> temp2;
        vector<double> temp3;
        psi_c.push_back(temp);
        psi_e.push_back(temp2);
        source.push_back(temp3);
        for(unsigned int j=0; j<mu_n.size();j++){
            psi_c[i].push_back(0);
            psi_e[i].push_back(0);
            source[i].push_back(0);
        }
    }
}

int SourceIteration::iterate(){
    unsigned int N = data->getN();
    int *bc = data->getbc();
    
    vector<int> discret = data->getdiscret();
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
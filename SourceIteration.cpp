//
//  SourceIteration.cpp
//  
//
//  Created by Ben Yee on 9/21/13.
//
//

#include "SourceIteration.h"

SourceIteration::SourceIteration(InputDeck *input,string outputfilename){
    data = input;
    outfilename = outputfilename;
    
    //Store useful information from input deck:
    phi_0 = data->getphi_0_0();
    phi_1 = data->getphi_1_0();
    J = phi_0.size();
    bc = data->getbc();
    alpha_mode = data->getalpha_mode();
    
    sigma_s0 = data->getsigma_s0();
    sigma_s1 = data->getsigma_s1();
    vector<double> sigma_a = data->getsigma_a();
    sigma_t = Utilities::vector_add(sigma_a,sigma_s0);
    
    N = data->getN();
    mu_n = Utilities::calc_mu_n(N);
    w_n = Utilities::calc_w_n(mu_n);
    
    hasLinearTerms = data->gethasLinearTerms();
    if(hasLinearTerms){
        phi_0_lin = data->getphi_0_0_lin();
        phi_1_lin = data->getphi_1_0_lin();
        for(unsigned int j = 0; j<J;j++){
            vector<double> temp;
            vector<double> temp2;
            psi_c_lin.push_back(temp);
            source_lin.push_back(temp2);
            for(unsigned int m=0;m<N;m++){
                psi_c_lin[j].push_back(0);
                source_lin[j].push_back(0);
            }
        }
    }
    
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
    ofstream outfile;
    outfile.open(outfilename.c_str());
    cout<<"Performing source iteration..."<<endl;
    outfile<<"Performing source iteration...\n";
    outfile<<setw(5)<<"it_num"<<setw(20)<<"||Change in Flux||"<<setw(20);
    outfile<<"Conv. Rate Est."<<setw(20)<<"Negative Fluxes"<<endl;;
    double tol = abs(data->gettol());
    double error;
    
    //Iterate until tolerance is achieved:
    it_num = 1;
    do{
        //Go either right then left or left then right:
        if(bc[0] == 1 && bc[1] == 0){
            leftIteration();
            rightIteration();
        }else{
            rightIteration();
            leftIteration();
        }
        
        if(alpha_mode < 4){
            finiteDifference();
        }
        
        error = updatePhi_calcSource();
        outfile<<setw(5)<<it_num<<setw(20)<<error<<setw(20);
        if(it_num ==1){
            outfile<<"---";
        }else{
            outfile<<error/old_error;
        }
        outfile<<setw(20)<<checkNegativeFlux()<<'\n';
        cout<<"For iteration number "<<it_num<<", the error is "<<error<<endl;
        old_error = error;
        it_num++;
    }while(error>tol && it_num < 1000);
    
    outfile.close();
    return 0;
}

void SourceIteration::printOutput(unsigned int tabwidth){
    ofstream outfile;
    outfile.open (outfilename.c_str(), ios::app);
    vector<double> phi_0e = calcEdgePhi(0);
    vector<double> phi_1e = calcEdgePhi(1);
    outfile<<setprecision(4);
    cout<<"Writing to file named "+outfilename<<endl;
    outfile<<'\n';
    outfile<<'\n';
    outfile<<'\n';
    outfile<<'\n';
    outfile<<'\n';
    outfile<<"<<<<<--------------OUTPUT:-------------->>>>>\n";
    outfile<<setw(5)<<"x"<<setw(tabwidth)<<"phi_0,c"<<setw(tabwidth)<<"phi_0,e"<<setw(tabwidth)<<"phi_1,c"<<setw(tabwidth)<<"phi_1,e"<<endl;
    cout<<setw(5)<<"x"<<setw(tabwidth)<<"phi_0,c"<<setw(tabwidth)<<"phi_0,e"<<setw(tabwidth)<<"phi_1,c"<<setw(tabwidth)<<"phi_1,e"<<endl;
    for(unsigned int j = 0; j<J; j++){
        outfile<<setw(5)<<fixed<<x_e[j]<<setw(tabwidth)<<"      "<<setw(tabwidth)<<scientific<<phi_0e[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1e[j]<<endl;
        cout<<setw(5)<<fixed<<x_e[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<scientific<<phi_0e[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1e[j]<<endl;
        outfile<<setw(5)<<fixed<<x[j]<<setw(tabwidth)<<scientific<<phi_0[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1[j]<<setw(tabwidth)<<"     "<<endl;
        cout<<setw(5)<<fixed<<x[j]<<setw(tabwidth)<<scientific<<phi_0[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1[j]<<setw(tabwidth)<<"     "<<endl;
    }
    outfile<<setw(5)<<fixed<<x_e[J]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<scientific<<phi_0e[J]<<setw(tabwidth)<<"    "<<setw(tabwidth)<<phi_1e[J]<<endl;
    cout<<setw(5)<<fixed<<x_e[J]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<scientific<<phi_0e[J]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1e[J]<<endl;
    outfile<<'\n';
    outfile<<'\n';
    outfile<<"<--psi_c-->\n";
    //cout<<"<--psi_c-->"<<endl;
    for(unsigned int j = 0; j<J; j++){
        for(unsigned int m = 0; m<N;m++){
            outfile<<setw(tabwidth)<<psi_c[j][m];
            //cout<<setw(tabwidth)<<psi_c[j][m];
        }
        outfile<<'\n';
        //cout<<endl;
    }
    outfile<<'\n';
    outfile<<'\n';
    outfile<<"<--psi_e-->\n";
    //cout<<"<--psi_e-->"<<endl;
    for(unsigned int j = 0; j<=J; j++){
        for(unsigned int m = 0; m<N;m++){
            outfile<<setw(tabwidth)<<psi_e[j][m];
      //      cout<<setw(tabwidth)<<psi_e[j][m];
        }
        outfile<<'\n';
        //cout<<endl;
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
            double halfhj = h[j]/2.0;
            if(alpha_mode==4){//For linear characteristic stuff
                double twosigt = 2*sigma_t[region];
                double musig = mu_n[m]/sigma_t[region];
                double sighmu = h[j]/musig;
                double srctwosigt = source[j][m]/twosigt;
                double srclintwosigt = source_lin[j][m]/twosigt;
                double C0 = psi_e[j+1][m] - srctwosigt - srclintwosigt*(halfhj-musig);
                double esighmu = exp(sighmu);
                psi_c[j][m] = C0/sighmu*(esighmu-1.0) + srctwosigt-srclintwosigt*musig;
                psi_c_lin[j][m] = srclintwosigt - 6.0*C0/sighmu/h[j]*(1.0 + esighmu-2.0/sighmu*(esighmu-1.0));
                psi_e[j][m] = C0*esighmu- srclintwosigt*(halfhj + musig) + srctwosigt;
            }else if(alpha_mode==5){
                double tau = sigma_t[region]*h[j]/2/mu_n[m];
                double numerator = h[j]*h[j]*tau*source_lin[j][m] + 2*h[j]*(3-tau)*source[j][m]-4*mu_n[m]*(3+2*tau)*psi_e[j+1][m];
                double denominator = 4*(4*mu_n[m]*tau-3*mu_n[m]-sigma_t[region]*h[j]*tau);
                psi_e[j][m] = numerator/denominator;
                psi_c[j][m] = source[j][m]/2/sigma_t[region] - (psi_e[j+1][m] - psi_e[j][m])/2/tau;
                psi_c_lin[j][m] = 2*(psi_c[j][m]-psi_e[j][m])/h[j];
            }else{
                double numerator = (-mu_n[m]-(sigma_t[region])*halfhj*(1.0+alpha[j][m]))*psi_e[j+1][m]+source[j][m]*halfhj;
                double denominator = -mu_n[m]+(sigma_t[region])*halfhj*(1.0-alpha[j][m]);
                psi_e[j][m] = numerator/denominator;
            }
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
            double halfhj = h[j]/2.0;
            if(alpha_mode==4){ //For Linear Characteristic stuff
                double twosigt = 2*sigma_t[region];
                double musig = mu_n[m]/sigma_t[region];
                double sighmu = h[j]/musig;
                double srctwosigt = source[j][m]/twosigt;
                double srclintwosigt = source_lin[j][m]/twosigt;
                double C0 = psi_e[j][m] - srctwosigt + srclintwosigt*(halfhj+musig);
                double esighmu = exp(-sighmu);
                psi_c[j][m] = C0/sighmu*(1.0 - esighmu) + srctwosigt-srclintwosigt*musig;
                psi_c_lin[j][m] = srclintwosigt - 6.0*C0/sighmu/h[j]*(1.0 + esighmu+2.0/sighmu*(esighmu-1));
                psi_e[j+1][m] = C0*esighmu+ srclintwosigt*(halfhj - musig) + srctwosigt;
            }else if(alpha_mode==5){
                double tau = sigma_t[region]*h[j]/2/mu_n[m];
                double numerator = h[j]*h[j]*tau*source_lin[j][m] + 2*h[j]*(3+tau)*source[j][m]+4*mu_n[m]*(3-2*tau)*psi_e[j][m];
                double denominator = 4*(4*mu_n[m]*tau+3*mu_n[m]+sigma_t[region]*h[j]*tau);
                psi_e[j+1][m] = numerator/denominator;
                psi_c[j][m] = source[j][m]/2/sigma_t[region] - (psi_e[j+1][m] - psi_e[j][m])/2/tau;
                psi_c_lin[j][m] = 2*(psi_e[j+1][m] - psi_c[j][m])/h[j];
            }else{
                double numerator = (mu_n[m]-(sigma_t[region])*halfhj*(1.0-alpha[j][m]))*psi_e[j][m]+source[j][m]*halfhj;
                double denominator = mu_n[m]+(sigma_t[region])*halfhj*(1.0+alpha[j][m]);
                psi_e[j+1][m] = numerator/denominator;
            }
            
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
            x_e.push_back(tempL+j*(tempR-tempL)/discret[i]);
            x.push_back(tempL+(j+0.5)*(tempR-tempL)/discret[i]);
            h.push_back((tempR-tempL)/discret[i]);
        }
        tempL = tempR;
    }
    x_e.push_back(X[X.size()-1]);
}

void SourceIteration::initializeAlpha(){
    if(alpha.size()){return;}
    if(alpha_mode > 3){return;}
    unsigned int region = 0;
    unsigned int within_region_counter = 0;
    for(unsigned int j = 0; j<=J;j++){
        vector<double> temp;
        alpha.push_back(temp);
        if(alpha_mode == 1){ //Step method
            for(unsigned int m=0; m<N/2;m++){
                alpha[j].push_back(1);
            }
            for(unsigned int m=N/2;m<N;m++){
                alpha[j].push_back(-1);
            }
        }else if(alpha_mode==0){ //Diamond Difference
            for(unsigned int m=0; m<N;m++){
                alpha[j].push_back(0);
            }
        }else if(alpha_mode==2){ //Step characteristic
            for(unsigned int m=0; m<N;m++){
                double tau = sigma_t[region]*h[j]/2/mu_n[m];
                alpha[j].push_back(1/tanh(tau)-1/tau);
                within_region_counter++;
            }
        }else if(alpha_mode==3){ //Characteristic alternative (3.1 of notes)
            for(unsigned int m=0;m<N;m++){
                double tau = sigma_t[region]*h[j]/2/mu_n[m];
                if(m < N/2){
                    alpha[j].push_back(tau+1);
                }else{
                    alpha[j].push_back(tau-1);
                }
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
    vector<double> X = data->getX();
    X.insert(X.begin(),0);
    vector<double> Q = data->getQ();
    vector<double> Q_lin = data->getQ_lin();
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
            //Note that for a linear source, Q[region]+Q_lin[region]*x[j] gives the average external source in that spatial cell
            source[j][m] = sigma_s0[region]*phi_0[j]+3*mu_n[m]*sigma_s1[region]*phi_1[j]+Q[region]+Q_lin[region]*(x[j]-X[region]);
        }
        within_region_counter++;
        if(within_region_counter==discret[region]){
            within_region_counter = 0;
            region++;
        }
    }
    
    if(hasLinearTerms){
        region = 0;
        within_region_counter = 0;
        for(unsigned int j = 0; j<J;j++){
            phi_0_lin[j] = 0;
            phi_1_lin[j] = 0;
            //Integrate over angle:
            for(unsigned int m = 0; m<N;m++){
                phi_0_lin[j]+= w_n[m]*psi_c_lin[j][m];
                phi_1_lin[j]+= mu_n[m]*w_n[m]*psi_c_lin[j][m];
            }
            //Calculate and update source term:
            for(unsigned int m=0;m<N;m++){
                source_lin[j][m] = sigma_s0[region]*phi_0_lin[j]+3*mu_n[m]*sigma_s1[region]*phi_1_lin[j]+Q_lin[region];
            }
            within_region_counter++;
            if(within_region_counter==discret[region]){
                within_region_counter = 0;
                region++;
            }
        }
    }
    return (Utilities::inf_norm(old_phi_0,phi_0));
    
}

//Update phi, calculate source:
vector<double> SourceIteration::calcEdgePhi(int num){
    vector<double> edgePhi;
    if(!num){
        for(unsigned int j = 0; j<J+1;j++){
            edgePhi.push_back(0);
            for(unsigned int m = 0; m<N;m++){
                edgePhi[j]+= w_n[m]*psi_e[j][m];
            }
        }
    }else{
        for(unsigned int j = 0; j<J+1;j++){
            edgePhi.push_back(0);
            for(unsigned int m = 0; m<N;m++){
                edgePhi[j]+= mu_n[m]*w_n[m]*psi_e[j][m];
            }
        }
    }
    return edgePhi;
    
}

unsigned int SourceIteration::checkNegativeFlux(){
    int temp = 0;
    for(unsigned int j = 0; j<J;j++){
        if(phi_0[j] < 0){
            temp++;
        }
    }
    return temp;
}
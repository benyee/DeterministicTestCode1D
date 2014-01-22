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
    accel_mode = data->getaccel_mode();
    
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
            source_lin.push_back(temp);
            for(unsigned int m=0;m<N;m++){
                source_lin[j].push_back(0);
            }
        }
        
        if(alpha_mode !=11){
            for(unsigned int j = 0; j<J;j++){
                vector<double> temp;
                psi_c_lin.push_back(temp);
                for(unsigned int m=0;m<N;m++){
                    psi_c_lin[j].push_back(0);
                }
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
    
    //Find the maximum value of c:
    c = 0;
    for(unsigned int i = 0; i<sigma_t.size();i++){
        double temp = sigma_s0[i]/sigma_t[i];
        if(temp > c){
            c = temp;
        }
    }
    
    updatePhi_calcSource();
}

SourceIteration::~SourceIteration(){
    data = NULL;
}

int SourceIteration::iterate(){
    ofstream outfile;
    outfile.open(outfilename.c_str());
    cout<<"Performing source iteration..."<<endl;
    outfile<<"Performing source iteration...\n";
    outfile<<setw(5)<<"it_num"<<setw(20)<<"||Change in Flux||"<<setw(20);
    outfile<<"Conv. Rate Est."<<setw(20)<<"Negative Fluxes"<<endl;;
    double tol = (1-c)*abs(data->gettol());
    double error = 0;
    
    //Iterate until tolerance is achieved:
    it_num = 0;
    do{
        it_num++;
        //Go either right then left or left then right:
        if(bc[0] == 1 && bc[1] == 0){
            leftIteration();
            rightIteration();
        }else{
            rightIteration();
            leftIteration();
        }
        
        //Special case:
        if(alpha_mode < 10){
            finiteDifference();
        }
        
        
        error = updatePhi_calcSource();
        
        //CMFD acceleration:
        if(accel_mode == 1){
            cmfd();
            updatePhi_calcSource(false);
        }
        
        outfile<<setw(5)<<it_num<<setw(20)<<error<<setw(20);
        if(it_num ==1){
            outfile<<"---";
        }else{
            outfile<<error/old_error;
        }
        outfile<<setw(20)<<checkNegativeFlux()<<'\n';
        cout<<"For iteration number "<<it_num<<", the error is "<<error<<endl;
        old_error = error;
    }while(error>tol && it_num < MAX_IT);
    if(it_num < MAX_IT){
        cout<<"Source iteration converged in "<<it_num<<" iterations"<<endl;
    }else{
        cout<<"Source iteration did NOT converge in "<<MAX_IT<<" iterations"<<endl;
    }
    
    outfile.close();
    return 0;
}

void SourceIteration::printOutput(bool isPrintingToWindow,unsigned int tabwidth){
    ofstream outfile;
    outfile.open (outfilename.c_str(), ios::app);
    vector<double> phi_0e = calcEdgePhi(0);
    vector<double> phi_1e = calcEdgePhi(1);
    outfile<<setprecision(4);
    if(isPrintingToWindow){
        cout<<"Writing to file named "+outfilename<<endl;
    }
    outfile<<'\n';
    outfile<<'\n';
    outfile<<'\n';
    outfile<<'\n';
    outfile<<'\n';
    outfile<<"<<<<<--------------OUTPUT:-------------->>>>>\n";
    outfile<<setw(5)<<"x"<<setw(tabwidth)<<"phi_0,c"<<setw(tabwidth)<<"phi_0,e"<<setw(tabwidth)<<"phi_1,c"<<setw(tabwidth)<<"phi_1,e"<<endl;
    if(isPrintingToWindow){
        cout<<setw(5)<<"x"<<setw(tabwidth)<<"phi_0,c"<<setw(tabwidth)<<"phi_0,e"<<setw(tabwidth)<<"phi_1,c"<<setw(tabwidth)<<"phi_1,e"<<endl;
    }
    for(unsigned int j = 0; j<J; j++){
        outfile<<setw(5)<<fixed<<x_e[j]<<setw(tabwidth)<<"      "<<setw(tabwidth)<<scientific<<phi_0e[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1e[j]<<endl;
        if(isPrintingToWindow){
            cout<<setw(5)<<fixed<<x_e[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<scientific<<phi_0e[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1e[j]<<endl;
        }
        outfile<<setw(5)<<fixed<<x[j]<<setw(tabwidth)<<scientific<<phi_0[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1[j]<<setw(tabwidth)<<"     "<<endl;
        if(isPrintingToWindow){
            cout<<setw(5)<<fixed<<x[j]<<setw(tabwidth)<<scientific<<phi_0[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1[j]<<setw(tabwidth)<<"     "<<endl;
        }
    }
    outfile<<setw(5)<<fixed<<x_e[J]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<scientific<<phi_0e[J]<<setw(tabwidth)<<"    "<<setw(tabwidth)<<phi_1e[J]<<endl;
    if(isPrintingToWindow){
        cout<<setw(5)<<fixed<<x_e[J]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<scientific<<phi_0e[J]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1e[J]<<endl;
    }
    outfile<<'\n';
    outfile<<'\n';
    outfile<<"<--psi_c-->\n";
    if(isPrintingToWindow){
        cout<<"<--psi_c-->"<<endl;
    }
    for(unsigned int j = 0; j<J; j++){
        for(unsigned int m = 0; m<N;m++){
            outfile<<setw(tabwidth)<<psi_c[j][m];
            if(isPrintingToWindow){
                //cout<<setw(tabwidth)<<psi_c[j][m];
            }
        }
        outfile<<'\n';
        if(isPrintingToWindow){
            //cout<<endl;
        }
    }
    outfile<<'\n';
    outfile<<'\n';
    outfile<<"<--psi_e-->\n";
    if(isPrintingToWindow){
        //cout<<"<--psi_e-->"<<endl;
    }
    for(unsigned int j = 0; j<=J; j++){
        for(unsigned int m = 0; m<N;m++){
            outfile<<setw(tabwidth)<<psi_e[j][m];
            if(isPrintingToWindow){
      //      cout<<setw(tabwidth)<<psi_e[j][m];
            }
        }
        outfile<<'\n';
        if(isPrintingToWindow){
            //cout<<endl;
        }
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
            if(alpha_mode==10){//For linear characteristic stuff
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
            }else if(alpha_mode==11){
                double twosigt = 2*sigma_t[region];
                double musig = mu_n[m]/sigma_t[region];
                double sighmu = h[j]/musig;
                double srctwosigt = source[j][m]/twosigt;
                double srclintwosigt = source_lin[j][m]/twosigt;
                double C0 = psi_e[j+1][m] - srctwosigt - srclintwosigt*(halfhj-musig);
                double esighmu = exp(sighmu);
                psi_c[j][m] = C0/sighmu*(esighmu-1.0) + srctwosigt-srclintwosigt*musig;
                psi_e[j][m] = C0*esighmu- srclintwosigt*(halfhj + musig) + srctwosigt;
            }else if(alpha_mode==20){
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
            if(alpha_mode==10){ //For Linear Characteristic stuff
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
            }else if(alpha_mode==11){ //"tilting method" (linear alternative)
                double twosigt = 2*sigma_t[region];
                double musig = mu_n[m]/sigma_t[region];
                double sighmu = h[j]/musig;
                double srctwosigt = source[j][m]/twosigt;
                double srclintwosigt = source_lin[j][m]/twosigt;
                double C0 = psi_e[j][m] - srctwosigt + srclintwosigt*(halfhj+musig);
                double esighmu = exp(-sighmu);
                psi_c[j][m] = C0/sighmu*(1.0 - esighmu) + srctwosigt-srclintwosigt*musig;
                psi_e[j+1][m] = C0*esighmu+ srclintwosigt*(halfhj - musig) + srctwosigt;
            }else if(alpha_mode==20){ // LD
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

void SourceIteration::cmfd(){
    vector<double> phi_0_CM;
    vector<double> Q_CM;
    vector<double> phi_1_e_CM;

    vector<double> Q = data->getQ();
    
    unsigned int FM_index = 0; //This tells you where you are in terms of fine grid cells.
    unsigned int CM_index = 0; //This tells you where you are in terms of coarse grid cells.
    
    //Extrapolate to calculate corase mesh scalar fluxes, currents, and sources:
    //(Equations 2.11b,2.11c,2.11d of Blake/Larsen paper)
    
    //i tracks which material region you're in
    //j tracks the CM cell that you're in within region i
    //k tracks the FM cell that you're in within the j-th CM cell of region i
    for(unsigned int i = 0; i<discret_CM.size(); i++){
        for(unsigned int j = 0; j<discret_CM[i];j++){
            phi_1_e_CM.push_back(0);
            for(unsigned int m = 0; m<N;m++){
                phi_1_e_CM[CM_index] += w_n[m]*mu_n[m]*psi_e[FM_index][m];
            }
            
            phi_0_CM.push_back(0);
            Q_CM.push_back(0);
            
            for(; x_e[FM_index]<x_CM_e[CM_index+1];FM_index++){
                phi_0_CM[CM_index] += phi_0[FM_index]*h[FM_index];
                Q_CM[CM_index] += Q[i]*h[FM_index];
            }
            phi_0_CM[CM_index] /= h_CM[CM_index];
            Q_CM[CM_index] /= h_CM[CM_index];
            
            //std::cout<<"Q_CM["<<CM_index<<"] = "<<Q_CM[CM_index]<<endl;
            
            CM_index++;
        }
    }
    
    //Get the last current value (there's an extra because there's always one more edge than there are cells)
    phi_1_e_CM.push_back(0);
    for(unsigned int m = 0; m<N;m++){
        phi_1_e_CM[CM_index] += w_n[m]*mu_n[m]*psi_e[FM_index][m];
    }
    
    //
    //Calculate correction factors:
    //
    vector<double> D_c(phi_1_e_CM); //Correction factors
    unsigned int c_size = phi_1_e_CM.size()-1;
    unsigned int f_size = psi_e.size()-1;
    
    //Keep track of where we are:
    unsigned int counter = 0; //Number of coarse cells that we've gone through in this region
    unsigned int region_1 = 0; //region to the left of our current edge
    unsigned int region_2 = 0; //region to the right of our current edge
    
    //Left edge: (Equation 2.15a)
    double in_curr = mu_n[0]*psi_e[0][0]*w_n[0];
    for(unsigned int m = 1; m < N/2; m++){
        in_curr += mu_n[m]*psi_e[0][m]*w_n[m];
    }
    D_c[0] = 2*in_curr;
    D_c[0] = (D_c[0]-phi_1_e_CM[0])/phi_0_CM[0];
    
    
    //Middle of slab: (Equation 2.14)
    for(unsigned int i = 1; i < c_size;i++){
        D_c[i] = (phi_1_e_CM[i] + D_actual_CM[i-1]*(phi_0_CM[i]-phi_0_CM[i-1]))/(phi_0_CM[i]+phi_0_CM[i-1]);
        //std::cout<<"D_c["<<i<<"] = "<<D_c[i]<<endl;
    }
    
    //Right edge: (Equation 2.15b)
    double out_curr = -mu_n[N/2]*psi_e[f_size][N/2]*w_n[N/2];
    for(unsigned int m = N/2; m < N; m++){
        out_curr -= mu_n[m]*psi_e[f_size][m]*w_n[m];
    }
    D_c[c_size] = 2*out_curr;
    D_c[c_size] = (D_c[c_size]+phi_1_e_CM[c_size])/phi_0_CM[c_size-1];
    
//    cout<<"in_curr = "<<in_curr<<" and out_curr = "<<out_curr<<endl;
//    Utilities::print_dvector(D_c);
    
    //
    //Formulate tridiagonal matrix for new coarse mesh fluxes and currents:   (Equations 2.16a-d)
    unsigned int phi_0_CM_size = phi_0_CM.size();
    vector<vector<double> > A(phi_0_CM_size,vector<double>(3,0));
    vector<double> phi_0_CM_new(phi_0_CM_size,0);
    
    //Equation 5.20 of my notes:
    phi_0_CM_new[0] = Q_CM[0]*h_CM[0]+2*in_curr;
    //std::cout<<"b["<<0<<"] = "<<phi_0_CM_new[0]<<endl;
    for(unsigned int i = 1; i<phi_0_CM_size-1;i++){
        phi_0_CM_new[i] = Q_CM[i]*h_CM[i];
        //std::cout<<"b["<<i<<"] = "<<phi_0_CM_new[i]<<endl;
    }
    phi_0_CM_new[phi_0_CM_size-1] = Q_CM[phi_0_CM_size-1]*h_CM[phi_0_CM_size-1] + 2*out_curr;
    //std::cout<<"b["<<phi_0_CM_size-1<<"] = "<<phi_0_CM_new[phi_0_CM_size-1]<<endl;
    
    //Equations 5.21 of my notes:
    A[0][1] = D_actual_CM[0]+opt_CM_a[0]+D_c[1]+D_c[0];
    A[0][2] = D_c[1]-D_actual_CM[0];
    for(unsigned int k = 1; k<phi_0_CM_size-1;k++){
        A[k][0] = -D_actual_CM[k-1]-D_c[k];
        A[k][1] = D_actual_CM[k]+D_actual_CM[k-1]+D_c[k+1]-D_c[k]+opt_CM_a[k];
        A[k][2] = -D_actual_CM[k]+D_c[k+1];
    }
    A[phi_0_CM_size-1][0] = -D_actual_CM[phi_0_CM_size-2]-D_c[phi_0_CM_size-1];
    A[phi_0_CM_size-1][1] = D_actual_CM[phi_0_CM_size-2]+D_c[phi_0_CM_size]-D_c[phi_0_CM_size-1]+opt_CM_a[phi_0_CM_size-1];
    
//    Utilities::print_dmatrix(A);
//    Utilities::print_dvector(D_c);
    //Solve A*phi = b:
    phi_0_CM_new = Utilities::solve_tridiag(A,phi_0_CM_new);
    
    //Need to solve for phi_1:
    
    
    //Now we need to update phi_0:
    FM_index = 0; //This tells you where you are in terms of fine grid cells.
    CM_index = 0; //This tells you where you are in terms of coarse grid cells.
    
    //i tracks which material region you're in
    //j tracks the CM cell that you're in within region i
    //k tracks the FM cell that you're in within the j-th CM cell of region i
    for(unsigned int i = 0; i<discret_CM.size(); i++){
        for(unsigned int j = 0; j<discret_CM[i];j++){
            for(; x_e[FM_index]<x_CM_e[CM_index+1];FM_index++){
                phi_0[FM_index] *= phi_0_CM_new[CM_index]/phi_0_CM[CM_index];
//                std::cout<<"Scaling "<<FM_index<<" by "<<phi_0_CM_new[CM_index]/phi_0_CM[CM_index]<<endl;
            }
            CM_index++;
        }
    }
    
    //Update phi_1:
    
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
    
    if(accel_mode){
        discret_CM = data->getdiscret_CM();
        if(x_CM.size()){return;}
        tempL = 0;
        for(unsigned int i = 0; i < X.size();i++){
            tempR = X[i];
            for(unsigned int j = 0; j<discret_CM[i];j++){
                x_CM_e.push_back(tempL+j*(tempR-tempL)/discret_CM[i]);
                x_CM.push_back(tempL+(j+0.5)*(tempR-tempL)/discret_CM[i]);
                h_CM.push_back((tempR-tempL)/discret_CM[i]);
            }
            tempL = tempR;
         }
         x_CM_e.push_back(X[X.size()-1]);
                           
         //Keep track of where we are:
         unsigned int counter, region_1, region_2;
         region_1 = 0; region_2 = 0;
         counter = 1;
         if(counter == discret_CM[0]){
             region_2++;
             counter = 0;
         }
        
         vector<double> sigma_a = data->getsigma_a();
        
         //Middle of slab: (Equation 2.14)
         for(unsigned int i = 1; i < h_CM.size();i++){
             D_actual_CM.push_back(2.0/3.0/( (sigma_t[region_2]-sigma_s1[region_2])*h_CM[i] + (sigma_t[region_1]-sigma_s1[region_1])*h_CM[i-1] ) );
            
             opt_CM_a.push_back(sigma_a[region_2]*h_CM[i]);
             //Keep track of where we are:
             counter++;
             if(region_2 > region_1){
                 region_1 = region_2;
             }
             if(counter == discret_CM[region_2]){
                region_2++;
                counter = 0;
             }
       }
       opt_CM_a.push_back(sigma_a[sigma_a.size()-1]*h_CM[h_CM.size()-1]);
    }
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
double SourceIteration::updatePhi_calcSource(bool usePsi){
    vector<double> old_phi_0(phi_0);
    vector<double> X = data->getX();
    X.insert(X.begin(),0);
    vector<double> Q = data->getQ();
    vector<double> Q_lin = data->getQ_lin();
    int region = 0;
    unsigned int within_region_counter = 0;
    for(unsigned int j = 0; j<J;j++){
        if(usePsi){
            phi_0[j] = 0;
            phi_1[j] = 0;
            //Integrate over angle:
            for(unsigned int m = 0; m<N;m++){
                phi_0[j]+= w_n[m]*psi_c[j][m];
                phi_1[j]+= mu_n[m]*w_n[m]*psi_c[j][m];
            }
        }
        if(!hasLinearTerms){
            //Calculate and update source term if there are no linear terms.  Otherwise, we need to update phi_0_lin and phi_1_lin first...
            for(unsigned int m=0;m<N;m++){
                //Note that for a linear source, Q[region]+Q_lin[region]*x[j] gives the average external source in that spatial cell
                source[j][m] = sigma_s0[region]*phi_0[j]+3*mu_n[m]*sigma_s1[region]*phi_1[j]+Q[region]+Q_lin[region]*(x[j]-X[region]);
            }
        }
        within_region_counter++;
        if(within_region_counter==discret[region]){
            within_region_counter = 0;
            region++;
        }
    }
    
    if(hasLinearTerms){
        if(alpha_mode !=11){
            region = 0;
            within_region_counter = 0;
            for(unsigned int j = 0; j<J;j++){
                if(usePsi){
                    phi_0_lin[j] = 0;
                    phi_1_lin[j] = 0;
                    //Integrate over angle:
                    for(unsigned int m = 0; m<N;m++){
                        phi_0_lin[j]+= w_n[m]*psi_c_lin[j][m];
                        phi_1_lin[j]+= mu_n[m]*w_n[m]*psi_c_lin[j][m];
                    }
                }
                //Calculate and update source term:
                for(unsigned int m=0;m<N;m++){
                    source[j][m] = sigma_s0[region]*phi_0[j]+3*mu_n[m]*sigma_s1[region]*phi_1[j]+Q[region]+Q_lin[region]*(x[j]-X[region]);
                    source_lin[j][m] = sigma_s0[region]*phi_0_lin[j]+3*mu_n[m]*sigma_s1[region]*phi_1_lin[j]+Q_lin[region];
                }
                within_region_counter++;
                if(within_region_counter==discret[region]){
                    within_region_counter = 0;
                    region++;
                }
            }
        }else{
            region = 0;
            within_region_counter = 0;
            vector<double> edgePhi0 = calcEdgePhi(0);
            vector<double> edgePhi1 = calcEdgePhi(1);
            for(unsigned int j = 0; j<J;j++){
                if(edgePhi0[j+1]+edgePhi0[j] == 0){
                    phi_0_lin[j] = 0;
                }else{
                    phi_0_lin[j] = 2*phi_0[j]/h[j]*(edgePhi0[j+1]-edgePhi0[j])/(edgePhi0[j+1]+edgePhi0[j]);
                }
                
                if(edgePhi1[j+1]+edgePhi1[j] == 0){
                    phi_1_lin[j] = 0;
                }else{
                    phi_1_lin[j] = 2*phi_1[j]/h[j]*(edgePhi1[j+1]-edgePhi1[j])/(edgePhi1[j+1]+edgePhi1[j]);
                }
                
                //Calculate and update source term:
                for(unsigned int m=0;m<N;m++){
                    source[j][m] = sigma_s0[region]*phi_0[j]+3*mu_n[m]*sigma_s1[region]*phi_1[j]+Q[region]+Q_lin[region]*(x[j]-X[region]);
                    source_lin[j][m] = sigma_s0[region]*phi_0_lin[j]+3*mu_n[m]*sigma_s1[region]*phi_1_lin[j]+Q_lin[region];
                }
                within_region_counter++;
                if(within_region_counter==discret[region]){
                    within_region_counter = 0;
                    region++;
                }
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
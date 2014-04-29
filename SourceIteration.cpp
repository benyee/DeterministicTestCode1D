//
//  SourceIteration.cpp
//  
//
//  Created by Ben Yee on 9/21/13.
//
//

#include "SourceIteration.h"

//0000

SourceIteration::SourceIteration(InputDeck *input,string outputfilename){
    diverge = 1e5;
    
    data = input;
    outfilename = outputfilename;
    
    //Store useful information from input deck:
    phi_0 = data->getphi_0_0();
    phi_1 = data->getphi_1_0();
    edgePhi0 = data->getedgePhi0_0();
    edgePhi1 = data->getedgePhi1_0();
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
        source_edge.push_back(temp3);
        for(unsigned int m=0; m<N;m++){
            if(j!=J){
                psi_c[j].push_back(0);
                source[j].push_back(0);
            }
            psi_e[j].push_back(0);
            source_edge[j].push_back(0);
        }
    }
    initializeGrid();
    initializeAlpha();
    
    //Find the maximum value of c:
    c = 0;
    double c_eff = 0;
    vector<double> X = data->getX();
    double L = X[X.size()-1];
    for(unsigned int i = 0; i<sigma_t.size();i++){
        double temp = sigma_s0[i]/sigma_t[i];
        if(temp > c){
            c = temp;
        }
        
        double DB2 = (Utilities::PI/L)*(Utilities::PI/L)/3/(sigma_t[i]-sigma_s1[i]);
        temp = sigma_s0[i]/(sigma_t[i]+DB2);
        if(temp > c_eff){
            c_eff = temp;
        }
    }
    if(c >= 1){
        c = c_eff;
        cout<<"c corrected to "<<c<<endl;
    }
    
    updatePhi_calcSource();
}

SourceIteration::~SourceIteration(){
    data = NULL;
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//------------------------CONSTRUCTOR-----------------------------------------//
//------------------------CONSTRUCTOR-----------------------------------------//

//0IIIIIIII
//----------------------------------------------------------------------------//

//------------------------MAIN LOOP-------------------------------------------//
//------------------------MAIN LOOP-------------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
int SourceIteration::iterate(bool isPrintingToWindow,bool isPrintingToFile, bool falseConvCorrection){
    
    //Remove later
    ofstream interSoln;
    if(intermedSoln){
        string name = "InterSoln.txt";
        interSoln.open(name.c_str());
    }
    //Remove later
    
    isConverged = 0;
    if(isPrintingToWindow){
        cout<<"Performing source iteration..."<<endl;
    }
    ofstream outfile;
    if(isPrintingToFile){
        outfile.open(outfilename.c_str());
        outfile<<"Performing source iteration...\n";
        outfile<<setw(5)<<"it_num"<<setw(20)<<"||Change in Flux||"<<setw(20);
        outfile<<"Conv. Rate Est."<<setw(20)<<"Negative Fluxes"<<endl;
        
    }else if(isPrintingToWindow){
        cout<<"WARNING: printOutput function will append instead of overwrite unless otherwise told"<<endl;
    }
    double tol = (1-c*falseConvCorrection)*abs(data->gettol());
    double error = 0;;
    
    //Iterate until tolerance is achieved:
    it_num = 0.5;
    do{
        
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
        
        
        error=updatePhi_calcSource();
        it_num += 0.5;
        
        //CMFD acceleration:
        if(accel_mode == 1){
            cmfd();
            updatePhi_calcSource(false);
        }else if(accel_mode == 2 || accel_mode == 3){
            pcmfd();
            updatePhi_calcSource(false);
        }
        
        it_num += 0.5;
        
        
        if(it_num ==1){
            spec_rad = 0;
        }else{
            spec_rad = (spec_rad*(it_num-2)+error/old_error)/(it_num-1);
        }
        if(isPrintingToFile){
            outfile<<setw(5)<<it_num<<setw(20)<<error<<setw(20);
            if(it_num ==1){
                outfile<<"---";
            }else{
                outfile<<error/old_error;
            }
            outfile<<setw(20)<<checkNegativeFlux()<<'\n';
        }
        
        if(interSoln){
            for(unsigned int i = 0; i<x.size();i++){
                interSoln<<x_e[i]<<'\t'<<edgePhi0[i]<<'\n';
                interSoln<<x[i]<<'\t'<<phi_0[i]<<'\n';
            }
            interSoln<<x_e[x.size()]<<'\t'<<edgePhi0[x.size()]<<'\n';
        }
        
        if(isPrintingToWindow){
            cout<<"For iteration number "<<it_num<<", the error is "<<error<<endl;
        }
        old_error = error;
        if(it_num == 1){
            init_error = error;
        }
    }while(error>tol && it_num < MAX_IT && (it_num < MAX_IT_accel || accel_mode==0) && (it_num < 5 || (error/init_error)<diverge));
    if(Utilities::nan_checker(phi_0)){
        cout<<"The flux has NaN's in it!"<<endl;
    }else if(error/init_error >= diverge){
        cout<<"Source iterationd diverged in "<<it_num<<" iterations"<<endl;
    }else if(it_num < MAX_IT && (it_num < MAX_IT_accel || accel_mode==0)){
        cout<<"Source iteration converged in "<<it_num<<" iterations"<<endl;
        isConverged = 1;
    }else{
        if(accel_mode==0){
            cout<<"Source iteration did NOT converge in "<<MAX_IT<<" iterations"<<endl;
        }else{
            cout<<"Source iteration did NOT converge in "<<MAX_IT_accel<<" iterations"<<endl;
        }
        cout<<"The final error was "<<error<<endl;
    }
    
    if(isPrintingToWindow){
        outfile.close();
    }
    return 0;
}


//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//--------------------MAIN LOOP-----------------------------------------------//
//--------------------MAIN LOOP-----------------------------------------------//

//0GGGGGGG
//----------------------------------------------------------------------------//

//--------------------GET SOLUTION--------------------------------------------//
//--------------------GET SOLUTION--------------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
vector<vector<double> > SourceIteration::get_solution(){
    vector<vector<double> > out;
    
    //Compute the x coordinates of the phi vector:
    vector<double> x_total;
    for(unsigned int i = 0; i<x.size();i++){
        x_total.push_back(x_e[i]);
        x_total.push_back(x[i]);
    }
    x_total.push_back(x_e[x.size()]);
    
    out.push_back(x_total);
    
    //Compute the phi vector (combining both edge and cell averaged)
    vector<double> phi_total;
    vector<double> phi_edge = calcEdgePhi(0);
    for(unsigned int i = 0; i<phi_0.size();i++){
        phi_total.push_back(phi_edge[i]);
        phi_total.push_back(phi_0[i]);
    }
    phi_total.push_back(phi_edge[phi_0.size()]);
    
    out.push_back(phi_total);
    
    return out;
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//--------------------GET SOLUTION--------------------------------------------//
//--------------------GET SOLUTION--------------------------------------------//

//0PPPPPP
//----------------------------------------------------------------------------//

//--------------------PRINT STUFF--------------------------------------------//
//--------------------PRINT STUFF--------------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

void SourceIteration::printOutput(bool isPrintingToWindow,unsigned int tabwidth, bool newFile){
    ofstream outfile;
    if(newFile){
        outfile.open (outfilename.c_str());
    }else{
        outfile.open (outfilename.c_str(), ios::app);
    }
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

void SourceIteration::print_dictionary(){
    
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//--------------------PRINT STUFF--------------------------------------------//
//--------------------PRINT STUFF--------------------------------------------//

//1111111
//1CCCCCC
//----------------------------------------------------------------------------//

//--------------------CALCULATE EDGE PHI-------------------------------------//
//--------------------CALCULATE EDGE PHI-------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
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

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//--------------------CALCULATE EDGE PHI-------------------------------------//
//--------------------CALCULATE EDGE PHI-------------------------------------//

//----------------------------------------------------------------------------//

//--------------------CHECK FOR NEGATIVE FLUX--------------------------------//
//--------------------CHECK FOR NEGATIVE FLUX--------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
unsigned int SourceIteration::checkNegativeFlux(){
    unsigned int temp = 0;
    for(unsigned int j = 0; j<J;j++){
        if(phi_0[j] < 0){
            temp++;
        }
    }
    return temp;
}


//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//--------------------CHECK FOR NEGATIVE FLUX--------------------------------//
//--------------------CHECK FOR NEGATIVE FLUX--------------------------------//

//1CCCCCC
//----------------------------------------------------------------------------//

//--------------------CMFD/PCMFD----------------------------------------------//
//--------------------CMFD/PCMFD----------------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
void SourceIteration::cmfd(){
    vector<double> phi_0_CM;
    vector<double> Q_CM; //Note that this is Q_CM*h_CM
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
            phi_1_e_CM.push_back(edgePhi1[FM_index]);
            
            phi_0_CM.push_back(0);
            Q_CM.push_back(0);
            
            for(; x_e[FM_index]<x_CM_e[CM_index+1];FM_index++){
                phi_0_CM[CM_index] += phi_0[FM_index]*h[FM_index];
                Q_CM[CM_index] += Q[i]*h[FM_index];
            }
            phi_0_CM[CM_index] /= h_CM[CM_index];
            
            CM_index++;
        }
    }
    
    //Get the last current value (there's an extra because there's always one more edge than there are cells)
    phi_1_e_CM.push_back(edgePhi1[FM_index]);
    
    //
    //Calculate correction factors:
    //
    vector<double> D_c(phi_1_e_CM); //Correction factors
    unsigned int c_size = phi_1_e_CM.size()-1;
    unsigned int f_size = psi_e.size()-1;
    
    
    //Left edge: (Equation 2.15a)
    double in_curr = 0;
    if(bc[0] != 1){
        for(unsigned int m = 0; m < N/2; m++){
            in_curr += mu_n[m]*psi_e[0][m]*w_n[m];
        }
        D_c[0] = 2*in_curr;
        D_c[0] = (D_c[0]-phi_1_e_CM[0])/phi_0_CM[0];
    }else{
        D_c[0] = 0;
    }
    
    //Middle of slab: (Equation 2.14)
    for(unsigned int i = 1; i < c_size;i++){
        D_c[i] = (phi_1_e_CM[i] + D_actual_CM[i-1]*(phi_0_CM[i]-phi_0_CM[i-1]))/(phi_0_CM[i]+phi_0_CM[i-1]);
    }
    
    //Right edge: (Equation 2.15b)
    double out_curr = 0;
    if(bc[1] != 1){
        for(unsigned int m = N/2; m < N; m++){
            out_curr -= mu_n[m]*psi_e[f_size][m]*w_n[m];
        }
        D_c[c_size] = 2*out_curr;
        D_c[c_size] = (D_c[c_size]+phi_1_e_CM[c_size])/phi_0_CM[c_size-1];
    }else{
        D_c[c_size] = 0;
    }
    
    //
    //Formulate tridiagonal matrix for new coarse mesh fluxes and currents:   (Equations 2.16a-d)
    unsigned int phi_0_CM_size = phi_0_CM.size();
    vector<vector<double> > A(phi_0_CM_size,vector<double>(3,0));
    vector<double> phi_0_CM_new(phi_0_CM_size,0);
    
    //Equation 5.20 of my notes:
    phi_0_CM_new[0] = Q_CM[0]+2*in_curr;
//    phi_0_CM_new[0] = Q_CM[0]*h_CM[0]+2*in_curr;
    for(unsigned int i = 1; i<phi_0_CM_size-1;i++){
        phi_0_CM_new[i] = Q_CM[i];
//        phi_0_CM_new[i] = Q_CM[i]*h_CM[i];
    }
    phi_0_CM_new[phi_0_CM_size-1] = Q_CM[phi_0_CM_size-1] + 2*out_curr;
//    phi_0_CM_new[phi_0_CM_size-1] = Q_CM[phi_0_CM_size-1]*h_CM[phi_0_CM_size-1] + 2*out_curr;
    
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
    
    //Calculate new edge fluxes:
    vector<double> phi_1_e_CM_new(phi_1_e_CM);
    //Equations 2.16b-d:
    phi_1_e_CM_new[0] = 2*in_curr - D_c[0]*phi_0_CM_new[0];
    for(unsigned int i = 1; i<phi_0_CM_size;i++){
        phi_1_e_CM_new[i] = -D_actual_CM[i]*(phi_0_CM_new[i]-phi_0_CM_new[i-1])+D_c[i]*(phi_0_CM_new[i]+phi_0_CM_new[i-1]);
    }
    phi_1_e_CM_new[phi_0_CM_size] = D_c[phi_0_CM_size]*phi_0_CM_new[phi_0_CM_size-1] - 2*out_curr;
    
    
    //Now we need to update phi_0:
    FM_index = 0; //This tells you where you are in terms of fine grid cells.
    CM_index = 0; //This tells you where you are in terms of coarse grid cells.
    
    //Store a copy of old phi_0:
    vector<double> old_phi_0(phi_0);
    //Store a copy of the old edgePhi1:
    vector<double> old_edgePhi1(edgePhi1);
    
    //i tracks which material region you're in
    //j tracks the CM cell that you're in within region i
    //k tracks the FM cell that you're in within the j-th CM cell of region i
    for(unsigned int i = 0; i<discret_CM.size(); i++){
        for(unsigned int j = 0; j<discret_CM[i];j++){
            for(; x_e[FM_index]<x_CM_e[CM_index+1];FM_index++){
                phi_0[FM_index] *= phi_0_CM_new[CM_index]/phi_0_CM[CM_index];
                phi_1[FM_index] *= phi_0_CM_new[CM_index]/phi_0_CM[CM_index];
                if(phi_1_e_CM[CM_index] != 0){
                    edgePhi1[FM_index] *= phi_1_e_CM_new[CM_index]/phi_1_e_CM[CM_index];
                }
             }
            CM_index++;
        }
    }
    if(phi_1_e_CM[CM_index] != 0){
        edgePhi1[FM_index] *= phi_1_e_CM_new[CM_index]/phi_1_e_CM[CM_index];
    }
    
    unsigned int psize = phi_0.size();
    
    if(alpha_mode >=30 || alpha_mode==11){
        if(EDGE_ACCEL_MODE==1){
            Utilities::print_dvector(edgePhi0);
            if(old_phi_0[0] != 0){
                edgePhi0[0] *= phi_0[0]/old_phi_0[0];
            }
            //Update phi_0_edge:
            for(unsigned int i = 1; i < psize;i++){
                if(old_phi_0[i]+old_phi_0[i-1]!=0){
                    edgePhi0[i] *= (phi_0[i]+phi_0[i-1])/(old_phi_0[i]+old_phi_0[i-1]);
                }
            }
            if(old_phi_0[psize-1] != 0){
                edgePhi0[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
            }
        }else if(EDGE_ACCEL_MODE==2){
            edgePhi0[0] = edgePhi0[0] + phi_0[0]-old_phi_0[0];
            //Update phi_0_edge:
            for(unsigned int i = 1; i < psize;i++){
                if(old_phi_0[i]+old_phi_0[i-1]!=0){
                    edgePhi0[i] = edgePhi0[i]+ (phi_0[i]+phi_0[i-1])/2 - (old_phi_0[i]+old_phi_0[i-1])/2;
                }
            }
            edgePhi0[psize] = edgePhi0[psize] + phi_0[psize-1]-old_phi_0[psize-1];
        }else if(EDGE_ACCEL_MODE==3){
            if(old_phi_0[0] != 0){
                edgePhi0[0] *= phi_0[0]/old_phi_0[0];
            }
            //Update phi_0_edge:
            for(unsigned int i = 1; i < psize;i++){
                if(old_phi_0[i]+old_phi_0[i-1]!=0){
                    edgePhi0[i] *= (phi_0[i]/old_phi_0[i]+ phi_0[i-1]/old_phi_0[i-1])/2;
                }
            }
            if(old_phi_0[psize-1] != 0){
                edgePhi0[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
            }
        }else if(EDGE_ACCEL_MODE==4){//HERE I AM
            //This only works if the fine mesh and coarse mesh are the same.
            //This only works if the regions are all uniform
            double lnSth = log(1 + 2.0/sigma_t[0]/h[0]);
            double Sth = sigma_t[0]*h[0];
            double Ssh = sigma_s0[0]*h[0];
            double Aed = (1 - Sth/2*lnSth)/(1 - Ssh/2*lnSth);
            double Bed =3*(Sth*(Sth/2*lnSth-1)+1)/8/(1 - Ssh/2*lnSth);
            double Ced =h[0]/2*lnSth/(1 - Ssh/2*lnSth);
            double Ded =(Sth*(Sth*lnSth-2)+2)/4;
            double Eed =(3*Sth*(Sth*Sth*log(Sth/(2+Sth))+ 2*Sth-2)+8 )/8;
            double Fed =(2 - Sth*lnSth)/4;
            
            double const1 =(Aed-Bed*Ded/Eed)/2;
            double const2 = Bed/Eed;
            
            double denom = 1.0/(Bed*Fed*Ssh/Eed+1);
            
            double Q_mod =(Ced-Bed*Fed*h[0]/Eed)*Q[0];
            
            //First calculate correction factors:
            vector<double> D_cm(edgePhi0);
            if(bc[0] == 1){
                D_cm[0] = (edgePhi0[0] - ( 2*const1*old_phi_0[0] + Q_mod)/denom )/2/old_phi_0[0];
            }else{
                D_cm[0] = (4*in_curr-2*phi_1_e_CM[0])/edgePhi0[0];
            }
            for(unsigned int i = 1; i <psize;i++){
                D_cm[i] = (edgePhi0[i] - ( const1*(old_phi_0[i] + old_phi_0[i-1]) + const2*old_edgePhi1[i] + Q_mod)\
                           /denom )/(old_phi_0[i]+old_phi_0[i-1]);
            }
            if(bc[1] == 1){
                D_cm[psize] = (edgePhi0[psize] - ( 2*const1*old_phi_0[psize-1] + Q_mod)/denom )/2/old_phi_0[psize-1];
            }else{
                D_cm[psize] = (4*out_curr+2*phi_1_e_CM[psize])/edgePhi0[psize];
            }

            Utilities::print_dvector(D_cm);
//            cout<<"Params"<<endl; //HERE I AM
//            cout<<edgePhi0[2]<<endl;
//            cout<<const1*(old_phi_0[2] + old_phi_0[1])/denom<<endl;
//            cout<<const2*old_edgePhi1[2]/denom<<endl;
            
            //Correct:
            if(bc[0] ==1){
                edgePhi0[0] = ( 2*const1*phi_0[0] + Q_mod )/denom + D_cm[0]*phi_0[0]*2;
            }else{
                edgePhi0[0] = (4*in_curr-2*phi_1_e_CM[0])/D_cm[0];
            }
            for(unsigned int i = 1; i <psize;i++){
                edgePhi0[i] = (const1*( phi_0[i] + phi_0[i-1] ) + const2*edgePhi1[i] + Q_mod)\
                    /denom + D_cm[i]*( phi_0[i] + phi_0[i-1] );
            }
            if(bc[1] == 1){
                edgePhi0[psize] = ( 2*const1*phi_0[psize-1] + Q_mod )/denom + D_cm[psize]*phi_0[psize-1]*2;
            }else{
                edgePhi0[psize] = ( 4*out_curr+2*phi_1_e_CM[psize] )/D_cm[psize];
            }
        }
    }
    
}

void SourceIteration::pcmfd(){
    vector<double> phi_0_CM;
    vector<double> Q_CM;  //Note that this is Q_CM*h_CM
    vector<double> phi_1_e_CML; //Left edge current on coarse mesh
    vector<double> phi_1_e_CMR; //Right edge current on coarse mesh
    
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
            phi_1_e_CMR.push_back(0);
            for(unsigned int m = 0; m<N/2;m++){
                phi_1_e_CMR[CM_index] += w_n[m]*mu_n[m]*psi_e[FM_index][m];
            }
            
            phi_1_e_CML.push_back(0);
            for(unsigned int m = N/2; m<N;m++){
                phi_1_e_CML[CM_index] -= w_n[m]*mu_n[m]*psi_e[FM_index][m];
            }
            
            phi_0_CM.push_back(0);
            Q_CM.push_back(0);
            
            for(; x_e[FM_index]<x_CM_e[CM_index+1];FM_index++){
                phi_0_CM[CM_index] += phi_0[FM_index]*h[FM_index];
                Q_CM[CM_index] += Q[i]*h[FM_index];
            }
            phi_0_CM[CM_index] /= h_CM[CM_index];
            //Q_CM[CM_index] /= h_CM[CM_index];
            
            CM_index++;
        }
    }
    
    //Get the last current value (there's an extra because there's always one more edge than there are cells)
    phi_1_e_CMR.push_back(0);
    for(unsigned int m = 0; m<N/2;m++){
        phi_1_e_CMR[CM_index] += w_n[m]*mu_n[m]*psi_e[FM_index][m];
    }
    phi_1_e_CML.push_back(0);
    for(unsigned int m = N/2; m<N;m++){
        phi_1_e_CML[CM_index] -= w_n[m]*mu_n[m]*psi_e[FM_index][m];
    }
//    Utilities::print_dvector(phi_1_e_CML);
//    Utilities::print_dvector(phi_1_e_CMR);
    
    //
    //Calculate correction factors:
    //
    unsigned int c_size = phi_1_e_CML.size()-1;
    vector<double> D_cL(c_size+1,0); //Correction factors
    vector<double> D_cR(c_size+1,0); //Correction factors
    
    //Left edge: (Equation 2.15a)
    if(bc[0] == 1){ //reflective bc
        D_cL[0] = 0;
    }else{
        D_cL[0] = (phi_1_e_CMR[0]+phi_1_e_CML[0])/phi_0_CM[0];
    }
    
    //Middle of slab: (Equation 2.14)
    for(unsigned int i = 1; i < D_cL.size()-1;i++){
        if(accel_mode ==2){ //pCMFD
            D_cR[i] = (phi_1_e_CMR[i]+D_actual_CM[i-1]*(phi_0_CM[i]-phi_0_CM[i-1])/2.)/phi_0_CM[i-1]; //
            D_cL[i] = (phi_1_e_CML[i]-(D_actual_CM[i-1]*(phi_0_CM[i]-phi_0_CM[i-1])/2.))/phi_0_CM[i]; //
        } else{ //mpCMFD
            D_cR[i] = (2*phi_1_e_CMR[i]+D_actual_CM[i-1]*phi_0_CM[i])/phi_0_CM[i-1]/(0.5 + D_actual_CM[i-1]); //
            D_cL[i] = (2*phi_1_e_CML[i]+D_actual_CM[i-1]*phi_0_CM[i-1])/phi_0_CM[i]/(0.5 + D_actual_CM[i-1]); //
        }
        
//        Utilities::print_dvector(D_cL);
//        Utilities::print_dvector(D_cR);
    }
    
    //Right edge: (Equation 2.15b)
    if(bc[1] ==1){ //reflective bc
        D_cR[c_size] = 0;
    }else{
        D_cR[c_size] = (phi_1_e_CMR[c_size]+phi_1_e_CML[c_size])/phi_0_CM[c_size-1];
    }
    
//    Utilities::print_dvector(phi_0_CM);
//    Utilities::print_dvector(D_actual_CM);
//    Utilities::print_dvector(D_cL);
//    Utilities::print_dvector(D_cR);
    
    //
    //Formulate tridiagonal matrix for new coarse mesh fluxes and currents:   (Equations 2.16a-d)
    unsigned int phi_0_CM_size = phi_0_CM.size();
    vector<vector<double> > A(phi_0_CM_size,vector<double>(3,0));
    vector<double> phi_0_CM_new(phi_0_CM_size,0);
    
    //To deal with reflective bc's:
    double in_curr;
    if(bc[0] ==1){
        in_curr =0;
    }else{
        in_curr =phi_1_e_CMR[0];
    }
    double out_curr;
    if(bc[1] ==1){
        out_curr =0;
    }else{
        out_curr = phi_1_e_CML[c_size];
    }
    
    //Equation 5.20 of my notes:
    phi_0_CM_new[0] = Q_CM[0]+2*in_curr;
//    phi_0_CM_new[0] = Q_CM[0]*h_CM[0]+2*in_curr;
    for(unsigned int i = 1; i<phi_0_CM_size-1;i++){
        phi_0_CM_new[i] = Q_CM[i];
//        phi_0_CM_new[i] = Q_CM[i]*h_CM[i];
    }
    phi_0_CM_new[phi_0_CM_size-1] = Q_CM[phi_0_CM_size-1] + 2*out_curr;
//    phi_0_CM_new[phi_0_CM_size-1] = Q_CM[phi_0_CM_size-1]*h_CM[phi_0_CM_size-1] + 2*out_curr;
    
    //Equations 5.21 of my notes:
    if(accel_mode ==2){ //CMFD
        A[0][1] = D_actual_CM[0]+opt_CM_a[0]+D_cR[1]+D_cL[0];
        A[0][2] = -D_cL[1]-D_actual_CM[0];
        for(unsigned int k = 1; k<phi_0_CM_size-1;k++){
            A[k][0] = -D_actual_CM[k-1]-D_cR[k];
            A[k][1] = D_actual_CM[k]+D_actual_CM[k-1]+D_cR[k+1]+D_cL[k]+opt_CM_a[k];
            A[k][2] = -D_actual_CM[k]-D_cL[k+1];
        }
        A[phi_0_CM_size-1][0] = -D_actual_CM[phi_0_CM_size-2]-D_cR[phi_0_CM_size-1];
        A[phi_0_CM_size-1][1] = D_actual_CM[phi_0_CM_size-2]+D_cR[phi_0_CM_size]+D_cL[phi_0_CM_size-1]+opt_CM_a[phi_0_CM_size-1];
    }else{ //mpCMFD
        //Equations 5.21 of my notes:
        A[0][1] = D_cR[1]/4 + (1 + D_cR[1])*(D_actual_CM[0]/2) + D_cL[0] + opt_CM_a[0];
        A[0][2] = -(D_cL[1]/4+(1+D_cL[1])*(D_actual_CM[0]/2));
        for(unsigned int k = 1; k<phi_0_CM_size-1;k++){
            A[k][0] = -(D_cR[k]/4+(1+D_cR[k])*(D_actual_CM[k-1]/2));
            A[k][1] = D_cL[k]/4+(1+D_cL[k])*(D_actual_CM[k-1]/2)+D_cR[k+1]/4+(1+D_cR[k+1])*(D_actual_CM[k]/2)+opt_CM_a[k];
            A[k][2] = -(D_cL[k+1]/4+(1+D_cL[k+1])*(D_actual_CM[k]/2));
        }
        A[phi_0_CM_size-1][0] = -(D_cR[phi_0_CM_size-1]/4+(1+D_cR[phi_0_CM_size-1])*(D_actual_CM[phi_0_CM_size-2]/2));
        A[phi_0_CM_size-1][1] = D_cL[phi_0_CM_size-1]/4+(1+D_cL[phi_0_CM_size-1])*(D_actual_CM[phi_0_CM_size-2]/2)+D_cR[phi_0_CM_size]+opt_CM_a[phi_0_CM_size-1];
    }
    
//    Utilities::print_dmatrix(A);
//    Utilities::print_dvector(D_c);
    //Solve A*phi = b:
    phi_0_CM_new = Utilities::solve_tridiag(A,phi_0_CM_new);
    
    //Now we need to update phi_0:
    FM_index = 0; //This tells you where you are in terms of fine grid cells.
    CM_index = 0; //This tells you where you are in terms of coarse grid cells.
    
    //Store a copy of old phi_0:
    vector<double> old_phi_0(phi_0);
    
    //i tracks which material region you're in
    //j tracks the CM cell that you're in within region i
    //k tracks the FM cell that you're in within the j-th CM cell of region i
    for(unsigned int i = 0; i<discret_CM.size(); i++){
        for(unsigned int j = 0; j<discret_CM[i];j++){
            for(; x_e[FM_index]<x_CM_e[CM_index+1];FM_index++){
                phi_0[FM_index] *= phi_0_CM_new[CM_index]/phi_0_CM[CM_index];
                phi_1[FM_index] *= phi_0_CM_new[CM_index]/phi_0_CM[CM_index];
            }
            CM_index++;
        }
    }
    
    if(alpha_mode >=30 || alpha_mode==11){
        if(EDGE_ACCEL_MODE==1){
            unsigned int psize = phi_0.size();
            if(old_phi_0[0] != 0){
                edgePhi0[0] *= phi_0[0]/old_phi_0[0];
                edgePhi1[0] *= phi_0[0]/old_phi_0[0];
            }
            //Update phi_0_edge:
            for(unsigned int i = 1; i < psize;i++){
                if(old_phi_0[i]+old_phi_0[i-1]!=0){
                    edgePhi0[i] *= (phi_0[i]+phi_0[i-1])/(old_phi_0[i]+old_phi_0[i-1]);
                    //This is a cheap way to get edgePhi1.  In reality, you can calculate it from the CMFD equations.
                    edgePhi1[i] *= (phi_0[i]+phi_0[i-1])/(old_phi_0[i]+old_phi_0[i-1]);
                }
            }
            if(old_phi_0[psize-1] != 0){
                edgePhi0[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
                edgePhi1[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
            }
        }else if(EDGE_ACCEL_MODE==2){
            unsigned int psize = phi_0.size();
            edgePhi0[0] = edgePhi0[0] + phi_0[0]-old_phi_0[0];
            if(old_phi_0[0] != 0){
                edgePhi1[0] *= phi_0[0]/old_phi_0[0];
            }
            //Update phi_0_edge:
            for(unsigned int i = 1; i < psize;i++){
                if(old_phi_0[i]+old_phi_0[i-1]!=0){
                    edgePhi0[i] = edgePhi0[i]+ (phi_0[i]+phi_0[i-1])/2 - (old_phi_0[i]+old_phi_0[i-1])/2;
                    //This is a cheap way to get edgePhi1.  In reality, you can calculate it from the CMFD equations.
                    edgePhi1[i] *= (phi_0[i]+phi_0[i-1])/(old_phi_0[i]+old_phi_0[i-1]);
                }
            }
            edgePhi0[psize] = edgePhi0[psize] + phi_0[psize-1]-old_phi_0[psize-1];
            if(old_phi_0[psize-1] != 0){
                edgePhi1[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
            }
        }else if(EDGE_ACCEL_MODE==3){
            unsigned int psize = phi_0.size();
            if(old_phi_0[0] != 0){
                edgePhi0[0] *= phi_0[0]/old_phi_0[0];
                edgePhi1[0] *= phi_0[0]/old_phi_0[0];
            }
            //Update phi_0_edge:
            for(unsigned int i = 1; i < psize;i++){
                if(old_phi_0[i]+old_phi_0[i-1]!=0){
                    edgePhi0[i] *= (phi_0[i]/old_phi_0[i]+ phi_0[i-1]/old_phi_0[i-1])/2;
                    //This is a cheap way to get edgePhi1.  In reality, you can calculate it from the CMFD equations.
                    edgePhi1[i] *= (phi_0[i]+phi_0[i-1])/(old_phi_0[i]+old_phi_0[i-1]);
                }
            }
            if(old_phi_0[psize-1] != 0){
                edgePhi0[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
                edgePhi1[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
            }
        }
    }
    
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//--------------------CMFD/PCMFD----------------------------------------------//
//--------------------CMFD/PCMFD----------------------------------------------//


//1FFFFFFF
//----------------------------------------------------------------------------//

//----------------FINITE DIFFERENCE--------------------------------------------//
//----------------FINITE DIFFERENCE--------------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//


void SourceIteration::finiteDifference(){
    for(unsigned int j = 0; j<J;j++){
        for(unsigned int m = 0; m<N;m++){
            psi_c[j][m] = ((1.0+alpha[j][m])*psi_e[j+1][m]+(1.0-alpha[j][m])*psi_e[j][m])/2.0;
        }
    }
    return;
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//--------------------FINITE DIFFERENCE---------------------------------------//
//--------------------FINITE DIFFERENCE---------------------------------------//


//1IIIIIIIIII
//----------------------------------------------------------------------------//

//----------------INITIALIZE STUFF--------------------------------------------//
//----------------INITIALIZE STUFF--------------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

void SourceIteration::initializeAlpha(){
    if(alpha.size()){return;}
    if(alpha_mode > 3){return;}
    unsigned int region = 0;
    unsigned int within_region_counter = 0;
    for(unsigned int j = 0; j<J;j++){
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
                double tau = sigma_t[region]*h[j]/2.0/mu_n[m];
                alpha[j].push_back(1.0/tanh(tau)-1.0/tau);
            }
            within_region_counter++;
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

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//----------------INITIALIZE STUFF--------------------------------------------//
//----------------INITIALIZE STUFF--------------------------------------------//


//1IIIIIIII
//----------------------------------------------------------------------------//

//-----------LEFT AND RIGHT ITERATE------------------------------------------//
//-----------LEFT AND RIGHT ITERATE------------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

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
            }else if(alpha_mode ==30){
                double tau = -sigma_t[region]*h[j]/mu_n[m];
                double numerator = psi_e[j+1][m] - source[j][m]*h[j]/mu_n[m]/2 - source_edge[j][m]*h[j]*tau/4/mu_n[m];
                double denominator = 1 + tau + tau*tau/2;
                psi_e[j][m] = numerator/denominator;
                psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];
                /*
                 double twomu_h = mu_n[m]/halfhj;
                 double mu_hsigt = mu_n[m]/h[j]*sigma_t[region];
                 psi_e[j][m] = (twomu_h*mu_hsigt*psi_e[j+1][m]-mu_hsigt*source[j][m]+source_edge[j][m]/2)/(twomu_h*(mu_hsigt-1)+sigma_t[region]);
                 psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];
                 */
            }else if(alpha_mode==31){
                double tau = -sigma_t[region]*h[j]/mu_n[m];
                double numerator = psi_e[j+1][m] - source[j][m]*h[j]/mu_n[m]/2 - source_edge[j][m]*h[j]/mu_n[m]*tau*(0.25+tau/12);
                //double numerator = psi_e[j+1][m] - source[j][m]*h[j]/mu_n[m]/2 + source_edge[j][m]/sigma_t[region]*(tau*tau/4+tau*tau*tau/12);
                double denominator = 1 + tau + tau*tau/2 + tau*tau*tau/6;
                psi_e[j][m] = numerator/denominator;
                psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];
            }else if(alpha_mode==32){
                double tau = -sigma_t[region]*h[j]/mu_n[m];
                double etau = exp(tau);
                double numerator = psi_e[j+1][m] - source[j][m]*h[j]/mu_n[m]/2 + source_edge[j][m]/sigma_t[region]/2*(etau-1-tau);
                psi_e[j][m] = numerator/etau;
                psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];
            }else if(alpha_mode==33){
                double tau = -sigma_t[region]*h[j]/mu_n[m];
                double tau_a = (sigma_s0[region]-sigma_t[region])*h[j]/mu_n[m];
                double etau = exp(tau_a) - 1 - tau_a - tau_a*tau_a/2;
                double numerator = psi_e[j+1][m] - source[j][m]*h[j]/mu_n[m]/2 + source_edge[j][m]/sigma_t[region]*(tau*tau/4+etau/2);
                double denominator = 1 + tau + tau*tau/2 + etau;
                psi_e[j][m] = numerator/denominator;
                psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];
            }else{ //Regular finite difference
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
                double tau = sigma_t[region]*halfhj/mu_n[m];
                double numerator = h[j]*h[j]*tau*source_lin[j][m] + 2*h[j]*(3+tau)*source[j][m]+4*mu_n[m]*(3-2*tau)*psi_e[j][m];
                double denominator = 4*(4*mu_n[m]*tau+3*mu_n[m]+sigma_t[region]*h[j]*tau);
                psi_e[j+1][m] = numerator/denominator;
                psi_c[j][m] = source[j][m]/2/sigma_t[region] - (psi_e[j+1][m] - psi_e[j][m])/2/tau;
                psi_c_lin[j][m] = 2*(psi_e[j+1][m] - psi_c[j][m])/h[j];
            }else if(alpha_mode ==30){
                double tau = sigma_t[region]*h[j]/mu_n[m];
                double numerator = psi_e[j][m] + source[j][m]*h[j]/mu_n[m]/2 + source_edge[j+1][m]*h[j]/mu_n[m]*tau/4;
                double denominator = 1 + tau + tau*tau/2;
                psi_e[j+1][m] = numerator/denominator;
                psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];
                /*
                 double twomu_h = mu_n[m]/halfhj;
                 double mu_hsigt = mu_n[m]/h[j]*sigma_t[region];
                 psi_e[j+1][m] = (twomu_h*mu_hsigt*psi_e[j][m]+mu_hsigt*source[j][m]+source_edge[j+1][m]/2)/(twomu_h*(mu_hsigt+1)+sigma_t[region]);
                 psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];*/
            }else if(alpha_mode==31){
                double tau = sigma_t[region]*h[j]/mu_n[m];
                double numerator = psi_e[j][m] + source[j][m]*h[j]/mu_n[m]/2 + source_edge[j+1][m]*h[j]/mu_n[m]*tau*(0.25+tau/12);
                double denominator = 1 + tau + tau*tau/2 + tau*tau*tau/6;
                psi_e[j+1][m] = numerator/denominator;
                psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];
            }else if(alpha_mode==32){
                double tau = sigma_t[region]*h[j]/mu_n[m];
                double etau = exp(tau);
                double numerator = psi_e[j][m] + source[j][m]*h[j]/mu_n[m]/2 + source_edge[j+1][m]/sigma_t[region]/2*(etau-1-tau);
                psi_e[j+1][m] = numerator/etau;
                psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];
            }else if(alpha_mode==33){
                double tau = sigma_t[region]*h[j]/mu_n[m];
                double tau_a = (sigma_t[region]-sigma_s0[region])*h[j]/mu_n[m];
                double etau = exp(tau_a) - 1 - tau_a - tau_a*tau_a/2;
                double numerator = psi_e[j][m] + source[j][m]*h[j]/mu_n[m]/2 + source_edge[j+1][m]/sigma_t[region]*(tau*tau/4+etau/2);
                double denominator = 1 + tau + tau*tau/2 + etau;
                psi_e[j+1][m] = numerator/denominator;
                psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];
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

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//-----------LEFT AND RIGHT ITERATE------------------------------------------//
//-----------LEFT AND RIGHT ITERATE------------------------------------------//

//1UUUUUUU
//----------------------------------------------------------------------------//

//---------------UPDATE PHI, CALC SOURCE--------------------------------------//
//---------------UPDATE PHI, CALC SOURCE--------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

//Update phi, calculate source:
double SourceIteration::updatePhi_calcSource(bool usePsi){
    vector<double> old_phi_0(phi_0);
    vector<double> X = data->getX();
    X.insert(X.begin(),0);
    vector<double> Q = data->getQ();
    vector<double> Q_lin = data->getQ_lin();
    int region = 0;
    unsigned int within_region_counter = 0;
    
    if(usePsi){
       edgePhi0 = calcEdgePhi(0);
       edgePhi1 = calcEdgePhi(1);
    }
    for(unsigned int j = 0; j<J;j++){
        if(usePsi){ //If updating phi's with psi:
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
        //Need to calculate edge sources.  We are evaluating Q at the edge in order to preserve the linearity of the method.
        if(alpha_mode >=30){
            for(unsigned int m = 0; m < N/2;m++){
                source_edge[j+1][m] = sigma_s0[region]*edgePhi0[j+1]+3*mu_n[m]*sigma_s1[region]*edgePhi1[j+1] + Q[region] + Q_lin[region]*(x_e[j+1]);
            }
            for(unsigned int m = N/2; m<N;m++){
                source_edge[j][m] = sigma_s0[region]*edgePhi0[j]+3*mu_n[m]*sigma_s1[region]*edgePhi1[j] + Q[region] + Q_lin[region]*(x_e[j]);
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
            edgePhi0 = calcEdgePhi(0);
            edgePhi1 = calcEdgePhi(1);
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
    return (Utilities::p_norm(old_phi_0,phi_0,2));
    
}
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
    phi2_plus = phi_0 ;
    phi2_minus = phi_0;
    phi2e_plus = edgePhi0;
    phi2e_minus = edgePhi0;
    
    J = phi_0.size();
    bc = data->getbc();
    alpha_mode = data->getalpha_mode();
    if(alpha_mode == 35){
        edgePhi0_MB2_L = data->getedgePhi0_0();
        edgePhi0_MB2_R = data->getedgePhi0_0();
    }
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
        Qhat_edge.push_back(0);
    }
    initializeGrid();
    if( alpha_mode == 35 )
        initializew_n_MB2();
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
    
    kappa = Utilities::find_kappa(c,mu_n,w_n);
    
    rho = Utilities::zeros(J,1.0);
    if(alpha_mode >= 40)
        initializerho();
    
    setEdgePhi_toAvgPhi();
    updatePhi_calcSource(false);
    initializeDictionary();
}

SourceIteration::~SourceIteration(){
    data = NULL;
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//------------------------CONSTRUCTOR---------------------------------------//
//------------------------CONSTRUCTOR---------------------------------------//

//0IIIII
//--------------------------------------------------------------------------//

//-------------------------MAIN LOOP----------------------------------------//
//-------------------------MAIN LOOP----------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
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
//    double tol = fabs(data->gettol());
    double tol = fabs(data->gettol());
    if(!accel_mode){
        tol *= (1-c*falseConvCorrection);
    }
    double error = 0;
    
    //Iterate until tolerance is achieved:
    it_num = 0;
    do{
        old_phi_0 = phi_0;
        old_edgePhi0 = edgePhi0;
        vector<double> old_Qhat_edge(Qhat_edge);
        
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
        
//        Utilities::print_dmatrix(psi_e);
        
        variable_status["psi_e"] += 1;
        variable_status["psi_c"] += 1;
        
        updatePhi2();
        if((alpha_mode >= 40 && alpha_mode < 50) /*&& (error < tol*10 || Qhat_edge[0] != 0) */){
            updateQhat_edge();
        }
        
        
        error=updatePhi_calcSource();
        it_num += 0.5;
        
        //CMFD acceleration:
        if(accel_mode == 1){
            if(alpha_mode >= 40)
                accelerate_MB3();
            else
                cmfd();
            updatePhi_calcSource(false);
        }else if(accel_mode == 2 || accel_mode == 3){
            pcmfd();
            updatePhi_calcSource(false);
        }else{
//            variable_status["phi_0"] += 0.5;
//            variable_status["phi_1"] += 0.5;
//            variable_status["edgePhi0"] += 0.5;
//            variable_status["phi2_plus"] += 0.5;
//            variable_status["phi2_minus"] += 0.5;
//            variable_status["phi2e_plus"] += 0.5;
//            variable_status["phi2e_minus"] += 0.5;
//            variable_status["edgePhi1"] += 0.5;
//            variable_status["source"] += 0.5;
//            variable_status["source_edge"] += 0.5;
        }
        //This may or may not need to go inside the else statement above:
        //Depends on whether or not I feel like implementing CMFD with linear terms and whether they need to be accelerated.
        if(hasLinearTerms){
            variable_status["phi_0_lin"] += 0.5;
            variable_status["phi_1_lin"] += 0.5;
            variable_status["phi_c_lin"] += 0.5;
            variable_status["source_lin"] += 0.5;
        }
        
        it_num += 0.5;
        error = Utilities::p_norm_of_rel_error(old_phi_0,phi_0,2 , tol * tol);
//        error = Utilities::p_norm(old_phi_0,phi_0,2) / (Utilities::p_norm(phi_0,2) + tol*tol);
        // If we're using MB3, make sure the Qhats are converged too:
        if( alpha_mode >= 40 && alpha_mode < 50 && alpha_mode != 41 && accel_mode == 0){
            double Qhat_error = Utilities::p_norm_of_rel_error(old_Qhat_edge,Qhat_edge, 2 , tol);
//            Utilities::print_dvector(Qhat_edge);
            if( Utilities::p_norm(Qhat_edge,2) / (Utilities::p_norm(phi_0,2) + tol*tol) > tol*10)
                error = max( error , Qhat_error );
        }
        
        if(it_num == 1){
            spec_rad = 0;
        }else{
//            cout << "spec_rad = " << spec_rad << endl;
//            cout << "error/old_error = " error/old_error << endl; //ZZZZ
            double temp_spec_rad = pow( error/init_error, 1.0 / it_num );
            if (temp_spec_rad != 0)
                spec_rad = temp_spec_rad;//(spec_rad*(it_num-2)+error/old_error)/(it_num-1); //ZZZZ
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
            if(checkNegativeAngularFlux()){
                cout<<"WARNING: "<<checkNegativeAngularFlux()<<" negative angular flux values in iteration "<<it_num<<endl;
            }
            if(checkNegativeFlux()){
                cout<<"WARNING: "<<checkNegativeFlux()<<" negative flux values in iteration "<<it_num<<endl;
            }
        }
        old_error = error;
        if(it_num == 1){
            init_error = error;
            //If this is the first iteration, we skipped CMFD.  Continue at least one more iteration.
            continue;
        }
        
        
        //Check for stopping conditions:
        if(Utilities::nan_checker(phi_0)){ //NaN's in solution!
            cout<<"The flux has NaN's in it!"<<endl;
            break;
        }else if(error <= tol){ //if converged
            if(checkNegativeFlux()){ //Convergence to unphysical solution
                cout<<"Source iteration converged to a negative solution with " << checkNegativeFlux() << " negative scalar fluxes in "<<it_num<<" iterations"<<endl;
                it_num += 0.25;
            }else{  //Convergence to real solution
                cout<<"Source iteration converged in "<<it_num<<" iterations"<<endl;
                isConverged = 1;
            }
            break;
        }else if(it_num >= MAX_IT){
            //Max. # of iterations reached.  Did not converge.
            cout<<"Source iteration did NOT converge in "<<MAX_IT<<" iterations"<<endl;
            cout<<"The final error was "<<error<<endl;
            break;
        }else if(it_num >= MAX_IT_accel && accel_mode){
            //Max. # of iterations reached.  Did not converge.
            cout<<"Source iteration did NOT converge in "<<MAX_IT_accel<<" iterations"<<endl;
            cout<<"The final error was "<<error<<endl;
            break;
        }else if(it_num >= 5 && (error/init_error) >= diverge){
            //Solution has diverged significantly.
            cout<<"Source iteration diverged in "<<it_num<<" iterations"<<endl;
            break;
        }
        
    }while(true);
    //while(error>tol && it_num < MAX_IT && (it_num < MAX_IT_accel || accel_mode==0) && (it_num < 5 || (error/init_error)<diverge));
    if(isPrintingToWindow){
        outfile.close();
    }
    return 0;
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//-------------------MAIN LOOP--------------------------------------------//
//-------------------MAIN LOOP--------------------------------------------//

//0GGGGGGG
//------------------------------------------------------------------------//

//-------------------GET SOLUTION-----------------------------------------//
//-------------------GET SOLUTION-----------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
double SourceIteration::get_exitcurr(bool isExitOnRight){
    double exitcurr = 0;
    if(isExitOnRight){
        for(unsigned int m = 0 ; m < N/2 ; m++){
            exitcurr = psi_e[J][m]*w_n[m]*mu_n[m];
        }
    }else{
        for(unsigned int m = 0 ; m < N/2 ; m++){
            exitcurr = psi_e[J][m]*w_n[m]*mu_n[m];
        }
    }
    return exitcurr;
}

vector<vector<double> > SourceIteration::get_solution( bool oldSolution ){
    vector<vector<double> > out;
    
    //Compute the x coordinates of the phi vector:
    vector<double> x_total;
    for(unsigned int i = 0; i<x.size();i++){
        x_total.push_back(x_e[i]);
        x_total.push_back(x[i]);
    }
    x_total.push_back(x_e[x.size()]);
    
    out.push_back(x_total);
    
    vector<double> phi_total;
    if( oldSolution ){ // ZZZZ not sure if this makes sense if accel mode is on
        for(unsigned int i = 0; i<phi_0.size();i++){
            phi_total.push_back(old_edgePhi0[i]);
            phi_total.push_back(old_phi_0[i]);
        }
        phi_total.push_back(old_edgePhi0[phi_0.size()]);
    }else{
        vector<double> phi_edge;
        //Compute the phi vector (combining both edge and cell averaged)
        if(accel_mode != 1 || EDGE_ACCEL_MODE == 0){
            phi_edge = calcEdgePhi(0);
        }else{
            phi_edge = edgePhi0;
        }
        for(unsigned int i = 0; i<phi_0.size();i++){
            phi_total.push_back(phi_edge[i]);
            phi_total.push_back(phi_0[i]);
        }
        phi_total.push_back(phi_edge[phi_0.size()]);
    }
    
    out.push_back(phi_total);
    
    return out;
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//-------------------GET SOLUTION------------------------------------------//
//-------------------GET SOLUTION------------------------------------------//

//0PPPPPP
//------------------------------------------------------------------------//

//-------------------PRINT STUFF------------------------------------------//
//-------------------PRINT STUFF------------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

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
    outfile<<setw(5)<<"x"<<setw(tabwidth)<<"      "<<setw(tabwidth)<<"phi_0"<<setw(tabwidth)<<"      "<<setw(tabwidth)<<"phi_1"<<endl;
    if(isPrintingToWindow){
        cout<<setw(5)<<"x"<<setw(tabwidth)<<"      "<<setw(tabwidth)<<"phi_0"<<setw(tabwidth)<<"      "<<setw(tabwidth)<<"phi_1"<<endl;
    }
    for(unsigned int j = 0; j<J; j++){
        outfile<<setw(5)<<fixed<<x_e[j]<<setw(tabwidth)<<"      "<<setw(tabwidth)<<scientific<<phi_0e[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1e[j]<<endl;
        if(isPrintingToWindow){
            cout<<setw(5)<<fixed<<x_e[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<scientific<<phi_0e[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1e[j]<<endl;
        }
        outfile<<setw(5)<<fixed<<x[j]<<setw(tabwidth)<<"      "<<setw(tabwidth)<<scientific<<phi_0[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1[j]<<endl;
        if(isPrintingToWindow){
            cout<<setw(5)<<fixed<<x[j]<<setw(tabwidth)<<"      "<<setw(tabwidth)<<scientific<<phi_0[j]<<setw(tabwidth)<<"     "<<setw(tabwidth)<<phi_1[j]<<endl;
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

//0PPPPPPP
void SourceIteration::printDictionary(){
    for(Dict::const_iterator it(variable_status.begin());\
      it != variable_status.end(); it++){
        cout<<it->first<<" is at iteration "<<it->second<<endl;
    }
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//-------------------PRINT STUFF------------------------------------------//
//-------------------PRINT STUFF------------------------------------------//
//0PPPPP


//0CCCCC
//------------------------------------------------------------------------//

//-------------------CHECK FOR NEGATIVE FLUX-------------------------------//
//-------------------CHECK FOR NEGATIVE FLUX-------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
unsigned int SourceIteration::checkNegativeFlux(bool checkEdgeFluxes, bool printFluxes){
    unsigned int temp = 0;
    for(unsigned int j = 0; j<J;j++){
        if(phi_0[j] < 0){
            temp++;
            if(printFluxes)
                cout << "phi_0["<<j<<"] = " << phi_0[j] << endl;
        }
    }
    if(checkEdgeFluxes){
        for(unsigned int j = 0; j <= J; j++){
            if(edgePhi0[j] < 0){
                temp++;
                if(printFluxes)
                    cout << "edgePhi0[" << j << "] = " << edgePhi0[j] << endl;
            }
        }
    }
    return temp;
}


//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//-------------------CHECK FOR NEGATIVE FLUX-------------------------------//
//-------------------CHECK FOR NEGATIVE FLUX-------------------------------//

//0CCCCC
//------------------------------------------------------------------------//

//-------------------CHECK FOR NEGATIVE ANGULAR FLUX-------------------------------//
//-------------------CHECK FOR NEGATIVE ANGULAR FLUX-------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
unsigned int SourceIteration::checkNegativeAngularFlux(bool printFluxes){
    unsigned int temp = 0;
    for(unsigned int j = 0; j<J;j++){
        for(unsigned int m = 0; m < N; m++){
            if(psi_e[j][m] < 0){
                temp++;
                if(printFluxes)
                    cout << "psi_e[" << j << "]["<<m<<"] = " << psi_e[j][m] << endl;
            }
            if(psi_c[j][m] < 0){
                temp++;
                if(printFluxes)
                    cout << "psi_c[" << j << "]["<<m<<"] = " << psi_c[j][m] << endl;
            }
        }
    }
    for(unsigned int m = 0; m < N; m++){
        if(psi_e[J][m] < 0){
            temp++;
            if(printFluxes)
                cout << "psi_e[" << J << "]["<<m<<"] = " << psi_e[J][m] << endl;
        }
    }
    return temp;
}


//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//-------------------CHECK FOR NEGATIVE ANGULAR FLUX--------------------------//
//-------------------CHECK FOR NEGATIVE ANGULAR FLUX--------------------------//

//0AAAA
//------------------------------------------------------------------------//

//-------------------ACCELERATE EDGE SCALAR FLUX------------------------------//
//-------------------ACCELERATE EDGE SCALAR FLUX------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
void SourceIteration::accelerate_edgePhi0(vector<double> preaccel_phi_0){
    vector<double> Q = data->getQ();
    
    unsigned int region_L = 0;
    unsigned int region_R = 0;
    
    //Left edge:
    double left_coeff = 0;
    double right_coeff = 0;
    double Q_addition = 0;
    double phi_addition = 0;
    
    double Sth_L;
    double Sth_R = sigma_t[region_R]*h[0];
    double Q_left_edge = h[0]*Q[region_R];
    for(unsigned int m = 0; m< N/2; m++){
        left_coeff += psi_e[0][m] * w_n[m];
        right_coeff += 2 * fabs(mu_n[N-m-1]) * psi_c[0][N-m-1] * \
            w_n[N-m-1] / ( Sth_R + 2 * fabs(mu_n[N-m-1]) );
        
        if(alpha_mode >= 40 && alpha_mode < 50 && bc[0] != 1){
            Q_addition += (Q_left_edge + h[0]*Qhat_edge[0]*mu_n[N-m-1])* w_n[N-m-1] / ( Sth_R + 2 * fabs(mu_n[N-m-1]) )/2;
        }else{
            Q_addition += Q_left_edge * w_n[N-m-1] / ( Sth_R + 2 * fabs(mu_n[N-m-1]) )/2;
        }
        phi_addition += h[0]*sigma_s0[region_R] * w_n[N-m-1]/ 2 / \
        ( Sth_R + 2 * fabs(mu_n[N-m-1]) );
    }
    right_coeff /= preaccel_phi_0[0];
    
    if(bc[0] != 1){
        edgePhi0[0] = left_coeff+right_coeff*phi_0[0] + Q_addition;
        edgePhi0[0] /= 1 - phi_addition;
    }else{
        edgePhi0[0] = 2*(right_coeff*phi_0[0] + Q_addition);
        edgePhi0[0] /= 1 - 2*phi_addition;
    }
    
    unsigned int within_region_counter_L = 0;
    unsigned int within_region_counter_R = 1;
    if(within_region_counter_L == discret[region_L]){
        within_region_counter_L=0;
        region_L++;
    }
    if(within_region_counter_R == discret[region_R]){
        within_region_counter_R=0;
        region_R++;
    }
    
    for(unsigned int i = 1; i<edgePhi0.size()-1;i++){
        left_coeff = 0;
        right_coeff = 0;
        Q_addition = 0;
        phi_addition = 0;
        
        Sth_L = sigma_t[region_L]*h[i-1];
        Sth_R = sigma_t[region_R]*h[i];
        
        double Q_L = h[i-1]*Q[region_L];
        double Q_R = h[i]*Q[region_R];
        
        for(unsigned int m = 0; m< N/2; m++){
            left_coeff += 2 * fabs(mu_n[m]) * psi_c[i-1][m] * w_n[m] / \
            ( Sth_L + 2 * fabs(mu_n[m]) ) ;
            right_coeff += 2 * fabs(mu_n[N-m-1]) * psi_c[i][N-m-1] * \
            w_n[N-m-1] / ( Sth_R + 2 * fabs(mu_n[N-m-1]) );
            
            if(alpha_mode >=40 && alpha_mode < 50){
                Q_addition += (Q_L + h[i-1]*Qhat_edge[i] * mu_n[m])* w_n[m] / \
                ( Sth_L + 2 * fabs(mu_n[m]) )/2;
                Q_addition += (Q_R + h[i]*Qhat_edge[i] * mu_n[N-m-1])* w_n[N-m-1] / \
                ( Sth_R + 2 * fabs(mu_n[N-m-1]) )/2;
            }else{
                Q_addition += Q_L * w_n[m] / \
                    ( Sth_L + 2 * fabs(mu_n[m]) )/2;
                Q_addition += Q_R * w_n[N-m-1] / \
                    ( Sth_R + 2 * fabs(mu_n[N-m-1]) )/2;
            }
            
            phi_addition += h[i-1]*sigma_s0[region_L] * w_n[m]/ 2 / \
            ( Sth_L + 2 * mu_n[m] );
            phi_addition += h[i]*sigma_s0[region_R] * w_n[N-m-1]/ 2 / \
            ( Sth_R + 2 * fabs(mu_n[N-m-1]) );
        }
        left_coeff /= preaccel_phi_0[i-1];
        right_coeff /= preaccel_phi_0[i];
        
        edgePhi0[i] = left_coeff*phi_0[i-1]+right_coeff*phi_0[i] + Q_addition;
        edgePhi0[i] /= 1 - phi_addition;
        
        within_region_counter_L++;
        within_region_counter_R++;
        if(within_region_counter_L == discret[region_L]){
            within_region_counter_L=0;
            region_L++;
        }
        if(within_region_counter_R == discret[region_R]){
            within_region_counter_R=0;
            region_R++;
        }
        
    }
    
    //Right edge:
    left_coeff = 0;
    right_coeff = 0;
    Q_addition = 0;
    phi_addition = 0;
    unsigned int i = h.size();
    
    Sth_L = sigma_t[region_L]*h[i-1];
    double Q_right_edge = h[i-1]*Q[region_L];
    for(unsigned int m = 0; m< N/2; m++){
        left_coeff += 2 * mu_n[m] * psi_c[i-1][m] * w_n[m] / \
        ( Sth_L + 2 * mu_n[m] ) ;
        right_coeff += psi_e[i][N-m-1];
        
        if(alpha_mode >= 30 && alpha_mode < 50 && bc[1] != 1){
            Q_addition += (Q_right_edge + h[i-1]*Qhat_edge[J]*mu_n[m]) * w_n[m] /\
            ( Sth_L + 2 * fabs(mu_n[m]) ) / 2;
        }else{
            Q_addition += Q_right_edge * w_n[m] /\
            ( Sth_L + 2 * fabs(mu_n[m]) ) / 2;
        }
        phi_addition += h[i-1]*sigma_s0[region_L] * w_n[m]/ 2 / \
        ( Sth_L + 2 * mu_n[m] );
    }
    
    left_coeff /= preaccel_phi_0[i-1];
    
    if(bc[1]!=1){
        edgePhi0[i] = left_coeff*phi_0[i-1]+right_coeff + Q_addition;
        edgePhi0[i] /= 1 - phi_addition;
    }else{
        edgePhi0[i] = 2*(left_coeff*phi_0[i-1] + Q_addition);
        edgePhi0[i] /= 1 - 2*phi_addition;
    }
    

}
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//-------------------ACCELERATE EDGE SCALAR FLUX------------------------------//
//-------------------ACCELERATE EDGE SCALAR FLUX------------------------------//


//0AAAA
//------------------------------------------------------------------------//

//-------------------ACCELERATE EDGE SCALAR FLUX FOR MB2----------------------//
//-------------------ACCELERATE EDGE SCALAR FLUX FOR MB2----------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
void SourceIteration::accelerate_edgePhi0_MB2(vector<double> preaccel_phi_0){
    //NOTE THAT accelerate_edgePhi0 needs to be run IN ADDITION to this
    //  function for the MB2 case
    
    vector<double> Q = data->getQ();
    
    unsigned int region_L = 0;
    unsigned int region_R = 0;
    
    //Left edge:
    double left_coeff = 0;
    double right_coeff = 0;
    double Q_addition = 0;
    double phi_addition = 0;
    
    double Sth_L;
    double Sth_R = sigma_t[region_R]*h[0];
    for(unsigned int m = 0; m< N/2; m++){
        left_coeff += psi_e[0][m] * w_n_MB2[0][m];
        right_coeff += 2 * fabs(mu_n[N-m-1]) * psi_c[0][N-m-1] * \
        w_n_MB2[0][N-m-1] / ( Sth_R + 2 * fabs(mu_n[N-m-1]) );
        
        Q_addition += w_n_MB2[0][N-m-1] / ( Sth_R + 2 * fabs(mu_n[N-m-1]) )/2;
        phi_addition += h[0]*sigma_s0[region_R] * w_n_MB2[0][N-m-1]/ 2 / \
        ( Sth_R + 2 * fabs(mu_n[N-m-1]) );
    }
    right_coeff /= preaccel_phi_0[0];
    Q_addition *= h[0]*Q[region_R];
    
    if(bc[0] != 1){
        edgePhi0_MB2_L[0] = left_coeff+right_coeff*phi_0[0] + Q_addition;
        edgePhi0_MB2_L[0] /= 1 - phi_addition;
    }else{
        edgePhi0_MB2_L[0] = 2*(right_coeff*phi_0[0] + Q_addition);
        edgePhi0_MB2_L[0] /= 1 - 2*phi_addition;
    }
    
    unsigned int within_region_counter_L = 0;
    unsigned int within_region_counter_R = 1;
    if(within_region_counter_L == discret[region_L]){
        within_region_counter_L=0;
        region_L++;
    }
    if(within_region_counter_R == discret[region_R]){
        within_region_counter_R=0;
        region_R++;
    }
    
    //edgePhi0_MB2_L:
    for(unsigned int i = 1; i<edgePhi0.size()-1;i++){
        left_coeff = 0;
        right_coeff = 0;
        Q_addition = 0;
        phi_addition = 0;
        
        Sth_L = sigma_t[region_L]*h[i-1];
        Sth_R = sigma_t[region_R]*h[i];
        
        double Q_L = h[i-1]*Q[region_L];
        double Q_R = h[i]*Q[region_R];
        
        for(unsigned int m = 0; m< N/2; m++){
            left_coeff += 2 * fabs(mu_n[m]) * psi_c[i-1][m] * w_n_MB2[i][m] / \
            ( Sth_L + 2 * fabs(mu_n[m]) ) ;
            right_coeff += 2 * fabs(mu_n[N-m-1]) * psi_c[i][N-m-1] * \
            w_n_MB2[i][N-m-1] / ( Sth_R + 2 * fabs(mu_n[N-m-1]) );
            
            Q_addition += Q_L * w_n_MB2[i][m] / \
            ( Sth_L + 2 * fabs(mu_n[m]) )/2;
            Q_addition += Q_R * w_n_MB2[i][N-m-1] / \
            ( Sth_R + 2 * fabs(mu_n[N-m-1]) )/2;
            
            phi_addition += h[i-1]*sigma_s0[region_L] * w_n_MB2[i][m]/ 2 / \
            ( Sth_L + 2 * mu_n[m] );
            phi_addition += h[i]*sigma_s0[region_R] * w_n_MB2[i][N-m-1]/ 2 / \
            ( Sth_R + 2 * fabs(mu_n[N-m-1]) );
        }
        left_coeff /= preaccel_phi_0[i-1];
        right_coeff /= preaccel_phi_0[i];
        
        edgePhi0_MB2_L[i] = left_coeff*phi_0[i-1]+right_coeff*phi_0[i] + \
        Q_addition;
        edgePhi0_MB2_L[i] /= 1 - phi_addition;
        
        within_region_counter_L++;
        within_region_counter_R++;
        if(within_region_counter_L == discret[region_L]){
            within_region_counter_L=0;
            region_L++;
        }
        if(within_region_counter_R == discret[region_R]){
            within_region_counter_R=0;
            region_R++;
        }
        
    }
    within_region_counter_L = 0;
    within_region_counter_R = 1;
    region_L = 0;
    region_R = 0;
    if(within_region_counter_L == discret[region_L]){
        within_region_counter_L=0;
        region_L++;
    }
    if(within_region_counter_R == discret[region_R]){
        within_region_counter_R=0;
        region_R++;
    }
    
    //edgePhi0_MB2_R:
    for(unsigned int i = 1; i<edgePhi0.size()-1;i++){
        left_coeff = 0;
        right_coeff = 0;
        Q_addition = 0;
        phi_addition = 0;
        
        Sth_L = sigma_t[region_L]*h[i-1];
        Sth_R = sigma_t[region_R]*h[i];
        
        double Q_L = h[i-1]*Q[region_L];
        double Q_R = h[i]*Q[region_R];
        
        for(unsigned int m = 0; m< N/2; m++){
            left_coeff += 2 * fabs(mu_n[m]) * psi_c[i-1][m] * w_n_MB2[i-1][m] / \
            ( Sth_L + 2 * fabs(mu_n[m]) ) ;
            right_coeff += 2 * fabs(mu_n[N-m-1]) * psi_c[i][N-m-1] * \
            w_n_MB2[i-1][N-m-1] / ( Sth_R + 2 * fabs(mu_n[N-m-1]) );
            
            Q_addition += Q_L * w_n_MB2[i-1][m] / \
            ( Sth_L + 2 * fabs(mu_n[m]) )/2;
            Q_addition += Q_R * w_n_MB2[i-1][N-m-1] / \
            ( Sth_R + 2 * fabs(mu_n[N-m-1]) )/2;
            
            phi_addition += h[i-1]*sigma_s0[region_L] * w_n_MB2[i-1][m]/ 2 / \
            ( Sth_L + 2 * mu_n[m] );
            phi_addition += h[i]*sigma_s0[region_R] * w_n_MB2[i-1][N-m-1]/ 2 / \
            ( Sth_R + 2 * fabs(mu_n[N-m-1]) );
        }
        left_coeff /= preaccel_phi_0[i-1];
        right_coeff /= preaccel_phi_0[i];
        
        edgePhi0_MB2_R[i] = left_coeff*phi_0[i-1]+right_coeff*phi_0[i] + \
        Q_addition;
        edgePhi0_MB2_R[i] /= 1 - phi_addition;
        
        within_region_counter_L++;
        within_region_counter_R++;
        if(within_region_counter_L == discret[region_L]){
            within_region_counter_L=0;
            region_L++;
        }
        if(within_region_counter_R == discret[region_R]){
            within_region_counter_R=0;
            region_R++;
        }
        
    }
    
    //Right edge:
    left_coeff = 0;
    right_coeff = 0;
    Q_addition = 0;
    phi_addition = 0;
    unsigned int i = h.size();
    
    Sth_L = sigma_t[region_L]*h[i-1];
    
    for(unsigned int m = 0; m< N/2; m++){
        left_coeff += 2 * mu_n[m] * psi_c[i-1][m] * w_n_MB2[i-1][m] / \
        ( Sth_L + 2 * mu_n[m] ) ;
        right_coeff += psi_e[i][N-m-1];
        
        Q_addition += w_n_MB2[i-1][m] / ( Sth_L + 2 * fabs(mu_n[m]) ) / 2;
        phi_addition += h[i-1]*sigma_s0[region_L] * w_n_MB2[i-1][m]/ 2 / \
        ( Sth_L + 2 * mu_n[m] );
    }
    
    left_coeff /= preaccel_phi_0[i-1];
    Q_addition *= h[i-1]*Q[region_L];
    
    if(bc[1]!=1){
        edgePhi0_MB2_R[i] = left_coeff*phi_0[i-1]+right_coeff + Q_addition;
        edgePhi0_MB2_R[i] /= 1 - phi_addition;
    }else{
        edgePhi0_MB2_R[i] = 2*(left_coeff*phi_0[i-1] + Q_addition);
        edgePhi0_MB2_R[i] /= 1 - 2*phi_addition;
    }
    
}
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//-------------------ACCELERATE EDGE SCALAR FLUX FOR MB2----------------------//
//-------------------ACCELERATE EDGE SCALAR FLUX FOR MB2----------------------//


//1AAAAA
//-------------------ACCELERATION FOR MB3-------------------------------------//
//-------------------ACCELERATION FOR MB3-------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
void SourceIteration::accelerate_MB3(){
    
    vector<double> sigma_a = data->getsigma_a();
    
    vector<double> mu_n2(mu_n);
    for(unsigned int m = 0; m < N; m++)
        mu_n2[m] *= mu_n2[m];
    
    vector<double> phi_2(J,0);
    for(unsigned int j = 0; j < J; j++)
        phi_2[j] = phi2_plus[j] + phi2_minus[j];
    
    vector<double> eddington(phi_0);
    for(unsigned int j = 0; j < J ; j++)
        eddington[j] = phi_2[j]/phi_0[j];
    
    vector< vector<double> > A( J , vector<double>( 3 , 0 ) );
    vector<double> b( J , 0 );
    
    unsigned int region = 0;
    unsigned int within_region_counter = 1;
    vector<double> Q = data->getQ();
    
    
    //Deal with the left edge:
    double Sth2 = sigma_t[0] * h[0] * h[0]; //ZZZZ only works for homogeneous right now
    double Sth = sigma_t[0] * h[0];
    
    //Calculate a bunch of parameters needed to eliminate phi_{0,1/2} in terms of phi_{0,1}
    double phi_2e_L= 0;
    for(unsigned int m = 0; m < N ; m++ )
        phi_2e_L += mu_n2[m]*psi_e[0][m]*w_n[m];
    
    /*if(EDGE_ACCEL_MODE == 2){
        double A01temp = 0;
        double denom = 0;
        for(unsigned int m = N/2; m < N; m++)
            denom -= w_n[m] / ( Sth - 2 * mu_n[m] );
        denom *= sigma_s0[0] * h[0];
        denom /= 2;
        denom += 1;
        for( unsigned int m = N/2 ; m < N ; m++ )
            A01temp -= mu_n[m] * w_n[m] * psi_c[0][m] / ( Sth - 2 * mu_n[m] );
        A01temp *= 2;
        A01temp /= denom * phi_0[0];
        
        double in_flux = 0;
        for(unsigned int m = 0; m < N/2 ; m++ )
            in_flux += w_n[m] * psi_e[0][m];
        double Q_addition = 0;
        for(unsigned int m = N/2 ; m < N ; m++ ){
            Q_addition += ( Q[0] + mu_n[m] * Qhat_edge[0] ) * w_n[m] / \
                (Sth - 2 * mu_n[m]);
        }
        Q_addition *= h[0] / 2;
        
        if(bc[0] == 1)
            A[0][1] = (sigma_t[0]-sigma_s0[0]) * Sth2 + eddington[0];
        else
            A[0][1] = (sigma_t[0]-sigma_s0[0]) * Sth2 + 3 * eddington[0] - 2 * phi_2e_L/edgePhi0[0] * A01temp;
        A[0][2] = - eddington[1];
        b[0] = Q[0] * Sth2;
        if(bc[0] != 1 )
            b[0] += 2 * phi_2e_L / edgePhi0[0] * (in_flux + Q_addition) / denom;
    }else  //ZZZZ not functional at the moment */ if (EDGE_ACCEL_MODE == 1){
        A[0][1] = sigma_a[0] * h[0] + eddington[0] / Sth_avg[0] + 2 * ( eddington[0]-  phi_2e_L/phi_0[0] ) / Sth;
        A[0][2] = - eddington[1] / Sth_avg[0];
        b[0] = Q[0] * h[0];
    }else if (EDGE_ACCEL_MODE == 3){
        double eddington_half =phi_2e_L/edgePhi0[0];
        A[0][1] =  sigma_a[region] * h[0] + (2 * eddington[0] - eddington_half) / Sth + eddington[0]/Sth_avg[0];
        A[0][2] = - eddington[1] / Sth_avg[0];
        b[0] = Q[0] * h[0] +  (2 * phi_2e_L - eddington_half * phi_0[0])/Sth;
    }
    
    
    //Deal with the center of problem:
    for(unsigned int j = 1; j < J-1 ; j++){
        if( within_region_counter == discret[region] ){
            region++;
            within_region_counter = 0;
            Sth2 = sigma_t[region] * h[j] * h[j];
        }
        
        A[j][0] = -eddington[j-1] / Sth_avg[j-1] ;
        A[j][1] = sigma_a[region] * h[j] + eddington[j] / Sth_avg[j-1] + eddington[j]/Sth_avg[j];
        A[j][2] = -eddington[j+1] / Sth_avg[j];
        
        b[j] = Q[region] * h[j];
        
        within_region_counter++;
    }
    
    //Deal with the right edge:
    region = sigma_t.size()-1;
    unsigned int j = J-1;
    Sth2 = sigma_t[region] * h[j] * h[j];
    Sth = sigma_t[region] * h[j];
    
    double phi_2e_R= 0;
    for(unsigned int m = 0; m < N ; m++ )
        phi_2e_R += mu_n2[m]*psi_e[J][m]*w_n[m];
    
    /*if( EDGE_ACCEL_MODE == 2 ){
        double AJ1temp = 0;
        double denom = 0;
        for(unsigned int m = 0; m < N/2; m++)
            denom -= w_n[m] / ( Sth + 2 * mu_n[m] );
        denom *= sigma_s0[region] * h[j];
        denom /= 2;
        denom += 1;
        for( unsigned int m = 0 ; m < N/2 ; m++ )
            AJ1temp += mu_n[m] * w_n[m] * psi_c[j][m] / ( Sth + 2 * mu_n[m] );
        AJ1temp *= 2;
        AJ1temp /= denom * phi_0[j];
        
        double out_flux = 0;
        for(unsigned int m = N/2 ; m < N ; m++ )
            out_flux += w_n[m] * psi_e[J][m];
        double Q_addition = 0;
        for(unsigned int m = 0; m < N/2 ; m++ ){
            Q_addition += ( Q[region] + mu_n[m] * Qhat_edge[J] ) * w_n[m] / \
                (Sth + 2 * mu_n[m]);
        }
        Q_addition *= h[j] / 2;
        
        A[j][0] = -eddington[j-1];
        if( bc[1] == 1 )
           A[j][1] = (sigma_t[region] - sigma_s0[region]) * Sth2 + eddington[j];
        else
            A[j][1] = (sigma_t[region] - sigma_s0[region]) * Sth2 + 3 * eddington[j] - 2 * phi_2e_R / edgePhi0[J] * AJ1temp;
        b[j] = Q[region] * Sth2;
        if( bc[1] != 1 )
               b[j] += 2 * phi_2e_R / edgePhi0[J] * (out_flux + Q_addition) / denom;
    }else  //ZZZZ */ if( EDGE_ACCEL_MODE == 1 ){
        A[j][0] = -eddington[j-1] / Sth_avg[j-1];
        A[j][1] = sigma_a[region] * h[j] + eddington[j] / Sth_avg[j-1] + 2 * (eddington[j] - phi_2e_R / phi_0[j] )/ Sth;
        b[j] = Q[region] * h[j];
    }else if( EDGE_ACCEL_MODE == 3 ){
        double eddington_half =phi_2e_R/edgePhi0[J];
        A[j][1] = - eddington[j-1] / Sth_avg[j-1];
        A[j][2] =  sigma_a[region] * h[j] + (2 * eddington[j] - eddington_half) / Sth + eddington[j]/Sth_avg[j-1];
        b[0] = Q[region] * h[j] +  (2 * phi_2e_R - eddington_half * phi_0[j])/Sth;
    }
    
    vector<double> preaccel_phi_0(phi_0);
    vector<double> preaccel_edgePhi0(edgePhi0);
    phi_0 = Utilities::solve_tridiag(A,b);
    variable_status["phi_0"] += 0.5;
    variable_status["phi_1"] += 0.5;
    
    /* //ZZZZ
    if( EDGE_ACCEL_MODE == 2 )
        accelerate_edgePhi0( preaccel_phi_0 );
    else*/ if (EDGE_ACCEL_MODE == 1 ){
        edgePhi0[0] *= phi_0[0]/preaccel_phi_0[0];
        for(unsigned int j = 1; j < J; j++){
            edgePhi0[j] *= (phi_0[j-1] + phi_0[j])/(preaccel_phi_0[j-1] + preaccel_phi_0[j]);
        }
        edgePhi0[J] *= phi_0[J-1]/preaccel_phi_0[J-1];
    }else if (EDGE_ACCEL_MODE == 3){
        edgePhi0[0] += phi_0[0]-preaccel_phi_0[0];
        for(unsigned int j = 1; j < J; j++){
            edgePhi0[j] += (phi_0[j-1] + phi_0[j])/2-(preaccel_phi_0[j-1] + preaccel_phi_0[j])/2;
        }
        edgePhi0[J] += phi_0[J-1]-preaccel_phi_0[J-1];
    }
    variable_status["edgePhi0"] += 0.5;
    variable_status["edgePhi1"] += 0.5;
    
    
    /*//Accelerate phi_2:
    Utilities::print_dvector(phi2_plus); // ZZZZ
    for(unsigned int j = 0; j < J ; j++){
        double ratio = phi_0[j] / preaccel_phi_0[j];
        phi2_plus[j] *= ratio;
        phi2_minus[j] *= ratio;
        ratio = edgePhi0[j] / preaccel_edgePhi0[j];
        phi2e_minus[j] *= ratio;
        if(bc[0] == 1 && j == 0)
            phi2e_plus[j] = phi2e_minus[j];
        else if( j != 0 )
            phi2e_plus[j] *= ratio;
    }
    phi2e_plus[J] *= edgePhi0[J] / preaccel_edgePhi0[J];
    if(bc[1] == 1)
        phi2e_minus[J] = phi2e_plus[J];
    Utilities::print_dvector(phi2_plus); // ZZZZ
    cout << "===" <<endl;
    
    variable_status["phi2_plus"] += 0.5;
    variable_status["phi2e_minus"] += 0.5;
    variable_status["phi2_plus"] += 0.5;
    variable_status["phi2e_minus"] += 0.5;

//    Utilities::print_dvector(Qhat_edge); //ZZZZ
    updateQhat_edge();
//    Utilities::print_dvector(Qhat_edge); //ZZZZ
//    cout << "===" << endl;*/
    
    /*
    Utilities::print_dvector(Qhat_edge);
    Qhat_edge[0] *= phi_0[0]/preaccel_phi_0[0];
    for(unsigned int j = 1; j < J; j++){
        Qhat_edge[j] *= (phi_0[j-1] + phi_0[j])/(preaccel_phi_0[j-1] + preaccel_phi_0[j]);
    }
    Qhat_edge[J] *= phi_0[J-1]/preaccel_phi_0[J-1];
    Utilities::print_dvector(Qhat_edge);
    cout << "+===" << endl;*/ //ZZZZ
    
}
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//-------------------ACCELERATION FOR MB3-------------------------------------//
//-------------------ACCELERATION FOR MB3-------------------------------------//

//1BBBBBB
//1CCCCCC
//------------------------------------------------------------------------//

//-------------------CALCULATE EDGE PHI------------------------------------//
//-------------------CALCULATE EDGE PHI------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
vector<double> SourceIteration::calcEdgePhi(int num){
    vector<double> edgePhi;
    if(!num){
        for(unsigned int j = 0; j<J+1;j++){
            edgePhi.push_back(0);
            for(unsigned int m = 0; m<N;m++){
                edgePhi[j]+= w_n[m]*psi_e[j][m];
            }
        }
        variable_status["edgePhi0"] = variable_status["psi_e"];
    }else{
        for(unsigned int j = 0; j<J+1;j++){
            edgePhi.push_back(0);
            for( unsigned int m = 0; m < N ; m++ ){
                edgePhi[j]+= mu_n[m]*w_n[m]*psi_e[j][m];
            }
        }
        variable_status["edgePhi1"] = variable_status["psi_e"];
    }
    return edgePhi;
    
}
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//-------------------CALCULATE EDGE PHI------------------------------------//
//-------------------CALCULATE EDGE PHI------------------------------------//

//1CCCCCC
//------------------------------------------------------------------------//

//-------------------CMFD/PCMFD--------------------------------------------//
//-------------------CMFD/PCMFD--------------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
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
    for(unsigned int i = 1; i<phi_0_CM_size-1;i++){
        phi_0_CM_new[i] = Q_CM[i];
    }
    phi_0_CM_new[phi_0_CM_size-1] = Q_CM[phi_0_CM_size-1] + 2*out_curr;
    
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
    vector<double> preaccel_phi_0(phi_0);
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
    
    if( (alpha_mode >=30 && alpha_mode < 50) || alpha_mode==11){
        if(EDGE_ACCEL_MODE==1){
            accelerate_edgePhi0(preaccel_phi_0);
            if(alpha_mode == 35)
                accelerate_edgePhi0_MB2(preaccel_phi_0);
        }
    }
    
    variable_status["edgePhi0"] += 0.5;
    if(alpha_mode==35)
        variable_status["edgePhi0_MB2"] += 0.5;
    variable_status["edgePhi1"] += 0.5;
    variable_status["phi_0"] += 0.5;
    variable_status["phi_1"] += 0.5;
    
}

//1CCCCCC
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
    vector<double> preaccel_phi_0(phi_0);
    
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
    
    variable_status["edgePhi0"] += 0.5;
    if(alpha_mode==35)
        variable_status["edgePhi0_MB2"] += 0.5;
    variable_status["edgePhi1"] += 0.5;
    variable_status["phi_0"] += 0.5;
    variable_status["phi_1"] += 0.5;
    
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//-------------------CMFD/PCMFD--------------------------------------------//
//-------------------CMFD/PCMFD--------------------------------------------//


//1FFFFFFF
//------------------------------------------------------------------------//

//----------------FINITE DIFFERENCE------------------------------------------//
//----------------FINITE DIFFERENCE------------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//


void SourceIteration::finiteDifference(){
    for(unsigned int j = 0; j<J;j++){
        for(unsigned int m = 0; m<N;m++){
            psi_c[j][m] = ((1.0+alpha[j][m])*psi_e[j+1][m]+(1.0-alpha[j][m])*psi_e[j][m])/2.0;
        }
    }
    return;
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//-------------------FINITE DIFFERENCE-------------------------------------//
//-------------------FINITE DIFFERENCE-------------------------------------//


//1IIIIII
//------------------------------------------------------------------------//

//----------------INITIALIZE STUFF------------------------------------------//
//----------------INITIALIZE STUFF------------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

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

//1IIIIII
void SourceIteration::initializeDictionary(){
    variable_status["edgePhi0"] = 0;
    if(alpha_mode == 35)
        variable_status["edgePhi0_MB2"] = 0;
    variable_status["edgePhi1"] = 0;
    variable_status["phi_0"] = 0;
    variable_status["phi2_plus"] = 0;
    variable_status["phi2_minus"] = 0;
    variable_status["phi2e_plus"] = 0;
    variable_status["phi2e_minus"] = 0;
    variable_status["phi_1"] = 0;
    variable_status["psi_e"] = -0.5;
    variable_status["psi_c"] = -0.5;
    variable_status["source"] = 0;
    if(alpha_mode >=30 && alpha_mode < 50){
        variable_status["source_edge"] = 0;
    }else{
        variable_status["source_edge"] = -1;
    }
    if(alpha_mode >= 40 && alpha_mode < 50){
        variable_status["Qhat_edge"] = 0;
    }
    if(hasLinearTerms){
        variable_status["phi_0_lin"] = 0;
        variable_status["phi_1_lin"] = 0;
        variable_status["phi_c_lin"] = 0;
        variable_status["source_lin"] = 0;
    }else{
        variable_status["phi_0_lin"] = -1;
        variable_status["phi_1_lin"] = -1;
        variable_status["phi_c_lin"] = -1;
        variable_status["source_lin"] = -1;
    }
}


//1IIIIII
void SourceIteration::initializerho(){
    //Temporary value for rho.  This needs to be made a vector for non-uniform meshes.
    if( c == 1 || alpha_mode != 40 ){
        return; //rho is already 1, we want rho to stay at 1
    }
    unsigned int counter = 0;
    unsigned int region = 0;
    for(unsigned int j = 0; j < J; j++){
        if(counter == discret[region]){
            region++;
            counter = 0;
        }
        
        double term = sigma_t[region]*kappa*h[j];
        if (fabs(term) < 0.01) {
            double term2 = term*term;
            double term4 = term2*term2;
            rho[j] = 1. + term2/12 + term4/360 + term4*term2/20160; //Use Taylor expansion if term is too small
        }else{
            rho[j] = exp(-term) - 2. + exp(term);
            rho[j] /= (term*term);
        }
        
        counter++;
    }
}


//1IIIIII
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
    
    unsigned int within_region_counter_L = 0;
    unsigned int within_region_counter_R = 1;
    unsigned int region_L = 0;
    unsigned int region_R = 0;
    if(within_region_counter_R == discret[region_R]){
        region_R++;
        within_region_counter_R = 0;
    }
    //Initialize h_avg:
    for(unsigned int i = 0; i < h.size()-1; i++){
        h_avg.push_back( ( h[i] + h[i+1] ) / 2 );
        
        Sth_avg.push_back( ( sigma_t[region_L]*h[i] + sigma_t[region_R]*h[i+1] )/2 );
        
        within_region_counter_L++;
        within_region_counter_R++;
        if(within_region_counter_L == discret[region_L]){
            region_L++;
            within_region_counter_L = 0;
        }
        if(within_region_counter_R == discret[region_R]){
            region_R++;
            within_region_counter_R = 0;
        }
    }
    
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

//IIII
void SourceIteration::initializew_n_MB2(){
    
    unsigned int region = 0;
    unsigned int counter = 0;
    
    double denom = 0;
    for(unsigned int j = 0; j < J; j++){
        if(counter == discret[region]){
            region++;
            counter = 0;
        }
        if( counter == 0 ){
            denom = 0;
            for(unsigned int m = 0; m < N; m++){
                denom += fabs(mu_n[m]) * w_n[m] / ( fabs(mu_n[m]) + sigma_t[region]* h[j] / 2);
            }
        }
        
        vector<double> temp(N,0);
        for(unsigned int m = 0; m < N; m++){
            temp[m] = 2*fabs(mu_n[m]) * w_n[m] / ( fabs(mu_n[m]) + sigma_t[region]*h[j] / 2) / denom;
        }
        w_n_MB2.push_back(temp);
        
        counter++;
    }
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//----------------INITIALIZE STUFF------------------------------------------//
//----------------INITIALIZE STUFF------------------------------------------//


//1IIIIIIII
//------------------------------------------------------------------------//

//-----------LEFT AND RIGHT ITERATE----------------------------------------//
//-----------LEFT AND RIGHT ITERATE----------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

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
            }else if(alpha_mode ==30 || alpha_mode ==35){
                double tau = -sigma_t[region]*h[j]/mu_n[m];
                double numerator = psi_e[j+1][m] - source[j][m]*h[j]/mu_n[m]/2 - source_edge[j][m]*h[j]*tau/4/mu_n[m];
                double denominator = 1 + tau + tau*tau/2;
                psi_e[j][m] = numerator/denominator;
                psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];
            }else if(alpha_mode==31){
                double tau = -sigma_t[region]*h[j]/mu_n[m];
                double numerator = psi_e[j+1][m] - source[j][m]*h[j]/mu_n[m]/2 - source_edge[j][m]*h[j]/mu_n[m]*tau*(0.25+tau/12);
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
            }else if(alpha_mode==40  || alpha_mode == 42){
                double tau = sigma_t[region]*h[j]/mu_n[m];
                double numerator = psi_e[j+1][m] - source[j][m]*tau/2/sigma_t[region] + ( source_edge[j][m] + mu_n[m]*Qhat_edge[j] )*rho[j]*tau*tau/4./sigma_t[region];
                double denominator = 1 - tau + rho[j]*tau*tau/2;
                psi_e[j][m] = numerator/denominator;
                psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];
            }else if(alpha_mode==41){
                double tau = -sigma_t[region]*h[j]/mu_n[m];
                double numerator = psi_e[j+1][m] - source[j][m]*h[j]/mu_n[m]/2 - (source_edge[j][m] - mu_n[m]*Qhat_edge[j]*2)*h[j]*rho[j]*tau/4/mu_n[m];
                double denominator = 1 + tau + rho[j]*tau*tau/2;
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
            }else if(alpha_mode ==30 || alpha_mode ==35){
                double tau = sigma_t[region]*h[j]/mu_n[m];
                double numerator = psi_e[j][m] + source[j][m]*h[j]/mu_n[m]/2 + source_edge[j+1][m]*h[j]/mu_n[m]*tau/4;
                double denominator = 1 + tau + tau*tau/2;
                psi_e[j+1][m] = numerator/denominator;
                psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];
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
            }else if(alpha_mode==40 || alpha_mode == 42){
                double tau = sigma_t[region]*h[j]/mu_n[m];
                double numerator = psi_e[j][m] + source[j][m]*h[j]/mu_n[m]/2. + (source_edge[j+1][m]+mu_n[m]*Qhat_edge[j+1])*rho[j]*h[j]/mu_n[m]*tau/4.;
                double denominator = 1 + tau + rho[j]*tau*tau/2.;
                psi_e[j+1][m] = numerator/denominator;
                psi_c[j][m] = (source[j][m]/2- mu_n[m]/h[j]*(psi_e[j+1][m]-psi_e[j][m]))/sigma_t[region];
            }else if(alpha_mode==41){
                double tau = sigma_t[region]*h[j]/mu_n[m];
                double numerator = psi_e[j][m] + source[j][m]*h[j]/mu_n[m]/2 + (source_edge[j+1][m]-mu_n[m]*Qhat_edge[j+1]*2)*rho[j]*h[j]/mu_n[m]*tau/4;
                double denominator = 1 + tau + rho[j]*tau*tau/2;
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

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//-----------LEFT AND RIGHT ITERATE----------------------------------------//
//-----------LEFT AND RIGHT ITERATE----------------------------------------//

//1SSSSSSS
//------------------------------------------------------------------------//

//----------------CALC EDGE PHI AS AVG OF NEIGHBORING PHI--------------------//
//----------------CALC EDGE PHI AS AVG OF NEIGHBORING PHI--------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

void SourceIteration::setEdgePhi_toAvgPhi(){
    unsigned int phisize = phi_0.size();
    edgePhi0[0] = max((double)0.0, \
        phi_0[0] - (phi_0[1] - phi_0[0])/(x[1]-x[0])*(x[0]-x_e[0]));
    for(unsigned int i = 1; i < phisize;i++){
        edgePhi0[i] = (phi_0[i-1] + phi_0[i])/2;
    }
    edgePhi0[phisize] = max((double)0.0,phi_0[phisize-1] + (phi_0[phisize-1]-\
        phi_0[phisize-2]) / (x[phisize-1]-x[phisize-2]) * \
        (x_e[phisize]-x[phisize-1]));
}


//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//----------------CALC EDGE PHI AS AVG OF NEIGHBORING PHI--------------------//
//----------------CALC EDGE PHI AS AVG OF NEIGHBORING PHI--------------------//


//1UUUUUUU
//------------------------------------------------------------------------//

//---------------UPDATE PHI, CALC SOURCE------------------------------------//
//---------------UPDATE PHI, CALC SOURCE------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//


void SourceIteration::updateEdgePhi0_MB2(){
    for(unsigned int j = 0; j < J+1; j++){
        if( j == 0 ){
            edgePhi0_MB2_R[j] = 0;
            for( unsigned int m = 0; m < N; m++ ){
                edgePhi0_MB2_R[j] += psi_e[j][m] * w_n_MB2[j][m];
            }
        } else{
            edgePhi0_MB2_R[j] = 0;
            for( unsigned int m = 0; m < N; m++ ){
                edgePhi0_MB2_R[j] += psi_e[j][m] * w_n_MB2[j-1][m];
            }
        }
        
        if(j == 0){
            edgePhi0_MB2_L[0] = edgePhi0_MB2_R[0];
        }else if( j < J ){
            edgePhi0_MB2_L[j] = 0;
            for( unsigned int m = 0; m < N; m++ ){
                edgePhi0_MB2_L[j] += psi_e[j][m] * w_n_MB2[j][m];
            }
        }else{
            edgePhi0_MB2_L[J] = edgePhi0_MB2_R[J];
        }
        
    }
//    cout << "Left:";
//    Utilities::print_dvector(edgePhi0_MB2_L);
//    cout << "Right:";
//    Utilities::print_dvector(edgePhi0_MB2_R);
//    Utilities::print_dvector(edgePhi0);
    variable_status["edgePhi0_MB2"] = variable_status["psi_e"];
}


//Update phi, calculate source:
double SourceIteration::updatePhi_calcSource(bool usePsi){
    //Make copies:
    vector<double> loc_old_phi_0(phi_0);
    
    vector<double> X = data->getX();
    X.insert(X.begin(),0);
    vector<double> Q = data->getQ();
    vector<double> Q_lin = data->getQ_lin();
    int region = 0;
    unsigned int within_region_counter = 0;
    
    if(usePsi){
       edgePhi0 = calcEdgePhi(0);
       edgePhi1 = calcEdgePhi(1);
       if(alpha_mode == 35)
           updateEdgePhi0_MB2();
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
        if(alpha_mode >=30 && alpha_mode < 50){
            if(alpha_mode == 35){
                for(unsigned int m = 0; m < N/2;m++){
                    source_edge[j+1][m] = sigma_s0[region]*edgePhi0_MB2_R[j+1]+3*mu_n[m]*sigma_s1[region]*edgePhi1[j+1] + Q[region] + Q_lin[region]*(x_e[j+1]);
                }
                for(unsigned int m = N/2; m<N;m++){
                    source_edge[j][m] = sigma_s0[region]*edgePhi0_MB2_L[j]+3*mu_n[m]*sigma_s1[region]*edgePhi1[j] + Q[region] + Q_lin[region]*(x_e[j]);
                }
            }else{
                for(unsigned int m = 0; m < N/2;m++){
                    source_edge[j+1][m] = sigma_s0[region]*edgePhi0[j+1]+3*mu_n[m]*sigma_s1[region]*edgePhi1[j+1] + Q[region] + Q_lin[region]*(x_e[j+1]);
                }
                for(unsigned int m = N/2; m<N;m++){
                    source_edge[j][m] = sigma_s0[region]*edgePhi0[j]+3*mu_n[m]*sigma_s1[region]*edgePhi1[j] + Q[region] + Q_lin[region]*(x_e[j]);
                }
            }
        }
        within_region_counter++;
        if(within_region_counter==discret[region]){
            within_region_counter = 0;
            region++;
        }
    }
    if(!hasLinearTerms){
        variable_status["source"] += 0.5;
        if(usePsi){
            variable_status["phi_0"] = variable_status["psi_c"];
            variable_status["phi_1"] = variable_status["psi_c"];
        }
    }
    if(alpha_mode >=30 && alpha_mode < 50){
        variable_status["source_edge"] += 0.5;
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
            variable_status["source"] += 0.5;
            if(usePsi){
                variable_status["phi_0"] = variable_status["psi_c"];
                variable_status["phi_1"] = variable_status["psi_c"];
            }
            variable_status["source_lin"] += 0.5;
        }else{
            region = 0;
            within_region_counter = 0;
            edgePhi0 = calcEdgePhi(0);
            edgePhi1 = calcEdgePhi(1);
            if( alpha_mode==35 )
                updateEdgePhi0_MB2();
            for(unsigned int j = 0; j<J;j++){
                if(edgePhi0[j+1]+edgePhi0[j] == 0){
                    phi_0_lin[j] = 0;
                }else{
                    phi_0_lin[j] = 2*phi_0[j]/h[j]*(edgePhi0[j+1]-edgePhi0[j])/(edgePhi0[j+1]+edgePhi0[j]); //ZZZZ need to fix for MB-2?
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
            variable_status["source"] += 0.5;
            variable_status["source_lin"] += 0.5;
        }
    }
    return (Utilities::p_norm(loc_old_phi_0,phi_0,2));
    
}
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//---------------UPDATE PHI, CALC SOURCE------------------------------------//
//---------------UPDATE PHI, CALC SOURCE------------------------------------//

//1UUUUUUU
//------------------------------------------------------------------------//

//---------------UPDATE PHI2----------------------------------------------//
//---------------UPDATE PHI2----------------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
void SourceIteration::updatePhi2(){
    unsigned int lastind = phi2_plus.size();
    vector<double> mu_n2(mu_n.size(),0);
    for(unsigned int m = 0; m < N ; m++){
        mu_n2[m] = mu_n[m]*mu_n[m];
    }
    
    for(unsigned int j = 0; j < lastind; j++ ){
        phi2_plus[j] = 0;
        phi2_minus[j] = 0;
        phi2e_plus[j] = 0;
        phi2e_minus[j] = 0;
        for(unsigned int m = 0; m < N/2 ; m++){
            phi2_plus[j] += mu_n2[m]*w_n[m]*psi_c[j][m];
            phi2e_plus[j] += mu_n2[m]*w_n[m]*psi_e[j][m];
        }
        for(unsigned int m = N/2; m < N ; m++){
            phi2_minus[j] += mu_n2[m]*w_n[m]*psi_c[j][m];
            phi2e_minus[j] += mu_n2[m]*w_n[m]*psi_e[j][m];
        }
    }
    phi2e_plus[lastind] = 0;
    phi2e_minus[lastind] = 0;
    for(unsigned int m = 0; m < N/2 ; m++){
        phi2e_plus[lastind] += mu_n2[m]*w_n[m]*psi_e[lastind][m];
    }
    for(unsigned int m = N/2; m < N ; m++){
        phi2e_minus[lastind] += mu_n2[m]*w_n[m]*psi_e[lastind][m];
    }
    variable_status["phi2_plus"] = variable_status["psi_c"];
    variable_status["phi2_minus"] = variable_status["psi_c"];
    variable_status["phi2e_plus"] = variable_status["psi_e"];
    variable_status["phi2e_minus"] = variable_status["psi_e"];

}
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
//---------------UPDATE PHI2----------------------------------------------//
//---------------UPDATE PHI2----------------------------------------------//

//1UUUUUUU
//------------------------------------------------------------------------//

//---------------UPDATE QHAT_EDGE-----------------------------------------//
//---------------UPDATE QHAT_EDGE-----------------------------------------//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

void SourceIteration::updateQhat_edge(){
    unsigned int lastind = Qhat_edge.size()-1;
    if(alpha_mode != 45){
        double int_mu_0_to_1 = 0; //Integral of mu_n from 0 to 1 via quadrature
        double int_mu_sq_0_to_1 = 0; //Integral of mu_n^2 from 0 to 1 via quadrature
        for(unsigned int m = 0; m < N/2; m ++){
            int_mu_0_to_1 += mu_n[m]*w_n[m];
            int_mu_sq_0_to_1 += mu_n[m]*mu_n[m]*w_n[m];
        }
        
        vector<double> Q = data->getQ();
        Qhat_edge[0] = 0;
        if( bc[0] != 1 ){
            for(unsigned int m = 0; m < N/2; m ++){
                Qhat_edge[0] += mu_n[m]*psi_e[0][m]*w_n[m];
            }
            Qhat_edge[0] *= -6.*sigma_t[0];
            Qhat_edge[0] += 3*int_mu_0_to_1*(sigma_s0[0]*edgePhi0[0]+Q[0]);
            Qhat_edge[0] -= 12./rho[0]/h[0]*(phi2_plus[0]-phi2e_plus[0]);
        }
        for(unsigned int j = 1; j < lastind; j++){
            Qhat_edge[j] = 2 * ( phi2e_plus[j] - phi2_plus[j-1] ) / h[j-1];
            Qhat_edge[j] -= ( phi2_plus[j] - phi2_plus[j-1] + phi2_minus[j] - phi2_minus[j-1] ) / h_avg[j-1];
            Qhat_edge[j] += 2 * ( phi2_minus[j] - phi2e_minus[j] ) / h[j];
            Qhat_edge[j] *= 3 / rho[j];
        }
        
        unsigned int last_region = sigma_t.size()-1;
        Qhat_edge[lastind] = 0;
        if( bc[1] != 1 ){
            for(unsigned int m = N/2; m < N ; m++){
                Qhat_edge[lastind] += mu_n[m]*psi_e[lastind][m]*w_n[m];
            }
            Qhat_edge[lastind] *= -6.*sigma_t[last_region];
            Qhat_edge[lastind] -= 3*int_mu_0_to_1*(sigma_s0[last_region]*edgePhi0[lastind]+Q[last_region]);
            Qhat_edge[lastind] -= 12./rho[J-1]/h[J-1]*(phi2e_minus[lastind]-phi2_minus[lastind-1]);
        }
        
    }else{
        //ZZZZ THIS ONLY WORKS FOR ONE ZONE UNIFORM MESH RIGHT NOW
        double term = - 9. / ( 4 * rho[0] * h[0] );
        
        double J_inc_L = 0; //left incoming current
        for(unsigned int m = 0; m < N/2; m++){
            J_inc_L += abs(mu_n[m]) * w_n[m] * psi_e[0][m];
        }
        double J_out_L = 0; //left outgoing current
        for(unsigned int m = N/2; m < N; m++){
            J_out_L += abs(mu_n[m]) * w_n[m] * psi_e[0][m];
        }
        Qhat_edge[0] = 0;//-term * ( phi_1[0] - edgePhi1[0] ) + 3 * sigma_t[0] * (J_inc_L + J_out_L)/2 - 3./4 * sigma_s0[0] * phi_0[0];
        
        cout << "First part = " <<  -term * ( phi_1[0] - edgePhi1[0] ) << endl;
        cout << "Second part = " << 3 * sigma_t[0] * (J_inc_L + J_out_L)/2 - 3./4 * sigma_s0[0] * phi_0[0] << endl;
        
        for(unsigned int i = 1; i < lastind; i++){
            Qhat_edge[i] = edgePhi1[i] - (phi_1[i-1]+phi_1[i])/2;
            Qhat_edge[i] *= term;
        }
        
        double J_inc_R = 0; //right incoming current
        for( unsigned int m = N/2; m < N; m++ ){
            J_inc_R += abs(mu_n[m]) * w_n[m] * psi_e[lastind][m];
        }
        double J_out_R = 0; //right incoming current
        for( unsigned int m = 0; m < N/2; m++ ){
            J_out_R += abs(mu_n[m]) * w_n[m] * psi_e[lastind][m];
        }
        Qhat_edge[lastind] = 0;//term * ( edgePhi1[lastind] - phi_1[lastind-1] ) - 3 * sigma_t[0] * (J_inc_R + J_out_R)/2 + 3./4 * sigma_s0[0] * phi_0[lastind-1];
    }
    
    variable_status["Qhat_edge"] = variable_status["phi2_plus"];
//    cout << "|Qhat_edge| = " << Utilities::p_norm(Qhat_edge,2) << endl; //ZZZZ
//    cout << "Qhat_edge[0] = " << Qhat_edge[0] << endl;
//    cout << "Qhat_edge[phi_0.size()] = " << Qhat_edge[phi_0.size()] << endl;
//    cout << "Qhat_edge[phi_0.size()/2] = " << Qhat_edge[phi_0.size()/2] << endl;
//    Utilities::print_dvector(Qhat_edge);
    
}
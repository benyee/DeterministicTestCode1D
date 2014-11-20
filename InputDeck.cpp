//
//  InputDeck.cpp
//  
//
//  Created by Ben Yee on 9/20/13.
//
//

#include "InputDeck.h"

InputDeck::InputDeck(string inputName){
    fileName = inputName;
}

int InputDeck::loadInputDeck(){
    //Read in the file:
    cout<<"Reading input from "<<fileName<<endl;
    
    ifstream inputFile;
    inputFile.open(fileName.c_str());
    string line;
    
    //Find the X variable:
    if(!searchForInput(inputFile,"X")){
        return 1;
    }
    while(getline(inputFile,line)){
    //Check if we have read in all the values yet:
        if(line == "."){
            break;
        }
    //Otherwise tack on the new value:
        X.push_back(atof(line.c_str()));
    }
    //-------------------------------------------
    if(!searchForInput(inputFile,"discret")){
        return 1;
    }
    while(getline(inputFile,line)){
        if(line == "." || discret.size() == X.size()){
            break;
        }
        discret.push_back(atoi(line.c_str()));
    }
    //-------------------------------------------
    if(!searchForInput(inputFile,"sigma_s0")){
        return 1;
    }
    while(getline(inputFile,line)){
        if(line == "." || sigma_s0.size() == X.size()){
            break;
        }
        sigma_s0.push_back(atof(line.c_str()));
    }
    //If only one entry is given, assume the value is constant for the entire domain
    while( sigma_s0.size() < discret.size() )
        sigma_s0.push_back(sigma_s0[0]);
    //-------------------------------------------
    if(!searchForInput(inputFile,"sigma_s1")){
        return 1;
    }
    while(getline(inputFile,line)){
        if(line == "." || sigma_s1.size() == X.size()){
            break;
        }
        sigma_s1.push_back(atof(line.c_str()));
    }
    //If only one entry is given, assume the value is constant for the entire domain
    while( sigma_s1.size() < discret.size() )
        sigma_s1.push_back(sigma_s1[0]);
    //-------------------------------------------
    if(!searchForInput(inputFile,"sigma_a")){
        return 1;
    }
    while(getline(inputFile,line)){
        if(line == "." || sigma_a.size() == X.size()){
            break;
        }
        sigma_a.push_back(atof(line.c_str()));
    }
    //If only one entry is given, assume the value is constant for the entire domain
    while( sigma_a.size() < discret.size() )
        sigma_a.push_back(sigma_a[0]);
    //-------------------------------------------
    if(!searchForInput(inputFile,"Q")){
        return 1;
    }
    while(getline(inputFile,line)){
        if(line == "." || Q.size() == X.size()){
            break;
        }
        Q.push_back(atof(line.c_str()));
    }
    //If only one entry is given, assume the value is constant for the entire domain
    while( Q.size() < discret.size() )
        Q.push_back(Q[0]);
    //-------------------------------------------
    if(!searchForInput(inputFile,"Q_lin")){
        return 1;
    }
    while(getline(inputFile,line)){
        if(line == "." || Q_lin.size() == X.size()){
            break;
        }
        Q_lin.push_back(atof(line.c_str()));
    }
    //If only one entry is given, assume the value is constant for the entire domain
    while( Q_lin.size() < discret.size() )
        Q_lin.push_back(Q_lin[0]);
    //--------------------------------------------
    if(!searchForInput(inputFile,"N")){
        return 1;
    }
    getline(inputFile,line);
    N = atoi(line.c_str());
    //-------------------------------------------
    if(!searchForInput(inputFile,"bc")){
        return 1;
    }
    getline(inputFile,line);
    bc[0] = atoi(line.c_str());
    getline(inputFile,line);
    bc[1] = atoi(line.c_str());
    if(bc[0]==2){
        if(!searchForInput(inputFile,"left")){
            return 1;
        }
        while(getline(inputFile,line)){
            if(line == "." || psi_bl.size() >= N/2){
                break;
            }
            psi_bl.push_back(atof(line.c_str()));
        }
    }
    if(bc[1]==2){
        if(!searchForInput(inputFile,"right")){
            return 1;
        }
        while(getline(inputFile,line)){
            if(line == "." || psi_br.size() >= N/2){
                break;
            }
            psi_br.push_back(atof(line.c_str()));
        }
    }
    //--------------------------------------------
    if(!searchForInput(inputFile,"alpha_mode")){
        return 1;
    }
    getline(inputFile,line);
    alpha_mode = atoi(line.c_str());
    //--------------------------------------------
    if(!searchForInput(inputFile,"tol")){
        return 1;
    }
    getline(inputFile,line);
    tol = atof(line.c_str());
    //--------------------------------------------
    //Get initial phi_0:
    if(!searchForInput(inputFile,"phi_0_0")){
        return 1;
    }
    getline(inputFile,line);
    if(line=="default"){
        for(unsigned int j = 0; j < discret.size(); j++){
            for(unsigned int i = 0; i < discret[j]; i++){
                phi_0_0.push_back(0);
            }
        }
    } else{
        phi_0_0.push_back(atof(line.c_str()));
        while(getline(inputFile,line)){
            if(line == "."){
                break;
            }
            phi_0_0.push_back(atof(line.c_str()));
        }
    }
    
    //--------------------------------------------
    //Get initial phi_1:
    if(!searchForInput(inputFile,"phi_1_0")){
        return 1;
    }
    getline(inputFile,line);
    if(line=="default"){
        for(unsigned int j = 0; j < discret.size(); j++){
            for(unsigned int i = 0; i < discret[j]; i++){
                phi_1_0.push_back(0);
            }
        }
    } else{
        phi_1_0.push_back(atof(line.c_str()));
        while(getline(inputFile,line)){
            if(line == "."){
                break;
            }
            phi_1_0.push_back(atof(line.c_str()));
        }
    }
    if(alpha_mode >=10 && alpha_mode < 30){
        hasLinearTerms = 1;
    }
    else{
        hasLinearTerms = 0;
    }
    //-------------------------------------------
    //Get initial phi_0_lin and phi_1_lin:
    if(hasLinearTerms){
        //--------------------------------------------
        //Get initial phi_0_lin:
        if(!searchForInput(inputFile,"phi_0_0_lin")){
            return 1;
        }
        getline(inputFile,line);
        if(line=="default"){
            for(unsigned int j = 0; j < discret.size(); j++){
                for(unsigned int i = 0; i < discret[j]; i++){
                    phi_0_0_lin.push_back(0);
                }
            }
        } else{
            phi_0_0_lin.push_back(atof(line.c_str()));
            while(getline(inputFile,line)){
                if(line == "."){
                    break;
                }
                phi_0_0_lin.push_back(atof(line.c_str()));
            }
        }
        
        //--------------------------------------------
        //Get initial phi_1:
        if(!searchForInput(inputFile,"phi_1_0_lin")){
            return 1;
        }
        getline(inputFile,line);
        if(line=="default"){
            for(unsigned int j = 0; j < discret.size(); j++){
                for(unsigned int i = 0; i < discret[j]; i++){
                    phi_1_0_lin.push_back(0);
                }
            }
        } else{
            phi_1_0_lin.push_back(atof(line.c_str()));
            while(getline(inputFile,line)){
                if(line == "."){
                    break;
                }
                phi_1_0_lin.push_back(atof(line.c_str()));
            }
        }
    }
    if(!searchForInput(inputFile,"accel_mode",false)){
        accel_mode = 0;
        discret_CM = discret;
    }else{
        getline(inputFile,line);
        accel_mode = atoi(line.c_str());
        if(!searchForInput(inputFile,"discret_CM")){
            cout<<"Problem with discret_CM... setting discret_CM = discret"<<endl;
            discret_CM = discret;
        }else{
            unsigned int i = 0;
            while(getline(inputFile,line)){
                if(line == "." || discret_CM.size() == X.size()){
                    break;
                }
                discret_CM.push_back(atoi(line.c_str()));
                if(discret[i] % discret_CM[i] != 0){
                    cout<<"WARNING: discret_CM["<<i<<"] does not divide evenly into discret["<<i<<"].  Setting discret_CM = discret instead."<<endl;
                    discret_CM[i] = discret[i];
                }
                i = i+1;
            }
        }
    }

    
    inputFile.close();
    
    //Initialize edge fluxes:
    vector<double> temp(phi_0_0.size()+1,0.0);
    edgePhi0_0 = temp;
    vector<double> temp2(phi_0_0.size()+1,0.0);
    edgePhi1_0 = temp2;
    
    //Check to make sure all these vectors are the same size:
    int vectorSizes[] = {X.size(), discret.size(),sigma_s0.size(),sigma_s1.size(),sigma_a.size(), Q.size(),Q_lin.size()};
    for(unsigned int i = 0; i < 6; i++){
        if(vectorSizes[i] != vectorSizes[i+1]){
            std::cout<<"There's an issue with the material definition vector sizes!"<<endl;
            return 1;
        }
    }
    
    //Make sure the initial values is the proper size:
    unsigned int expectedSize = 0;
    for(unsigned int i = 0; i<X.size(); i++){
        expectedSize += discret[i];
    }
    if(expectedSize != phi_0_0.size()){
        std::cout<<"There's an issue with the vector sizes of the initial conditions!"<<endl;
        return 1;
    }else if(expectedSize!=phi_1_0.size()){
        std::cout<<"There's an issue with the vector sizes of the initial conditions!"<<endl;
        return 1;
    }
    
    
    //Repeat if there are linear terms:
    if(hasLinearTerms){
        if(expectedSize != phi_0_0_lin.size()){
            std::cout<<"There's an issue with the vector sizes of the linear part of the initial conditions!"<<endl;
            return 1;
        }else if(expectedSize!=phi_1_0_lin.size()){
            std::cout<<"There's an issue with the vector sizes of the linear part of the initial conditions!"<<endl;
            return 1;
        }
    }
    
    return 0;
}

void InputDeck::readValues(){
    cout<<"X = [0 ";
    for(unsigned int i = 0; i < X.size(); i++){
        cout<<X[i]<<" ";
    }
    cout<<"]"<<endl;
    
    cout<<"discret = ";
    Utilities::print_uivector(discret);
    if(accel_mode){
        cout<<"discret_CM = ";
        Utilities::print_uivector(discret_CM);
    }
    cout<<"sigma_s0 = ";
    Utilities::print_dvector(sigma_s0);
    cout<<"sigma_s1 = ";
    Utilities::print_dvector(sigma_s1);
    cout<<"sigma_a = ";
    Utilities::print_dvector(sigma_a);
    cout<<"Q = ";
    Utilities::print_dvector(Q);
    cout<<"Q_lin = ";
    Utilities::print_dvector(Q_lin);
    
    cout<<"N = "<<N<<endl;
    cout<<"bc = ["<<bc[0]<<' '<<bc[1]<<']'<<endl;
    if(bc[0] ==2){
        cout<<"psi_bl = ";
        Utilities::print_dvector(psi_bl);
    }
    if(bc[1]==2){
        cout<<"psi_br = ";
        Utilities::print_dvector(psi_br);
    }
    cout<<"alpha_mode = "<<alpha_mode<<endl;
    cout<<"accel_mode = "<<accel_mode<<endl;
    cout<<"tol = "<<tol<<endl;
    
    cout<<"phi_0_0 = [";
    unsigned int temp = 0;
    unsigned int temp2 = 0;
    for(unsigned int i = 0; i < phi_0_0.size(); i++){
        //Divide the phi_0_0 readout by material sections
        if(temp < discret.size() && temp2 == (unsigned int)discret[temp] && i!= phi_0_0.size()-1){
            cout<<endl;
            temp++;
            temp2 = 0;
        }
        cout<<phi_0_0[i]<<" ";
        temp2++;
    }
    
    cout<<"]"<<endl;
    cout<<"phi_1_0 = [";
    temp = 0;
    temp2 = 0;
    for(unsigned int i = 0; i < phi_1_0.size(); i++){
        //Divide the phi_1_0 readout by material sections
        if(temp < discret.size() && temp2 == (unsigned int)discret[temp] && i!= phi_1_0.size()-1){
            cout<<endl;
            temp++;
            temp2 = 0;
        }
        cout<<phi_1_0[i]<<" ";
        temp2++;
    }
    cout<<"]"<<endl;
    
    if(hasLinearTerms){
        cout<<"]"<<endl;
        cout<<"phi_0_0_lin = [";
        unsigned int temp = 0;
        unsigned int temp2 = 0;
        for(unsigned int i = 0; i < phi_0_0_lin.size(); i++){
            //Divide the phi_0_0_lin readout by material sections
            if(temp < discret.size() && temp2 == (unsigned int)discret[temp] && i!= phi_0_0_lin.size()-1){
                cout<<endl;
                temp++;
                temp2 = 0;
            }
            cout<<phi_0_0_lin[i]<<" ";
            temp2++;
        }
        
        cout<<"]"<<endl;
        cout<<"phi_1_0_lin = [";
        temp = 0;
        temp2 = 0;
        for(unsigned int i = 0; i < phi_1_0_lin.size(); i++){
            //Divide the phi_1_0_lin readout by material sections
            if(temp < discret.size() && temp2 == (unsigned int)discret[temp] && i!= phi_1_0_lin.size()-1){
                cout<<endl;
                temp++;
                temp2 = 0;
            }
            cout<<phi_1_0_lin[i]<<" ";
            temp2++;
        }
        cout<<"]"<<endl;
    }
}

bool InputDeck::searchForInput(ifstream &file, string inp, bool hasErrorOutput){
    string line;
    while(getline(file,line)){
        if(line == "<--"+inp+"-->"){
            return 1;
        }
    }
    if(hasErrorOutput){
        cout<<"We couldn't find "<<inp<<"!"<<endl;
    }
    return 0;
}

void InputDeck::diffusionSolve(){
    if(discret.size() > 1){
        cout<<"InputDeck::diffusionSolve() doesn't work for multi-region problems yet."<<endl;
        return;
    }
    
    vector<double> temp(discret[0]+1,0);
    vector<double> temp2(discret[0]+1,0);
    
    edgePhi0_0 = temp;
    edgePhi1_0 = temp2;
    
    double dx = X[0]/discret[0];
    double D = 1.0/3/(sigma_s0[0]-sigma_s1[0]+sigma_a[0]);
    vector<double> mu_n = Utilities::calc_mu_n(N);
    vector<double> w_n = Utilities::calc_w_n(mu_n);
    
    vector<vector<double> > A(2*discret[0]+1,vector<double>(3,0));
    vector<double> b(2*discret[0]+1,Q[0]);
    b[0] = 0;
    if(bc[0] == 1){
        A[0][1] = 1;
        A[0][2] = -1;
    }else{
        A[0][1] = (1 + 4*D/dx);
        A[0][2] = -4*D/dx;
        if(bc[0] == 2){
            for(unsigned int m = 0; m < psi_bl.size();m++){
                b[0] += mu_n[m]*w_n[m]*psi_bl[m];
            }
            b[0] *=4;
        }
    }
    
    for(unsigned int i = 1; i<A.size()-1;i++){
        A[i][0] = -4*D/dx/dx;
        A[i][1] = (8*D/dx/dx + sigma_a[0]);
        A[i][2] = -4*D/dx/dx;
    }
    
    unsigned int bsizem1 = b.size()-1;
    b[bsizem1] = 0;
    if(bc[1] == 1){
        A[bsizem1][0] = -1;
        A[bsizem1][1] = 1;
    }else{
        A[bsizem1][0] = -4*D/dx;
        A[bsizem1][1] = (1 + 4*D/dx);
        if(bc[1]==2){
            for(unsigned int m = N/2; m < N;m++){
                b[bsizem1] -= mu_n[m]*w_n[m]*psi_br[m-N/2];
            }
            b[bsizem1] *=4;
        }
    }
//    Utilities::print_dmatrix(A);
//    Utilities::print_dvector(b);
    
    vector<double> temp_phi = Utilities::solve_tridiag(A,b);
//    Utilities::print_dvector(temp_phi); //ZZZZ
    
    //Make sure nothing is negative:
    for(unsigned int i = 0; i<temp_phi.size();i++){
        if(temp_phi[i] < 0){
            cout << "temp_phi["<<i<<"] = " << temp_phi[i] << endl;
            temp_phi[i] == 0;
        }
    }
//    Utilities::print_dvector(temp_phi);
    Utilities::split_Phi(temp_phi,edgePhi0_0,phi_0_0);
    vector<double> temp_phi_1(temp_phi);
    
    //Estimate initial currents:
    temp_phi_1[0] = 2*D*(temp_phi[0]-temp_phi[1])/dx;
    for(unsigned int i = 1; i<temp_phi.size()-1;i++){
        temp_phi_1[i] = (temp_phi[i-1]-temp_phi[i+1])/dx*D;
    }
    temp_phi_1[temp_phi.size()-1] = 2*(temp_phi[temp_phi.size()-2] - temp_phi[temp_phi.size()-1])*D/dx;
    Utilities::split_Phi(temp_phi_1,edgePhi1_0,phi_1_0);
    
    cout<<"End of diffusion solve!"<<endl;
//    Utilities::print_dvector(temp_phi);
//    Utilities::print_dvector(phi_0_0);
//    Utilities::print_dvector(edgePhi0_0);
}


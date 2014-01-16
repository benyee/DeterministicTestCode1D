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
    if(alpha_mode >=10){
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
    }else{
        getline(inputFile,line);
        accel_mode = atoi(line.c_str());
    }

    
    inputFile.close();
    
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


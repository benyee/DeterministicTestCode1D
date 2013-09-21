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
    if(!searchForInput(inputFile,"bc")){
        return 1;
    }
    getline(inputFile,line);
    bc[0] = atoi(line.c_str());
    getline(inputFile,line);
    bc[1] = atoi(line.c_str());
    //--------------------------------------------
    if(!searchForInput(inputFile,"N")){
        return 1;
    }
    getline(inputFile,line);
    N = atoi(line.c_str());
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
            for(int i = 0; i < discret[j]; i++){
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
            for(int i = 0; i < discret[j]; i++){
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
    
    inputFile.close();
    
    //Check to make sure all these vectors are the same size:
    int vectorSizes[] = {X.size(), discret.size(),sigma_s0.size(),sigma_s1.size(),sigma_a.size(), Q.size()};
    for(unsigned int i = 0; i < 5; i++){
        if(vectorSizes[i] != vectorSizes[i+1]){
            return 1;
        }
    }
    
    //Make sure the initial values is the proper size:
    unsigned int expectedSize = 0;
    for(unsigned int i = 0; i<X.size(); i++){
        expectedSize += discret[i];
    }
    if(expectedSize != phi_0_0.size()){
        return 1;
    }else if(expectedSize!=phi_1_0.size()){
        return 1;
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
    Utilities::print_ivector(discret);
    cout<<"sigma_s0 = ";
    Utilities::print_dvector(sigma_s0);
    cout<<"sigma_s1 = ";
    Utilities::print_dvector(sigma_s1);
    cout<<"sigma_a = ";
    Utilities::print_dvector(sigma_a);
    cout<<"Q = ";
    Utilities::print_dvector(Q);
    
    cout<<"bc = ["<<bc[0]<<' '<<bc[1]<<']'<<endl;
    cout<<"N = "<<N<<endl;
    cout<<"alpha_mode = "<<alpha_mode<<endl;
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
        //Divide the phi_0_0 readout by material sections
        if(temp < discret.size() && temp2 == (unsigned int)discret[temp] && i!= phi_1_0.size()-1){
            cout<<endl;
            temp++;
            temp2 = 0;
        }
        cout<<phi_1_0[i]<<" ";
        temp2++;
    }
    cout<<"]"<<endl;
}

bool InputDeck::searchForInput(ifstream &file, string inp){
    string line;
    while(getline(file,line)){
        if(line == "<--"+inp+"-->"){
            return 1;
        }
    }
    cout<<"We couldn't find it!"<<endl;
    return 0;
}


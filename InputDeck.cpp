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
    if(!searchForInput(inputFile,"phi_0_0")){
        return 1;
    }
    getline(inputFile,line);
    if(line=="default"){
        for(int j = 0; j < discret.size(); j++){
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
    //
    
    inputFile.close();
    
    //Check to make sure all these vectors are the same size:
    int vectorSizes[] = {X.size(), discret.size(),sigma_s0.size(),sigma_s1.size(),sigma_a.size(), Q.size()};
    for(int i = 0; i < 5; i++){
        if(vectorSizes[i] != vectorSizes[i+1]){
            return 1;
        }
    }
    
    //Make sure the initial values is the proper size:
    int expectedSize = 0;
    for(int i = 0; i<X.size(); i++){
        expectedSize += discret[i];
    }
    if(expectedSize != phi_0_0.size()){
        return 1;
    }
    
    return 0;
}

void InputDeck::readValues(){
    cout<<"X = [0 ";
    for(int i = 0; i < X.size(); i++){
        cout<<X[i]<<" ";
    }
    cout<<"]"<<endl;
    
    cout<<"discret = [ ";
    for(int i = 0; i < discret.size(); i++){
        cout<<discret[i]<<" ";
    }
    cout<<"]"<<endl;
    
    cout<<"sigma_s0 = [";
    for(int i = 0; i < sigma_s0.size(); i++){
        cout<<sigma_s0[i]<<" ";
    }
    cout<<"]"<<endl;
    
    cout<<"sigma_s1 = [";
    for(int i = 0; i < sigma_s1.size(); i++){
        cout<<sigma_s1[i]<<" ";
    }
    cout<<"]"<<endl;
    
    cout<<"sigma_a = [";
    for(int i = 0; i < sigma_a.size(); i++){
        cout<<sigma_a[i]<<" ";
    }
    cout<<"]"<<endl;
    
    cout<<"Q = [";
    for(int i = 0; i < Q.size(); i++){
        cout<<Q[i]<<" ";
    }
    cout<<"]"<<endl;
    
    cout<<"phi_0_0 = [";
    int temp = 0;
    int temp2 = 0;
    for(int i = 0; i < phi_0_0.size(); i++){
        //Divide the phi_0_0 readout by material sections
        if(temp < discret.size() && temp2 == discret[temp] && i!= phi_0_0.size()-1){
            cout<<endl;
            temp++;
            temp2 = 0;
        }
        cout<<phi_0_0[i]<<" ";
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


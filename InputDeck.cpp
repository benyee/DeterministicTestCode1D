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
    cout<<fileName<<endl;
    
    X.push_back(1);
    
    discret.push_back(5);
    
    sigma_s0.push_back(0);
    
    sigma_s1.push_back(0);
    
    sigma_a.push_back(1);
    
    Q.push_back(1);
    
    for(int j = 0; j < discret.size(); j++){
        for(int i = 0; i < discret[j]; i++){
            phi_0_0.push_back(0);
        }
    }
    //
    
    
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
    for(int i = 0; i < phi_0_0.size(); i++){
        cout<<phi_0_0[i]<<" ";
        //Divide the phi_0_0 readout by material sections
        if(temp < discret.size() && i == discret[temp] && i!= phi_0_0.size()-1){
            cout<<endl;
            temp += 1;
        }
    }
    cout<<"]"<<endl;
}


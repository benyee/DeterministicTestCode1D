//
//  Utilities.cpp
//  
//
//  Created by Ben Yee on 9/21/13.
//
//

#include "Utilities.h"

vector<double> Utilities::calc_mu_n(int N){
    vector<double> out;
    return out;
}
vector<double> Utilities::calc_w_n(int N){
    vector<double> out;
    return out;
}
vector<double> Utilities::lege_coef(int N){
    vector< vector<double> > temp;
    for(int i = 0; i<=N; i++){
        vector<double> temp2;
        temp.push_back(temp2);
    }
    //P_0:
    temp[0].push_back(1);
    if(N==0){return temp[0];}
    
    //P_1:
    temp[1].push_back(0);
    temp[1].push_back(1);
    if(N==1){return temp[1];}
    
    //Fill in zeros:
    temp[0].push_back(0);
    for(int i = 2; i<=N; i++){
        temp[0].push_back(0);
        temp[1].push_back(0);
    }
    
    //All the other ones:
    for(int i = 2; i <=N; i++){
        temp[i].push_back(-(i-1)*temp[i-2][0]/i);
        for(int j =1; j<=i;j++){
            temp[i].push_back(((2*i-1)*temp[i-1][j-1]-(i-1)*temp[i-2][j])/i);
        }
        for(int j=i+1;j<=N;j++){
            temp[i].push_back(0);
        }
    }
    return temp[N];
}
vector<double> Utilities::lege_roots(int N){
    vector<double> out;
    return out;
}

double Utilities::lege_eval(vector<double> coeff,double x){
    double out = 0;
    for(unsigned int i = 0; i<coeff.size();i++){
        out += coeff[i]*pow(x,i);
    }
    return out;
}
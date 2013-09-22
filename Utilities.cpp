//
//  Utilities.cpp
//  
//
//  Created by Ben Yee on 9/21/13.
//
//

#include "Utilities.h"

vector<double> Utilities::calc_mu_n(int N,double tol){
    vector<double> out;
    double x,x1;
    for(int i = 1; i<=N;i++){
        x = cos(PI*(i-0.25)/(N+0.25));
        do{
            x1 = x;
            x -= lege_eval(N,x)/lege_eval_diff(N,x);
        }while(abs(x1-x)>tol);
        out.push_back(x);
    }
    return out;
}
vector<double> Utilities::calc_w_n(vector<double> mu_n){
    vector<double> out;
    double diff;
    for(int i = 0; i< (int)mu_n.size(); i++){
        diff = lege_eval_diff(mu_n.size(),mu_n[i]);
        out.push_back(2/((1-pow(mu_n[i],2))*pow(diff,2)));
    }
    return out;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
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

double Utilities::lege_eval(vector<double> coeff,double x){
    double out = 0;
    for(int i = 0; i<(int)coeff.size();i++){
        out += coeff[i]*pow(x,i);
    }
    return out;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void Utilities::print_uivector(vector<unsigned int> input_vector,char space){
    cout<<"[";
    for(unsigned int i = 0; i<input_vector.size(); i++){
        cout<<input_vector[i]<<space;
    }
    cout<<"]"<<endl;
}
void Utilities::print_ivector(vector<int> input_vector,char space){
    cout<<"[";
    for(unsigned int i = 0; i<input_vector.size(); i++){
        cout<<input_vector[i]<<space;
    }
    cout<<"]"<<endl;
}
void Utilities::print_dvector(vector<double> input_vector,char space){
    cout<<"[";
    for(unsigned int i = 0; i<input_vector.size(); i++){
        cout<<input_vector[i]<<space;
    }
    cout<<"]"<<endl;
}
double Utilities::inf_norm(vector<double> v1, vector<double> v2){
    double inf_norm = 0;
    double temp;
    for(unsigned int i = 0; i<min(v1.size(),v2.size()); i++){
        temp = abs(v2[i]-v1[i]);
        if(temp>inf_norm){inf_norm = temp;}
    }
    return inf_norm;
}
vector<double> Utilities::vector_add(vector<double> v1, vector<double> v2){
    vector<double> sum;
    for(unsigned int i = 0; i<min(v1.size(),v2.size()); i++){
        sum.push_back(v1[i]+v2[i]);
    }
    return sum;
}
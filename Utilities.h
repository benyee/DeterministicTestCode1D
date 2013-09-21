//
//  Utilities.h
//  
//
//  Created by Ben Yee on 9/21/13.
//
//

#ifndef ____Utilities__
#define ____Utilities__

#include <iostream>
#include <vector>

using namespace std;

class Utilities{
public:
    vector<double> calc_mu_n(int N);
    vector<double> calc_w_n(int N);
    vector<double> lege_coef(int N);
    vector<double> lege_roots(int N);

    double lege_eval(int N,double x){return lege_eval(lege_coef(N),x);}
    double lege_eval_diff(int N,double x){return lege_eval_diff(lege_coef(N),x);}
    double lege_eval(vector<double> coeff,double x);
    double lege_eval_diff(vector<double> coeff,double x);
    
private:
    Utilities();
};


#endif /* defined(____Utilities__) */

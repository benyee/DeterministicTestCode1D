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
#include <math.h>

using namespace std;

class Utilities{
public:
    static vector<double> calc_mu_n(int N);
    static vector<double> calc_w_n(int N);
    static vector<double> lege_coef(int N);
    static vector<double> lege_roots(int N);

    static double lege_eval(int N,double x){return lege_eval(lege_coef(N),x);}
    static double lege_eval_diff(int N,double x){return N*(x*lege_eval(N,x)-lege_eval(N-1,x))/(x*x-1);}
    static double lege_eval(vector<double> coeff,double x);
    
private:
    Utilities();
};


#endif /* defined(____Utilities__) */

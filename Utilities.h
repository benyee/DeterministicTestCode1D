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
#include <cmath>

using namespace std;

class Utilities{
public:
    static const double PI = 3.14159265359;
    
    //Functions for calculating Gaussian quadrature constants:
    static vector<double> calc_mu_n(int N,double tol = 1e-10);
    static vector<double> calc_w_n(int N){return calc_w_n(calc_mu_n(N));}
    static vector<double> calc_w_n(vector<double> mu_n);
    
    //Functions for calculating and evaluating Legendre Polynomials:
    static vector<double> lege_coef(int N);
    static double lege_eval(int N,double x){return lege_eval(lege_coef(N),x);}
    static double lege_eval_diff(int N,double x){return N*(x*lege_eval(N,x)-lege_eval(N-1,x))/(x*x-1);}
    static double lege_eval(vector<double> coeff,double x);
    
    //Miscellaneous:
    static void print_uivector(vector<unsigned int> input_vector,char space = ' ');
    static void print_ivector(vector<int> input_vector,char space = ' ');
    static void print_dvector(vector<double> input_vector,char space = ' ');

    static double inf_norm(vector<double> &v1, vector<double> &v2);
    static double p_norm(vector<double> v1, vector<double> v2, unsigned int p);
    static vector<double> vector_add(vector<double> v1, vector<double> v2);
    static vector<double> vector_subtract(vector<double> v1, vector<double> v2);
    
    
    //Remove the last four characters of a string:
    static string remove_ext(string input_string){return input_string.substr(0,input_string.length());}
    
private:
    Utilities();
};


#endif /* defined(____Utilities__) */

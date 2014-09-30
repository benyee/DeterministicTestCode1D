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
    static void print_dmatrix(vector<vector<double> > A,char space = ' ');

    static double inf_norm(vector<double> &v1, vector<double> &v2);
    static double p_norm(vector<double> v1, vector<double> v2, unsigned int p);
    
    //phi_error computes the relative error between two solutions.  norm = p > 0 for p_norm, norm = -1 for inf norm
    //the solution vectors should be 2xN and 2xM vectors respectively where N and M are the number of spatial points.  The first component should be the x-coordinates while the 2nd component should be the scalar flux values at those coordinates.
    static double phi_error(vector<vector<double> > &ref_soln, vector<vector<double> > &soln, int norm);
    
    //Vector math:
    static vector<double> vector_add(vector<double> v1, vector<double> v2);
    static vector<double> vector_subtract(vector<double> v1, vector<double> v2);
    
    static double symmetry_checker(vector<double> &v1);
    static bool nan_checker(vector<double> &v1);
    
    //Find out if a vector is blowing up
    static int blowUpChecker(vector<double> &test, double thres = 10^8);
    static vector<int> blowUpChecker(vector<vector<double> > &test, double thres = 10^8);
    
    //Linear algebra functions:
    static vector<double> solve_tridiag(vector<vector<double> > A, const vector<double> &b);
    
    //Remove the last four characters of a string:
    static string remove_ext(string input_string){return input_string.substr(0,input_string.length());}
    
    //Split up a phi vector into cell + edge fluxes
    static void split_Phi(const vector<double> &phi_all, vector<double> &phi_edge, vector<double> &phi_cent);
    
    //Solve for kappa where kappa is the constant in the solution psi_n(x) = a_n e^{-\Sigma_t \kappa x}
    static double find_kappa(double c, const vector<double> &mu_n, const vector<double> &w_n, double tol = 1e-6);
    
private:
    static double kappa_fun(double c, const vector<double> &mu_n, const vector<double> &w_n, double kappa); //Helper function for find_kappa()
    static double kappa_fun_deriv(double c, const vector<double> &mu_n, const vector<double> &w_n, double kappa); //Helper function for find_kappa(), derivative of kappa_fun()
    
    Utilities();
};


#endif /* defined(____Utilities__) */

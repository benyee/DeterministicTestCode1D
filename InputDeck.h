//
//  InputDeck.h
//  
//
//  Created by Ben Yee on 9/20/13.
//
//

#ifndef ____InputDeck__
#define ____InputDeck__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Utilities.h"


using namespace std;

class InputDeck{
public:
    InputDeck(string inputName="input.txt");
    ~InputDeck();
    
    int loadInputDeck();
    
    vector<double> getX(){return X;}
    vector<unsigned int> getdiscret(){return discret;}
    vector<double> getsigma_s0(){return sigma_s0;}
    vector<double> getsigma_s1(){return sigma_s1;}
    vector<double> getsigma_a(){return sigma_a;}
    vector<double> getQ(){return Q;}
    vector<double> getQ_lin(){return Q_lin;}
    
    void setX(vector<double> input){X = input;}
    void setdiscret(vector<unsigned int> input){discret=input;}
    void setsigma_s0(vector<double> input){sigma_s0=input;}
    void setsigma_s1(vector<double> input){sigma_s1=input;}
    void setsigma_a(vector<double> input){sigma_a=input;}
    void setQ(vector<double> input){Q=input;}
    void setQ_lin(vector<double> input){Q_lin=input;}
    
    void setbc(int input[2]){bc[0] = input[0]; bc[1]=input[1];}
    void setN(int input){N=input;}
    void setaccel_mode(unsigned int input){accel_mode = input;}
    void setalpha_mode(unsigned int input){alpha_mode = input;}
    void settol(double input){tol=input;}
    void sethasLinearTerms(bool input){hasLinearTerms=input;}
    
    int* getbc(){return bc;}
    int getN(){return N;}
    unsigned int getaccel_mode(){return accel_mode;}
    unsigned int getalpha_mode(){return alpha_mode;}
    double gettol(){return tol;}
    bool gethasLinearTerms(){return hasLinearTerms;}
    
    string getfileName(){return fileName;}
    void setfileName(string input){fileName=input;}
    
    vector<double> getphi_0_0(){return phi_0_0;}
    vector<double> getphi_1_0(){return phi_0_0;}
    vector<double> getphi_0_0_lin(){return phi_0_0_lin;}
    vector<double> getphi_1_0_lin(){return phi_0_0_lin;}
    vector<double> getpsi_bl(){return psi_bl;}
    vector<double> getpsi_br(){return psi_br;}
    
    void setphi_0_0(vector<double> input){phi_0_0=input;}
    void setphi_1_0(vector<double> input){phi_0_0=input;}
    void setphi_0_0_lin(vector<double> input){phi_0_0_lin=input;}
    void setphi_1_0_lin(vector<double> input){phi_0_0_lin=input;}
    void setpsi_bl(vector<double> input){psi_bl=input;}
    void setpsi_br(vector<double> input){psi_br=input;}
    
    void readValues(); //Read out input values
    
    void setDefaultValues();

private:
    string fileName; //input file name
    
    //The following vectors should always have the same length:
    vector<double> X; //location of material boundaries in cm (left boundary is always 0, DO NOT INCLUDE 0 AS THE FIRST ELEMENT OF X, IT IS IMPLIED)
    vector<unsigned int> discret;  //number of elements to divide each region into
    vector<double> sigma_s0; //isotropic scattering cross sections in cm^{-1}
    vector<double> sigma_s1; //anisotropic scattering cross sections in cm^{-1}
    vector<double> sigma_a; //absorption cross section in cm^{-1}
    vector<double> Q;  //isotropic source (in neutrons per cm^3)
    vector<double> Q_lin;
    //Note: Q and Q_lin are defined such that, for material section j with x
    // between x_{j-1} and x_{j}, the source in that region is given by
    // Q(x) = Q + Q_lin*(x - x_{j-1})
    
    vector<double> phi_0_0; //initial guess for source iteration, should be of size (\sum_{i=0}^{size(xborders)} discret[i] )
    vector<double> phi_1_0; //initial guess for source iteration, should be of size (\sum_{i=0}^{size(xborders)} discret[i] )
    vector<double> phi_0_0_lin; //initial guess for source iteration, should be of size (\sum_{i=0}^{size(xborders)} discret[i] )
    vector<double> phi_1_0_lin; //initial guess for source iteration, should be of size (\sum_{i=0}^{size(xborders)} discret[i] )
    
    vector<double> psi_bl; //user specified left boundary, should be of size N/2
    vector<double> psi_br; //user specified right boundary, should be of size N/2
    
    int bc[2]; // boundary conditions for left and right end of domain.  0 = vaccuum, 1 = reflective, 2 = user specified
    
    unsigned int N;  // number of angular ordinates
    unsigned int accel_mode; //type of acceleration used.
    /* 0 = no acceleration
     1 = CMFD
     2 = pCMFD
     3 = mpCMFD
    */
    unsigned int alpha_mode; //type of finite difference scheme.
    /*0 = diamond difference
     1 = step up
     2 = step characteristic
     3 = characteristic alternative (3.1 of notes)
     10 = linear characteristic (3.2 of notes)
     11 = alternative linear characteristic (3.3 of notes)
     20 = linear discontinuous (4.1 of notes)
     */
    double tol; //tolerance for convergence
    bool hasLinearTerms;
    
    bool searchForInput(ifstream &file, string inp);  //Search for the input named inp in file, returns 0 if not found, returns 1 if found
};

#endif /* defined(____InputDeck__) */

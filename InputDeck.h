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
    vector<int> getdiscret(){return discret;}
    vector<double> getsigma_s0(){return sigma_s0;}
    vector<double> getsigma_s1(){return sigma_s1;}
    vector<double> getsigma_a(){return sigma_a;}
    vector<double> getQ(){return Q;}
    
    int* getbc(){return bc;}
    int getN(){return N;}
    int getalpha_mode(){return alpha_mode;}
    double gettol(){return tol;}
    
    string getfileName(){return fileName;}
    
    vector<double> getphi_0_0(){return phi_0_0;}
    vector<double> getphi_1_0(){return phi_0_0;}
    
    void readValues(); //Read out input values

private:
    string fileName; //input file name
    
    //The following vectors should always have the same length:
    vector<double> X; //location of material boundaries in cm (left boundary is always 0)
    vector<int> discret;  //number of elements to divide each region into
    vector<double> sigma_s0; //isotropic scattering cross sections in cm^{-1}
    vector<double> sigma_s1; //anisotropic scattering cross sections in cm^{-1}
    vector<double> sigma_a; //absorption cross section in cm^{-1}
    vector<double> Q;  //isotropic source (in neutrons per cm^3)
    
    vector<double> phi_0_0; //initial guess for source iteration, should be of size (\sum_{i=0}^{size(xborders)} discret[i] )
    vector<double> phi_1_0; //initial guess for source iteration, should be of size (\sum_{i=0}^{size(xborders)} discret[i] )
    
    int bc[2]; // boundary conditions for left and right end of domain.  0 = vaccuum, 1 = reflective, 2 = user specified
    
    int N;  // number of angular ordinates
    int alpha_mode; //type of finite difference scheme.  0 = diamond difference, 1 = step up, 2 = step characteristic
    double tol; //tolerance for convergence
    
    bool searchForInput(ifstream &file, string inp);  //Search for the input named inp in file, returns 0 if not found, returns 1 if found
};

#endif /* defined(____InputDeck__) */

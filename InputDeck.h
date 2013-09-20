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
    
    vector<double> getphi_0_0(){return phi_0_0;}
    
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
    
    int bc[2]; // boundary conditions for left and right end of domain
    
    bool searchForInput(ifstream &file, string inp);  //Search for the input named inp in file, returns 0 if not found, returns 1 if found
};

#endif /* defined(____InputDeck__) */

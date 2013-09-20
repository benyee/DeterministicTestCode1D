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
#include <string>
#include <vector>

using namespace std;

class InputDeck{
public:
    InputDeck(string inputName="input.txt");
    ~InputDeck();
    
    int loadInputDeck();
    
    vector<int> getX(){return X;}
    vector<int> getdiscret(){return discret;}
    vector<int> getsigma_s0(){return sigma_s0;}
    vector<int> getsigma_s1(){return sigma_s1;}
    vector<int> getsigma_a(){return sigma_a;}
    vector<int> getQ(){return Q;}
    
    vector<int> getphi_0_0(){return phi_0_0;}
    
    void readValues(); //Read out input values

private:
    string fileName; //input file name
    
    //The following vectors should always have the same length:
    vector<int> X; //location of material boundaries in cm (left boundary is always 0)
    vector<int> discret;  //number of elements to divide each region into
    vector<int> sigma_s0; //isotropic scattering cross sections in cm^{-1}
    vector<int> sigma_s1; //anisotropic scattering cross sections in cm^{-1}
    vector<int> sigma_a; //absorption cross section in cm^{-1}
    vector<int> Q;  //isotropic source (in neutrons per cm^3)
    
    vector<int> phi_0_0; //initial guess for source iteration, should be of size (\sum_{i=0}^{size(xborders)} discret[i] )
    
    int bc[2]; // boundary conditions for left and right end of domain
};

#endif /* defined(____InputDeck__) */

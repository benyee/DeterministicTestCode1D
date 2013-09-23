//
//  SourceIteration.h
//  
//
//  Created by Ben Yee on 9/21/13.
//
//

#ifndef ____SourceIteration__
#define ____SourceIteration__

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

#include "Utilities.h"
#include "InputDeck.h"

using namespace std;

class SourceIteration{
public:
    SourceIteration(InputDeck *input);
    ~SourceIteration();
    int iterate();
    void printOutput(string outfilename="output.txt");
    
private:
    InputDeck *data;
    
    //Grid:
    vector<double> x;
    vector<double> h;
    vector<unsigned int> discret;
    vector<double> x_e;
    
    //Neutron info:
    vector<double> phi_0;
    vector<double> phi_1;
    vector<vector<double> > psi_e;
    vector<vector<double> > psi_c;
    vector< vector<double> > source;
    
    vector<double> sigma_s0; //isotropic scattering cross sections in cm^{-1}
    vector<double> sigma_s1; //anisotropic scattering cross sections in cm^{-1}
    vector<double> sigma_a; //absorption cross section in cm^{-1}
    
    vector< vector<double> > alpha;
    
    vector<double> mu_n;
    vector<double> w_n;
    
    unsigned int J,N;
    int* bc;
    
    void rightIteration();
    void leftIteration();
    void finiteDifference();
    void initializeAlpha();
    void initializeGrid();
    double updatePhi_calcSource();
};

#endif /* defined(____SourceIteration__) */

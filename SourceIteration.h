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
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

#include "Utilities.h"
#include "InputDeck.h"

using namespace std;

class SourceIteration{
public:
    static const unsigned int MAX_IT = 20;
    
    SourceIteration(InputDeck *input,string outputfilename="output.txt");
    ~SourceIteration();
    int iterate();
    void printOutput( bool isPrintingToWindow = false,unsigned int tabwidth=20);
    
private:
    InputDeck *data; //Input deck
    string outfilename; //Name of outpile file
    
    //Fine Grid:
    vector<double> x; //List of x-values at cell centers
    vector<double> h; //List of widths (dx's)
    vector<unsigned int> discret; //Number of spatial cells for each region
    vector<double> x_e; //List of x-values at cell edges
    //Coarse Grid:
    vector<double> x_CM; //List of x-values at cell centers
    vector<double> h_CM; //List of widths (dx's)
    vector<unsigned int> discret_CM; //Number of spatial cells for each region
    vector<double> x_CM_e; //List of x-values at cell edges
    
    //Neutron info:
    vector<double> phi_0; //Cell-averaged scalar flux
    vector<double> phi_1; //Cell-averaged net current in x direction
    vector<vector<double> > psi_e; //Edge angular flux
    vector<vector<double> > psi_c; //Cell-averaged angular flux
    vector< vector<double> > source;  //Source term for source iteration
    
    //Linear neutron info for higher order approximations:
    bool hasLinearTerms;
    vector<double> phi_0_lin; //Cell-averaged scalar flux
    vector<double> phi_1_lin; //Cell-averaged net current in x direction
    vector<vector<double> > psi_c_lin; //Cell-averaged angular flux
    vector< vector<double> > source_lin;  //Source term for source iteration
    
    //Cross sections:
    vector<double> D_actual_CM; //"averaged" diffusion coefficient for coarse grid
    vector<double> opt_CM_a; //"averaged" absorption optical lengths for coarse grid
    
    vector<double> sigma_s0; //isotropic scattering cross sections in cm^{-1}
    vector<double> sigma_s1; //anisotropic scattering cross sections in cm^{-1}
    vector<double> sigma_t; //absorption cross section in cm^{-1}
    double c;  //maximum value of sigma_s0/sigma_t
    
    vector< vector<double> > alpha; //Finite difference coefficients
    
    vector<double> mu_n; //mu values for S_N approximation
    //Convention: mu_n[0] > mu_n[1] > ... mu_n[N-1]
    vector<double> w_n; //weights for S_N approximation
    
    unsigned int J,N; //number of spatial cells, order of S_N approximation
    int* bc; //boundary conditions
    unsigned int alpha_mode;
    unsigned int accel_mode;
    
    int it_num; //Iteraiton number
    double old_error; //stores ||phi_{i}-phi_{i-1}|| from the previous iteration
    
    void rightIteration(); //Sweep left to right
    void leftIteration(); //Sweep right to left
    void finiteDifference(); //Calculate cell-averaged angular fluxes
    void cmfd(); //Perform cmfd acceleration
    void pcmfd(); //Perform pcmfd acceleration
    void initializeAlpha(); //Calculate alpha values
    void initializeGrid(); //Calculate values associated with grid locations
    double updatePhi_calcSource(bool usePsi = true); //Update fluxes and currents, calculate new source, calculate difference between new and old scalar flux
    vector<double> calcEdgePhi(int num); //Integrate the edge fluxes to get scalar edge fluxes and edge currents.  This is not necessary for the source iteration procedure and can be performed at the end.  num = 0 for scalar edge flux, num = 1 for edge current
    unsigned int checkNegativeFlux(); //Returns the number of negative values in phi_0
};

#endif /* defined(____SourceIteration__) */

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
    static const unsigned int MAX_IT = 100000;
    static const unsigned int MAX_IT_accel = 5000;
    static const unsigned int EDGE_ACCEL_MODE = 0;
    static const bool intermedSoln = 0;
    
    SourceIteration(InputDeck *input,string outputfilename="output.txt");
    ~SourceIteration();
    int iterate(bool isPrintingToWindow = true, bool isPrintingToFile = true, bool falseConvCorrection = true);
    void printOutput( bool isPrintingToWindow = false,unsigned int tabwidth=20, bool newFile = false);
    
    bool isConverged;
    
    int get_it_num(){if(isConverged){return it_num;} return -it_num;}
    
    double get_error(){return old_error;}
    double get_spec_rad(){return spec_rad;}
    
    vector<vector<double> > get_psi_e(){return psi_e;}
    vector<vector<double> > get_solution(); //Returns the phi and x values from the solution.
    
    void set_diverge(double div){diverge = div;}
    
    void set_outputfilename(string filename){outfilename = filename;}
    
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
    vector<double> edgePhi0;
    vector<double> edgePhi1;
    vector<double> phi_0; //Cell-averaged scalar flux
    vector<double> old_phi_0;
    vector<double> phi_1; //Cell-averaged net current in x direction
    vector<vector<double> > psi_e; //Edge angular flux
    vector<vector<double> > psi_c; //Cell-averaged angular flux
    vector< vector<double> > source;  //Source term for source iteration
    vector< vector<double> > source_edge;  //Source term for source iteration
    
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
    double c;  //This is max(sigma_s0/sigma_t) unless this gives c >= 1.
    //In that case, c = max( sigma_s0/(sigma_t+DB^2) ))
    
    vector< vector<double> > alpha; //Finite difference coefficients
    
    vector<double> mu_n; //mu values for S_N approximation
    //Convention: mu_n[0] > mu_n[1] > ... mu_n[N-1]
    vector<double> w_n; //weights for S_N approximation
    
    unsigned int J,N; //number of spatial cells, order of S_N approximation
    int* bc; //boundary conditions
    unsigned int alpha_mode;
    unsigned int accel_mode;
    
    int it_num; //Iteraiton number
    double init_error; //Stores the first error, used to check for divergence.
    double diverge; //Divergence criterion.  Algorithm diverges if error/init_error > diverge.
    double old_error; //stores ||phi_{i}-phi_{i-1}|| from the previous iteration
    double spec_rad; //stores the latest estimate for the spectral radius (average of all the previous spectral radii)
    
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

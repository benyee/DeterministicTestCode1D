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
#include <map>

#include "Utilities.h"
#include "InputDeck.h"

using namespace std;

typedef std::map<string, double> Dict;

class SourceIteration{
public:
    static const unsigned int MAX_IT = 10000;
    static const unsigned int MAX_IT_accel = 5000;
    static const unsigned int EDGE_ACCEL_MODE = 1;
    static const unsigned int ERROR_HISTORY_LENGTH = 100;
    static const bool intermedSoln = 0;
    
    bool isConverged;
    
    SourceIteration(InputDeck *input,string outputfilename="output.txt");
    ~SourceIteration();
    int iterate(bool isPrintingToWindow = true, bool isPrintingToFile = true, bool falseConvCorrection = true);
    vector<vector<double> > get_solution(bool oldSolution = 0); //Returns the phi and x values from the solution.
    void printOutput( bool isPrintingToWindow = false,unsigned int tabwidth=20, bool newFile = false);
    void printDictionary();
    
    unsigned int checkNegativeFlux(bool checkEdgeFluxes = false , bool printFluxes = false); //Returns the number of negative values in phi_0
    unsigned int checkNegativeAngularFlux(bool printFluxes = false); //Returns the number of negative values in psi_c and psi_e
    
    double get_error(){return old_error;}
    double get_exitcurr(bool isExitOnRight=true);
    double get_it_num(){if(isConverged){return it_num;} return -it_num;}
    double get_kappa(){return kappa;}
    double get_Qhatnorm(){if(alpha_mode >= 40 && alpha_mode != 41 && alpha_mode < 40){return 0;} return Utilities::p_norm(Qhat_edge,2);}
    vector<double> get_Qhat_edge(){return Qhat_edge;}
    vector<double> get_rho(){return rho;}
    vector<vector<double> > get_psi_e(){return psi_e;}
    double get_spec_rad(){return spec_rad;}
    void set_diverge(double div){diverge = div;}
    void set_outputfilename(string filename){outfilename = filename;}
    void set_rho(vector<double> rho_input){rho = rho_input;}
    void set_kappa(double kappa_input){kappa = kappa_input;}
    
    
private:
    InputDeck *data; //Input deck
    string outfilename; //Name of outpile file
    
    Dict variable_status; //Keeps track of which iteration number each variable is at
    
    //Fine Grid:
    vector<double> x; //List of x-values at cell centers
    vector<double> h; //List of widths (dx's)
    vector<double> h_avg; //List of h_{j+1/2} values.  h_{j+1/2} is the average of h_j and h_{j+1}
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
    vector<double> edgePhi0_MB2_R;
    vector<double> edgePhi0_MB2_L;
    vector<double> phi_0; //Cell-averaged scalar flux
    vector<double> old_phi_0;
    vector<double> old_edgePhi0;
    vector<double> phi_1; //Cell-averaged net current in x direction
    vector<vector<double> > psi_e; //Edge angular flux
    vector<vector<double> > psi_c; //Cell-averaged angular flux
    vector< vector<double> > source;  //Source term for source iteration
    vector< vector<double> > source_edge;  //Source term for source iteration
    vector<double> Qhat_edge; //Correction to multiple balance method that helps enforce Fick's Law
    //mu_n*Qhat_edge is subtracted from the source_edge terms in the multiple balance auxiliary equations
    vector<double> phi2_plus;
    vector<double> phi2_minus;
    vector<double> phi2e_plus;
    vector<double> phi2e_minus;
    
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
    vector<double> Sah_avg; //Sah_avg_{j-1/2} = avg of sigma_t in cells j and j+1 weighted by h_j, times h_avg[j].
    vector<double> Sth_avg; //Sth_avg_{j-1/2} = avg of sigma_t in cells j and j+1 weighted by h_j, times h_avg[j].
    double c;  //This is max(sigma_s0/sigma_t) unless this gives c >= 1.
    //In that case, c = max( sigma_s0/(sigma_t+DB^2) ))
    
    vector<double> Q_avg; //Q_avg_{j-1/2} = Q_avg[j] = avg of sigma_t in cells j and j+1
    
    double kappa; //Constant in psi_n(x) = a_n e^{-\Sigma_t \kappa x}
    vector<double> rho; //dx becomes rho*dx in the multiple balance aux. equations
    
    vector< vector<double> > alpha; //Finite difference coefficients
    
    vector<double> mu_n; //mu values for S_N approximation
    //Convention: mu_n[0] > mu_n[1] > ... mu_n[N-1]
    vector<double> w_n; //weights for S_N approximation
    vector< vector<double> > w_n_MB2; //special weights for MB-2 mehtod
    
    unsigned int J,N; //number of spatial cells, order of S_N approximation
    int* bc; //boundary conditions
    unsigned int alpha_mode;
    unsigned int accel_mode;
    
    double it_num; //Iteration number.  Can have half iterations.
    double init_error; //Stores the first error, used to check for divergence.
    double diverge; //Divergence criterion.  Algorithm diverges if error/init_error > diverge.
    double old_error; //stores ||phi_{i}-phi_{i-1}|| from the previous iteration
    double spec_rad; //stores the latest estimate for the spectral radius (average of all the previous spectral radii)
    vector<double> error_history;
    
    void accelerate_edgePhi0(vector<double> preaccel_phi_0); //accelerate edge fluxes
    void accelerate_edgePhi0_MB2(vector<double> preaccel_phi_0); //accelerate edge fluxes for the MB2 case
    void accelerate_MB3(); //Diffusion based acceleration scheme for MB3
    vector<double> calcEdgePhi(int num); //Integrate the edge fluxes to get scalar edge fluxes and edge currents.  This is not necessary for the source iteration procedure and can be performed at the end.  num = 0 for scalar edge flux, num = 1 for edge current
    void cmfd(); //Perform cmfd acceleration
        void pcmfd(); //Perform pcmfd acceleration
    void finiteDifference(); //Calculate cell-averaged angular fluxes
    void initializeAlpha(); //Calculate alpha values
    void initializeDictionary(); //Initialize dictionary.
    void initializerho(); //Initialize rho for MB-3
    void initializeGrid(); //Calculate values associated with grid locations
    void initializew_n_MB2(); //Initialize w_n_MB2
    void rightIteration(); //Sweep left to right
        void leftIteration(); //Sweep right to left
    void setEdgePhi_toAvgPhi(); //Used to make a crude improvement of the initial edge flux guess
    void updateEdgePhi0_MB2(); //Updates the special edge fluxes in the MB-3 case
    double updatePhi_calcSource(bool usePsi = true); //Update fluxes and currents, calculate new source, calculate difference between new and old scalar flux
    void updatePhi2(); //
    void updateQhat_edge(); //For MB-3 (modified MB method)
};

#endif /* defined(____SourceIteration__) */

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
#include <math>
#include <vector>

#include "Utilities.h"
#include "InputDeck.h"

class SourceIteration{
public:
    SourceIteration(InputDeck input);
    int iterate();
    void printOutput();
    
private:
    InputDeck data;
    vector<double> phi_0;
    vector<double> phi_1;
    vector<vector<double> > psi_e;
    vector<vector<double> > psi_e;
};

#endif /* defined(____SourceIteration__) */

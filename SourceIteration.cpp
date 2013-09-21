//
//  SourceIteration.cpp
//  
//
//  Created by Ben Yee on 9/21/13.
//
//

#include "SourceIteration.h"

SourceIteration::SourceIteration(InputDeck input){
    data = input;
}

void SourceIteration::iterate(){
    int bc[2] = data.getbc();
    return;
}

void SourceIteration::printOutput(){
    cout<<"Printing to a file named"<<Utilities::remove_ext(input.getfileName())<<endl;
    cout<<"Output goes here"<<endl;
}
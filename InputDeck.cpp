//
//  InputDeck.cpp
//  
//
//  Created by Ben Yee on 9/20/13.
//
//

#include "InputDeck.h"

InputDeck::InputDeck(string inputName){
    fileName = inputName;
    loadInputDeck();
}

int InputDeck::loadInputDeck(){
    cout<<fileName<<endl;
    return 0;
}

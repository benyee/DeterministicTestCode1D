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

using namespace std;

class InputDeck{
public:
    InputDeck(string inputName);
    ~InputDeck();

private:
    string fileName;
    int loadInputDeck();
};

#endif /* defined(____InputDeck__) */

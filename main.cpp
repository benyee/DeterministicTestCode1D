//
//  main.cpp
//  
//
//  Created by Ben Yee on 9/20/13.
//
//

#include <iostream>
#include <string>

#include "InputDeck.h"

using namespace std;

int main ()
{
    string ifilename = "input.txt";
    InputDeck *input = new InputDeck(ifilename);
    cout << "Hello World!"<<endl;
    return 0;
}

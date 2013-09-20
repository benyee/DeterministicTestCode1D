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
    cout << "Hello world!"<<endl;
    InputDeck *input = new InputDeck();
    int debug =input->loadInputDeck();
    input->readValues();
    if(debug){
        cout << "There is an issue with the sizes of the input vectors!" << endl;
        return 1;
    }
    cout << "Vector sizes look good!"<<endl;
    cout << "Goodbye world!"<<endl;
    return 0;
}

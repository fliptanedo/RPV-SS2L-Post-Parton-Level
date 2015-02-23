// FlipCommandFileFixer.h
// For modifying Pythia Command Files
// INCLUDE GUARD
// #ifndef __FLIPLHE_H_INCLUDED__
// #define __FLIPLHE_H_INCLUDED__

#include <string>              
#include <sstream>              // for string stream
#include <iostream>             // for i don't know
#include <iomanip>              // for setting precision?
#include <fstream>              // for file in/out
//
#include <algorithm>            // These four are all from
#include <functional>           //  http://stackoverflow.com/
#include <cctype>               //  questions/216823/whats-the-
#include <locale>               //  best-way-to-trim-stdstring
using namespace std;            

int getnevents(std::string &lhefile);
//
// Usage: feed in lhe filename, returns nevents value



// END INCLUDE GUARD
// #endif // __FLIPLHE_H_INCLUDED__


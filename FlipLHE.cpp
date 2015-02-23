/******************************************************************************** 
*   FlipCommandFileFixer.cpp by Flip Tanedo (pt267@cornell.edu)                 *
*   Code for PartonRPV.cc, Aug 2012                                             *
*   Contains functions used for modifying spc and cmnd files                    *
********************************************************************************/

#include "FlipLHE.h"

int getnevents(std::string &lhefile){
    
    // FILE STREAMS
    ifstream instream;
    instream.open(lhefile.c_str());
    
    int result = 0;
    
    string line;
    bool keepgoing(true);
    bool foundit(false);
    
    while (keepgoing){
        getline(instream, line);
        if (line.find("nevents") != string::npos){
           keepgoing = false;
           foundit = true;
        }
        if (instream.eof())
            keepgoing = false;
    }
    
    if (foundit){
        size_t startpos = line.find_first_not_of(" \t");
    
        string test2 = line.substr(startpos, line.length());
        string numberwang = test2.substr(0,test2.find(' '));
    
        result = atoi(numberwang.c_str());
    }

    return result;
}




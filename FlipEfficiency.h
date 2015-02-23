// FlipEfficiency.h
// INCLUDE GUARD
// #ifndef __FLIPEFFICIENCY_H_INCLUDED__
// #define __FLIPEFFICIENCY_H_INCLUDED__


#include "Pythia.h"                         // Include Pythia headers
#include <fastjet/ClusterSequence.hh>       // fastjet clustering
#include <cmath>                            // for error function
#include <sstream>                          // for string stream
#include <time.h>                           // for random number seed
#include <iostream>                         // for i don't know
#include <iomanip>                          // for setting precision?
#include <fstream>                          // for file in/out
using namespace std;

struct signalregion{
    // this is just a data structure to hold the signal region cuts
    // see table 2 of SUS-12-017-pas
    unsigned int minJets;
    unsigned int minbJets;
    double minMET;
    double minHT;
    bool plusplus;          // allow same sign + charge leptons
    bool minusminus;        // allow same sign - charge leptons
};

double signal_efficiency(string, vector< pair<string, int> >&, int);
    // This is our main workhorse, it's defined in a separate file
    // FlipEfficiencySignal.cpp
    // Inputs: command file, intermediate count vector, signal region index

double BG_efficiency(Pythia8::Pythia&, vector< pair<string, int> >&, int, int);
    // This is just a straight copy of signal_efficiency, but now it takes
    // in an lhe file in addition to the command file
    // Defined in FlipEfficiencySignal.cpp
    // Inputs: pythia object, count vector, signal region index, # event
    // ** eventually I shuold merge these two functions

double signal_efficiency_b(string, vector< pair<string, int> >&, int);
    // Same as signal_efficiency, but with b-tagging!
    // Oct 15 2013
    // This is our main workhorse, it's defined in a separate file
    // FlipEfficiencySignal.cpp
    // Inputs: command file, intermediate count vector, signal region index


    
/******************************************************************************** 
*   Helper functions that calculate intermediate steps, output, etc.            *
********************************************************************************/

void read_count(vector< pair<string, int> >);
void fill_vector(vector< pair<string, int> > &, string, int);
double get_deltaR(fastjet::PseudoJet, fastjet::PseudoJet);

bool lepton_kinematic_cut(pair<int, fastjet::PseudoJet>);
bool jet_kinematic_cut(pair<int, fastjet::PseudoJet>);
bool lepton_selection_cut(pair<int, fastjet::PseudoJet>);
bool lepton_ID_eff(pair<int, fastjet::PseudoJet>);
bool lepton_iso_eff(    pair<int, fastjet::PseudoJet>, 
                        vector< pair<int, fastjet::PseudoJet> >);
bool b_selection_efficiency(pair<int, fastjet::PseudoJet>);
bool lepton_trig_efficiency(vector< pair<int, fastjet::PseudoJet> >);
bool METefficiency(double, double);
bool HTefficiency(double, double);

void fill_signalregions(vector<signalregion>&);




// END INCLUDE GUARD
// #endif // __FLIPEFFICIENCY_H_INCLUDED__


/******************************************************************************** 
*   RPVscanFlip.cc by Flip Tanedo (pt267@cornell.edu)                           *
*   Scan of RPV stop/gluino parameter space to determine the SS2L reach         *
*   31 Aug 2012                                                                 *
*   - no longer depends on Peter Skands' hacked Pythia 8.165 code               *
*   - uses FlipCommandFileFixer.h for 
*   - uses FlipEfficiency.h for doing scan
********************************************************************************/

// Inputs: command file, stop mass, gluino mass, signal region, output filename
//  For example:
//  ./RPVHscanFlip RPVgluinoScan.cmnd 400.0 800.0 8



#include "FlipEfficiency.h"         // all of my functions
#include "FlipCommandFileFixer.h"   // all of my functions
#include <sstream>                  // for string stream
#include <fstream>                  // for file in/out


using namespace std;


int main(int argc, char *argv[]) { 

    // INITIALIZE
    // ----------
    srand((unsigned)time(0));               // Initialize random numbers
    string outfile = "output.dat";          // Output filename
    vector< pair<string, int> > counts;     // counts @ each cut w/ descriptions


    // A bunch of definitions for setting the stop and gluion masses
    // -------------------------------------------------------------
    string mgluino      = "800";                // default Gluino mass        
    string mstop        = "300";                // default stop mass
    string cmndtemp     = "CmndTemp.cmnd";      // default cmnd file template
    string cmndrun      = "CommandRun.cmnd";    // output cmnd file for run
    string spctemp      = "template.spc";       // default spc template
    string spcint       = "intermediate.spc";   // intermediate spc file
    string spcRun       = "spcRun.spc";         // output spc file for run
    //
    string cmndspc      = "SLHA:file = ";       // line to change in cmnd file
    string cmndspcnew   = cmndspc + spcRun;     // ... replace with this
    //
    string blockmass    = "BLOCK MASS";         
    string blockdiv     = "BLOCK";
    string gluinoID     = "1000021";
    string stopID       = "1000006";
    
    
    // Other definitions for the run
    // -----------------------------
    int iSR = 8;        // Signal region #, defined in SUS-12-017
    
    
    // Take in external values
    // -----------------------
    if (argc > 1)  mstop    = argv[1];       // stop mass
    if (argc > 2)  mgluino  = argv[2];       // gluino mass
    if (argc > 3)  iSR      = atoi(argv[3]); // signal region
    if (argc > 4)  cmndtemp = argv[4];       // template command file
    if (argc > 5)  outfile  = argv[5];       // output filename
    if (argc > 6)  spctemp  = argv[6];       // template spectrum file


    // Lines for updating the spectrum
    // --------------------------------
    string gluinoNew    = "   1000021   " + mgluino;    // line replacement
    string stopNew      = "   1000006   " + mstop;      // line replacement



    // UPDATE SPECTRUM
    // ---------------

    // Make sure command file is using the same spc file that we're creating
    if(!FixCommand(cmndtemp, cmndrun, cmndspc, cmndspcnew)) 
        cout << endl << " ERROR in FixCommand, setting " << cmndspcnew << endl;

    // Using template spectrum, set gluino mass. Save to intermediate spectrum
    if(!FixSpectrum(spctemp, spcint, blockmass, blockdiv, gluinoID, gluinoNew))
        cout << endl << "ERROR: FixSpectrum, setting mass " << gluinoNew << endl;
    
    // Using intermediate spectrum, set stop mass. Save to final spectrum
    if(!FixSpectrum(spcint, spcRun, blockmass, blockdiv, stopID, stopNew))
        cout << endl << "ERROR: FixSpectrum, setting mass " << stopNew << endl;
        
        
    // OUTPUT FILE STREAM
    // ------------------

    ofstream outstream;
    outstream.open(outfile.c_str(), ios::app); // append to end of file
    outstream.precision(6); 
    outstream.setf(ios::fixed);
    outstream.setf(ios::showpoint);


    /****************************************************************************
    *   THIS PART DOES THE CALCULATION                                          *
    *****************************************************************************/

    outstream << mstop << "\t" << mgluino << "\t" << iSR << "\t" 
        << (signal_efficiency_b(cmndrun, counts, iSR) * 0.10608) << endl;
        // << (signal_efficiency(cmndrun, counts, iSR) * 0.10608) << endl;
        // 0.10608 = 0.3257^2 from W decays forced to go to leptons (for stats)
        // Why? Because we assume you forced the W to decay leptonically in
        // the command file.

    // // IF YOU WANT VERBOSE SCREEN OUTPUT:
    // cout << "STOP: " << mstop << endl;
    // cout << "GLUINO: " << mgluino << endl;
    // cout << "Signal Region " << iSR << endl; 
    // read_count(counts); // gives intermediate steps
    // cout << endl << endl;   



    /****************************************************************************
    *   CLEAN UP: close filestreams, etc.                                       *
    *****************************************************************************/
    
    outstream.close();
        
    
    return 0;
        
}
    

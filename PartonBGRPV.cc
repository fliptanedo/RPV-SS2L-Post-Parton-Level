/******************************************************************************** 
*   PartonBGPRV.cc by Flip Tanedo (pt267@cornell.edu)                           *
*   Scan of RPV stop/gluino parameter space to determine the SS2L reach         *
*   Focuses on checking background generated from an LHE file                   *
*   7 Oct 2012                                                                  *
*   - debugging                                                                 *
********************************************************************************/

#include "FlipEfficiency.h"         // all of my functions
#include "FlipCommandFileFixer.h"   // all of my functions
#include "FlipLHE.h"
#include <sstream>                  // for string stream
#include <fstream>                  // for file in/out


using namespace std;


int main(int argc, char *argv[]) { 

    // INITIALIZE
    // ----------
    srand((unsigned)time(0));                   // Initialize random numbers

    string command_file = "background.cmnd";    // cmnd file for run
    string input_lhe    = "eventsplus.lhe";         // input LHE file
    string outfile      = "output.dat";         // Output filename
    int iSR             = 8;                    // Signal region #
                                                //  defined in SUS-12-017
    vector< pair<string, int> > counts;         // counts @ each cut 
                                                //  with descriptions
    int nEvents         = getnevents(input_lhe);// read number of events
    
    // TAKE IN EXTERNAL VALUES
    // -----------------------
    if (argc > 1)  input_lhe    = argv[1];       // input LHE file
    if (argc > 2)  command_file = argv[2];       // command file
    if (argc > 3)  iSR          = atoi(argv[3]); // signal region
    if (argc > 4)  outfile      = argv[4];       // output file
    
    
    // SET UP THE PYTHIA OBJECT
    // ------------------------
    Pythia8::Pythia pythia;                 // Declare Pythia object
    // Pythia8::Event& event = pythia.event;   // Declare event as a shortcut
    pythia.readFile(command_file);          // Read in command file

    // int nEvent = pythia.mode("Main:numberOfEvents");
    // int nAbort = pythia.mode("Main:timesAllowErrors");

    pythia.init(input_lhe);                 // Initialize in LHE file

    cout << endl << BG_efficiency(pythia, counts, iSR, nEvents) << endl << endl;
    
    
        
    // // OUTPUT FILE STREAM
    // // ------------------
    // 
    // ofstream outstream;
    // outstream.open(outfile.c_str(), ios::app); // append to end of file
    // outstream.precision(6); 
    // outstream.setf(ios::fixed);
    // outstream.setf(ios::showpoint);


    /****************************************************************************
    *   THIS PART DOES THE CALCULATION                                          *
    *****************************************************************************/

        cout << input_lhe << "\t"
    // outstream << input_lhe << "\t" 
        // << BG_efficiency(cmndrun, input_lhe, counts, iSR) 
        << "\t (times W decay prefactor)" 
        << endl;
        // << (BG_efficiency(cmndrun, input_lhe counts, iSR) * 0.10608) << endl;
        // 0.10608 = 0.3257^2 from W decays forced to go to leptons (for stats)
        // Why? Because we assume you forced the W to decay leptonically in
        // the command file.

    // // IF YOU WANT VERBOSE SCREEN OUTPUT:
    // cout << "STOP: " << mstop << endl;
    // cout << "GLUINO: " << mgluino << endl;
    // cout << "Signal Region " << iSR << endl; 
    read_count(counts); // gives intermediate steps
    // cout << endl << endl;   



    /****************************************************************************
    *   CLEAN UP: close filestreams, etc.                                       *
    *****************************************************************************/
    
    // outstream.close();
        
    
    return 0;
        
}
    

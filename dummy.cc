#include "Pythia.h"                         // Include Pythia headers
#include <string>              
#include "FlipEfficiency.h"         // all of my functions
#include "FlipLHE.h"
using namespace std;


int main(){ 

    // Pythia8::Pythia pythia;                 // Declare Pythia object
    // Pythia8::Event& event = pythia.event;   // Declare event as a shortcut
    // pythia.readFile("background.cmnd");          // Read in command file
    // pythia.init("eventsminiplus.lhe");
    // pythia.next();
    
    string blah = "testtext.txt";

    // cout << "  100000 = nevents ! Number of unweighted events requested " << endl;
    // 
    // string test = "  100000 = nevents ! Number of unweighted events requested ";
    // 
    // size_t startpos = test.find_first_not_of(" \t");
    // 
    // string test2 = test.substr(startpos, test.length());
    // string numberwang = test2.substr(0,test2.find(' '));
    // 
    // cout << endl << numberwang << endl << endl;
    // int numberthang = atoi(numberwang.c_str());
    cout << endl << endl << getnevents(blah) << endl;
    
 
    return 0;
        
}
    
// int main()
// {
//   string fullname,
//          firstname,
//          middlename,
//          lastname;
// 
//   cout << "Enter your first, middle, and last names: ";
//   cin >> fullname;
// 
//   firstname = fullname.substr(0,fullname.find(' ')); // Get all chars from index 0 to index of first space found.
//   fullname = fullname.substr(fullname.find(' ')+1,fullname.length()); // Set fullname to contain only chars found past the first space.
//   middlename = fullname.substr(0,fullname.find(" "));
//   fullname = fullname.substr(fullname.find(" ")+1,fullname.length());
//   lastname = fullname.substr(0,fullname.length());
//   cout << firstname << " " << middlename << " " << lastname;
// 
//   return 0;
// }
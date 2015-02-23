THIS FILE NEEDS TO BE UPDATED. -7 Oct 2012

WHAT DOES IT DO: This code calculates the signal efficiency for SS2L production
at the LHC in a simplified model with stops and gluinos and an RPV coupling
W ~ UDD. The output is stored in an output file (output.dat) and the code
is written to make it easy to scan over parameters using a bash script.

PartonRPV: parton level signal, unchanged
PartonBGRPV: background from LHE file

This code has been tested to work with Pythia 8.165.

USAGE:

1. Go to the Makefile and modify the first two non-commented lines
    These correspond to the path to your PYTHIA and FASTJET folders respectively
    Everything else in this file should work out-of-the-box.
    You might want to change the compiler/flags. Whatever floats your boat, dude.
   
    
2. In the terminal, compile this program using the command 'make'
    If there are errors, then it's either your fault or my fault.
   
    
3. If successful, make will output instructions for how to use the compiled code.
    In particular, there are 'template' command and spectrum files which
    the program runs. It doesn't read these directly, but instead modifies
    them with a revised stop and gluino mass. It then generates a NEW cmnd and
    spc file (actually, it also generates an intermediate spc file) and feeds
    these new files into Pythia. 
    
    If you want to modify something in the Pythia run, go ahead and modify the
    command file template (or even better, create a new one). For example, the
    current code assumes that we're forcing leptons to decay leptonically. 
    (The code compensates for this with an overall prefactor to the efficiency
    that is sent to the output file. Note that the efficiency that is listed
    upon running does not yet have this prefactor.) You might want to change this
    in the template command file.
    
4. Anyway, go ahead and run the program following the sample call from the
    instructions listed in the makefile. For example, you can run with default
    values without any options:
    
        ./PartonRPV
    
    You can also put in a bunch of options:

        ./PartonRPV 300 800 8 CmndShort.cmnd output.dat template.spc
    
    These correspond to (in order):
        stop mass
        gluino mass
        signal region (for comparison with CMS SUS-12-017)
        command file template (CmndShort.cmnd only runs 100 events)
        output file
        spectrum file template
    
    There are default values for each of these, but if you specify any options,
    you have to make sure all possible options before are also specified. In 
    other words, if you only wanted to change the signal region, you would still
    have to specify the stop and gluino masses, since
    
        ./PartonRPV 8
        
    would be interpreted as setting the stop mass to 8.
    
5. Scanning with a batch script: this was the raison d'etre for this code. 
    This should be fairly straightforward since you can just scan over the
    options for the program. 
    
    
Good scanning,
Flip, Sept 2012
    
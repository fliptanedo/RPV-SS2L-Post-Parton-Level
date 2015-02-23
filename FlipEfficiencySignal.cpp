/******************************************************************************** 
*   FlipEfficiencySignal.cpp by Flip Tanedo (pt267@cornell.edu)                 *
*   Code for FlipEfficiencySignal.cc, Aug 2012                                  *
*   Modified for BG, b-tagging Oct 2012                                         *
*   Contains routine for calculating signal efficiency                          *
*   To do: fix this all up so there's just one signal efficiency function       *
********************************************************************************/

#include "FlipEfficiency.h"

double signal_efficiency(
    string command_file,                    // Pythia data
    vector< pair<string, int> > &counts,    // intermediate data (for checking)
    int iSR                                 // Signal Region #
    ){
    // For a given parameter space point, outputs the signal efficiency
    // Fills the vector with a list of intermediate counts    
    
    /****************************************************************************
    *   SET UP GENERATION                                                       *
    ****************************************************************************/
    
    Pythia8::Pythia pythia;                 // Declare Pythia object
    Pythia8::Event& event = pythia.event;   // Declare event as a shortcut
    pythia.readFile(command_file);          // Read in command file

    int nEvent = pythia.mode("Main:numberOfEvents");
    int nAbort = pythia.mode("Main:timesAllowErrors");

    pythia.init();
    
    
    // SIGNAL REGIONS
    vector<signalregion> signal_region;    // as defined in SUS-12-017 Table 1
    fill_signalregions(signal_region);     // fills data from above paper

    
    
    /****************************************************************************
    *   SET UP COUNTERS FOR SANITY CHECK COUNTS                                 *
    ****************************************************************************/
    
    int nGenerated  = 0; // # generated events 
    int nKinematic  = 0; // # events that pass kinematic cuts on leptons
    // int nLepSelect  = 0; // # events that pass lepton selection efficiencies
    int nLepID      = 0; // # events that pass lepton ID efficiencies
    int nLepIso     = 0; // # events that pass lepton isolation efficiencies
    int nbjetSelect = 0; // # events that pass bJet selection efficiencies
    int nDilepton   = 0; // # events that pass dilepton requirement
    int nDilepTrig  = 0; // # events that pass dilep req & trig efficiency
    int nSS2L       = 0; // # events that pass same sign leptons requirement
    int nJets       = 0; // # events with mininum number of jets
    int nbJets      = 0; // # events with minimum number of tagged b jets
    int nMET        = 0; // # events that pass minimum MET requirement
    int nHT         = 0; // # events that pass minimum HT requirement
    int nCharge     = 0; // # events that pass ++ or -- requirement
    int nPassed     = 0; // # events that passed all cuts
    
    
    /****************************************************************************
    *   GENERATE EVENTS                                                         *
    ****************************************************************************/
    
    int iAbort = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) { // event loop
        
        // Quit if too many aborts
        if (!pythia.next()) {                   // if no new event
            if (++iAbort < nAbort) continue;    // if not over abort limit
            cout << " Event generation aborted prematurely, owing to error!\n"; 
            break;
        } // End of 'if no new event'        
        
            
    
        /************************************************************************
        * SET UP EVENT DATA FOR LATER INSPECTION                                *
        ************************************************************************/
        
        vector< pair<int, fastjet::PseudoJet> > preleptons;     // generated
        vector< pair<int, fastjet::PseudoJet> > prepartons;
        fastjet::PseudoJet METvec (0.0, 0.0, 0.0, 0.0);
        double MET (0.0);
        double HT (0.0);
        
        
        
        /************************************************************************
        * LOOP THROUGH EVENT PARTICLES                                          *
        ************************************************************************/


        // LOOP THROUGH THE TOTAL EVENT
        for (int iPart = 0; iPart < event.size(); iPart++){
            
            // Skip things that we can't see
            if (!event[iPart].isFinal()) continue;
            if (!event[iPart].isVisible()) continue;   
            if (abs(event[iPart].eta()) >= 5.0) continue;
            
            fastjet::PseudoJet momentum(event[iPart].px(),
                                        event[iPart].py(),
                                        event[iPart].pz(),
                                        event[iPart].e());
            
            // Fill Missing ET vector                            
            METvec -= momentum;
            
            // Keep track of leptons
            if ((abs(event[iPart].id()) == 11) || 
                (abs(event[iPart].id()) == 13) ) {
                //
                pair<int,fastjet::PseudoJet> cur_lept;
                cur_lept.first = event[iPart].id();     // PDG code
                cur_lept.second = momentum;             // 4-vector
                preleptons.push_back(cur_lept);            // put in list
                continue;
            } // End "if this is an identfiable lepton"  
            
            
            // Anything left is a parton
            pair<int,fastjet::PseudoJet> cur_parton;
            cur_parton.first = event[iPart].id();
            cur_parton.second = momentum;
            prepartons.push_back(cur_parton); 
            HT += momentum.pt();
            
            
        } // End loop through event particles
        
        // Increment counter
        nGenerated++;        
        
        
        
        /************************************************************************
        * IMPOSE CUTS                                                           *
        ************************************************************************/
        
        // Lepton Kinematics
        // -----------------
        
        vector< pair<int, fastjet::PseudoJet> > leptons_kin;
        for(unsigned int iLep = 0; iLep < preleptons.size(); iLep++){
            if (lepton_kinematic_cut(preleptons[iLep]))
                leptons_kin.push_back(preleptons[iLep]);
        } // end for loop over leptons
        
        if (leptons_kin.size() > 1) nKinematic++;
        else continue;
        
        
        
        // Parton Kinematics
        // -----------------
        
        vector< pair<int, fastjet::PseudoJet> > partons;
        for(unsigned int iJet = 0; iJet < prepartons.size(); iJet++){
            if (jet_kinematic_cut(prepartons[iJet]))
                partons.push_back(prepartons[iJet]);
        } // end for loop over partons
                

        
        // Selection Cuts
        // --------------

        // // UNUSED: separate selection eff = ID eff * isolation eff
        // vector< pair<int, fastjet::PseudoJet> > leptons;
        // for(unsigned int iLep = 0; iLep < leptons_kin.size(); iLep++){
        //     if (lepton_selection_cut(leptons_kin[iLep]))
        //         leptons.push_back(leptons_kin[iLep]);
        // } // end for loop over leptons
        // 
        // if (leptons.size() > 1) nLepSelect++;
        // else continue;
        
        vector< pair<int, fastjet::PseudoJet> > leptons_ID;
        for(unsigned int iLep = 0; iLep < leptons_kin.size(); iLep++){
            if (lepton_ID_eff(leptons_kin[iLep]))
                leptons_ID.push_back(leptons_kin[iLep]);
        } // end for loop over leptons
        
        if (leptons_ID.size() > 1) nLepID++;
        else continue;
        
        
        
        vector< pair<int, fastjet::PseudoJet> > leptons;
        for(unsigned int iLep = 0; iLep < leptons_ID.size(); iLep++){
            if (lepton_iso_eff(leptons_ID[iLep], partons))
                leptons.push_back(leptons_ID[iLep]);
        } // end for loop over leptons
        
        if (leptons.size() > 1) nLepIso++;
        else continue;
        
        
        
        vector< pair<int, fastjet::PseudoJet> > bJets;
        for(unsigned int iJet = 0; iJet < partons.size(); iJet++)
            if( ( abs(partons[iJet].first)==5 ) && 
                b_selection_efficiency(partons[iJet]) )
                bJets.push_back(partons[iJet]);

        if (bJets.size() > 1) nbjetSelect++;
        else continue;
        
        
        
        // Event cuts
        // ----------
        
        // Exactly two leptons
        if (leptons.size() != 2) continue;
        else nDilepton++;
        
        // Trigger efficiency for dilepton
        if (lepton_trig_efficiency(leptons)) continue;
        else nDilepTrig++;
        
        // Same-sign dileptons
        if (leptons[0].first/abs(leptons[0].first) != 
            leptons[1].first/abs(leptons[1].first)) continue;
        else nSS2L++;
        
        
        
        // Signal region cuts: from input
        // ------------------------------
        
        if (partons.size() < signal_region[iSR].minJets) continue;        
        else nJets++;
        
        if (bJets.size() < signal_region[iSR].minbJets) continue;
        else nbJets++;
        
        MET = METvec.pt(); 
        // if (MET < signal_region[iSR].minMET) continue;
        if (!METefficiency(MET,signal_region[iSR].minMET)) continue;
        else nMET++;
        
        // if (HT < signal_region[iSR].minHT) continue;
        if (!HTefficiency(HT,signal_region[iSR].minHT)) continue;
        else nHT++;
        
        bool minmin = (leptons[0].first > 0) && signal_region[iSR].minusminus;
        bool pluplu = (leptons[0].first < 0) && signal_region[iSR].plusplus;
        
        if (!(minmin || pluplu)) continue;
        else nCharge++;
        
        // Made it this far? YOU PASS
        nPassed++;
        
    } // end for loop, going through Events
    
    
    
    // Fill counts
    // -----------
    fill_vector(counts, "Generated events \t", nGenerated);
    fill_vector(counts, ">1 lep. kin. cuts\t", nKinematic);
    // fill_vector(counts, ">1 lep. sel. eff.\t", nLepSelect);
    fill_vector(counts, ">1 lep. ID. eff.\t", nLepID);
    fill_vector(counts, ">1 lep. Iso. eff.\t", nLepIso);
    fill_vector(counts, ">1 bjets tagged \t", nbjetSelect);
    fill_vector(counts, "exactly two leptons \t", nDilepton);
    fill_vector(counts, "triggered two leptons \t", nDilepTrig);
    fill_vector(counts, "same sign dileptons \t", nSS2L);
    
    // The following cuts depend on the signal region, so we have to
    //  "dynamically" generate their labels
    
    stringstream nJetComment;
    nJetComment << "at least " << signal_region[iSR].minJets << " jets \t";
    fill_vector(counts, nJetComment.str(), nJets);
    
    stringstream nbJetComment;
    nbJetComment << "at least " << signal_region[iSR].minbJets << " b jets \t";
    fill_vector(counts, nbJetComment.str(), nbJets);
    
    stringstream nMETComment;
    nMETComment << "at least " << signal_region[iSR].minMET << " GeV MET \t";
    fill_vector(counts, nMETComment.str(), nMET);
    
    stringstream HTComment;
    HTComment << "at least " << signal_region[iSR].minHT << " GeV HT \t";
    fill_vector(counts, HTComment.str(), nHT);
    
    stringstream nChargeComment;
    if ( signal_region[iSR].minusminus && !signal_region[iSR].plusplus)
        nChargeComment << "only -- leptons \t";
    else if ( !signal_region[iSR].minusminus && signal_region[iSR].plusplus)
        nChargeComment << "only ++ leptons \t";
    else if ( signal_region[iSR].minusminus && signal_region[iSR].plusplus)
        nChargeComment << "either ++ or -- leptons";
    else nChargeComment << "You fucked up, neither ++ or -- leptons ";
    
    fill_vector(counts, nChargeComment.str(), nCharge);
    
    
    
    
    cout << endl << endl << "STOP: " << pythia.particleData.m0(1000006) << endl;
    cout << "GLUINO: " << pythia.particleData.m0(1000021) << endl;
    cout << "Signal Region " << iSR << endl; 
    cout << "Efficiency: " << double(nPassed) / double(nEvent) << endl;

    return double(nPassed) / double(nEvent); 
    

    
} // end void signal_efficiency(...)









double BG_efficiency(
    Pythia8::Pythia& pythia,                // pythia object
    vector< pair<string, int> > &counts,    // intermediate data (for checking)
    int iSR,                                // Signal Region #
    int nEvent                              // # events
    ){
    // For a given parameter space point, outputs the signal efficiency
    // Fills the vector with a list of intermediate counts    
    
    Pythia8::Event& event = pythia.event;   // Declare event as a shortcut
    // int nEvent = pythia.mode("Main:numberOfEvents"); // THIS IS A PROBLEM
    int nAbort = pythia.mode("Main:timesAllowErrors");

    
    // SIGNAL REGIONS
    vector<signalregion> signal_region;    // as defined in SUS-12-017 Table 1
    fill_signalregions(signal_region);     // fills data from above paper
    
    
    
    /****************************************************************************
    *   SET UP COUNTERS FOR SANITY CHECK COUNTS                                 *
    ****************************************************************************/
    
    int nGenerated  = 0; // # generated events 
    int nKinematic  = 0; // # events that pass kinematic cuts on leptons
    // int nLepSelect  = 0; // # events that pass lepton selection efficiencies
    int nLepID      = 0; // # events that pass lepton ID efficiencies
    int nLepIso     = 0; // # events that pass lepton isolation efficiencies
    int nbjetSelect = 0; // # events that pass bJet selection efficiencies
    int nDilepton   = 0; // # events that pass dilepton requirement
    int nDilepTrig  = 0; // # events that pass dilep req & trig efficiency
    int nSS2L       = 0; // # events that pass same sign leptons requirement
    int nJets       = 0; // # events with mininum number of jets
    int nbJets      = 0; // # events with minimum number of tagged b jets
    int nMET        = 0; // # events that pass minimum MET requirement
    int nHT         = 0; // # events that pass minimum HT requirement
    int nCharge     = 0; // # events that pass ++ or -- requirement
    int nPassed     = 0; // # events that passed all cuts
    
    
    /****************************************************************************
    *   GENERATE EVENTS                                                         *
    ****************************************************************************/
    
    int iAbort = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) { // event loop
        
        // Quit if too many aborts
        if (!pythia.next()) {                   // if no new event
            if (++iAbort < nAbort) continue;    // if not over abort limit
            cout << " Event generation aborted prematurely, owing to error!\n"; 
            break;
        } // End of 'if no new event'        
        
            
    
        /************************************************************************
        * SET UP EVENT DATA FOR LATER INSPECTION                                *
        ************************************************************************/
        
        vector< pair<int, fastjet::PseudoJet> > preleptons;     // generated
        vector< pair<int, fastjet::PseudoJet> > prepartons;
        fastjet::PseudoJet METvec (0.0, 0.0, 0.0, 0.0);
        double MET (0.0);
        double HT (0.0);
        
        
        
        /************************************************************************
        * LOOP THROUGH EVENT PARTICLES                                          *
        ************************************************************************/
        
        for (int iPart = 0; iPart < event.size(); iPart++){
            
            // Skip things that we can't see
            if (!event[iPart].isFinal()) continue;
            if (!event[iPart].isVisible()) continue;   
            if (abs(event[iPart].eta()) >= 5.0) continue;
            
            fastjet::PseudoJet momentum(event[iPart].px(),
                                        event[iPart].py(),
                                        event[iPart].pz(),
                                        event[iPart].e());
            
            // Fill Missing ET vector                            
            METvec -= momentum;
            
            // Keep track of leptons
            if ((abs(event[iPart].id()) == 11) || 
                (abs(event[iPart].id()) == 13) ) {
                //
                pair<int,fastjet::PseudoJet> cur_lept;
                cur_lept.first = event[iPart].id();     // PDG code
                cur_lept.second = momentum;             // 4-vector
                preleptons.push_back(cur_lept);            // put in list
                continue;
            } // End "if this is an identfiable lepton"  
            
            
            // Anything left is a parton
            pair<int,fastjet::PseudoJet> cur_parton;
            cur_parton.first = event[iPart].id();
            cur_parton.second = momentum;
            prepartons.push_back(cur_parton); 
            HT += momentum.pt();
            
            
        } // End loop through event particles
        
    
        // Increment counter
        nGenerated++;        
        
        
        
        /************************************************************************
        * IMPOSE CUTS                                                           *
        ************************************************************************/
        
        // Lepton Kinematics
        // -----------------
        
        vector< pair<int, fastjet::PseudoJet> > leptons_kin;
        for(unsigned int iLep = 0; iLep < preleptons.size(); iLep++){
            if (lepton_kinematic_cut(preleptons[iLep]))
                leptons_kin.push_back(preleptons[iLep]);
        } // end for loop over leptons
        
        if (leptons_kin.size() > 1) nKinematic++;
        else continue;
        
        
        
        // Parton Kinematics
        // -----------------
        
        vector< pair<int, fastjet::PseudoJet> > partons;
        for(unsigned int iJet = 0; iJet < prepartons.size(); iJet++){
            if (jet_kinematic_cut(prepartons[iJet]))
                partons.push_back(prepartons[iJet]);
        } // end for loop over partons
                
    
        
        // Selection Cuts
        // --------------
    
        // // UNUSED: separate selection eff = ID eff * isolation eff
        // vector< pair<int, fastjet::PseudoJet> > leptons;
        // for(unsigned int iLep = 0; iLep < leptons_kin.size(); iLep++){
        //     if (lepton_selection_cut(leptons_kin[iLep]))
        //         leptons.push_back(leptons_kin[iLep]);
        // } // end for loop over leptons
        // 
        // if (leptons.size() > 1) nLepSelect++;
        // else continue;
        
        vector< pair<int, fastjet::PseudoJet> > leptons_ID;
        for(unsigned int iLep = 0; iLep < leptons_kin.size(); iLep++){
            if (lepton_ID_eff(leptons_kin[iLep]))
                leptons_ID.push_back(leptons_kin[iLep]);
        } // end for loop over leptons
        
        if (leptons_ID.size() > 1) nLepID++;
        else continue;
        
        
        
        vector< pair<int, fastjet::PseudoJet> > leptons;
        for(unsigned int iLep = 0; iLep < leptons_ID.size(); iLep++){
            if (lepton_iso_eff(leptons_ID[iLep], partons))
                leptons.push_back(leptons_ID[iLep]);
        } // end for loop over leptons
        
        if (leptons.size() > 1) nLepIso++;
        else continue;
        
        
        
        vector< pair<int, fastjet::PseudoJet> > bJets;
        for(unsigned int iJet = 0; iJet < partons.size(); iJet++)
            if( ( abs(partons[iJet].first)==5 ) && 
                b_selection_efficiency(partons[iJet]) )
                bJets.push_back(partons[iJet]);
    
        if (bJets.size() > 1) nbjetSelect++;
        else continue;
        
        
        
        // Event cuts
        // ----------
        
        // Exactly two leptons
        if (leptons.size() != 2) continue;
        else nDilepton++;
        
        // Trigger efficiency for dilepton
        if (lepton_trig_efficiency(leptons)) continue;
        else nDilepTrig++;
        
        // Same-sign dileptons
        if (leptons[0].first/abs(leptons[0].first) != 
            leptons[1].first/abs(leptons[1].first)) continue;
        else nSS2L++;
        
        
        
        // Signal region cuts: from input
        // ------------------------------
        
        if (partons.size() < signal_region[iSR].minJets) continue;        
        else nJets++;
        
        if (bJets.size() < signal_region[iSR].minbJets) continue;
        else nbJets++;
        
        MET = METvec.pt(); 
        // if (MET < signal_region[iSR].minMET) continue;
        if (!METefficiency(MET,signal_region[iSR].minMET)) continue;
        else nMET++;
        
        // if (HT < signal_region[iSR].minHT) continue;
        if (!HTefficiency(HT,signal_region[iSR].minHT)) continue;
        else nHT++;
        
        bool minmin = (leptons[0].first > 0) && signal_region[iSR].minusminus;
        bool pluplu = (leptons[0].first < 0) && signal_region[iSR].plusplus;
        
        if (!(minmin || pluplu)) continue;
        else nCharge++;
        
        // Made it this far? YOU PASS
        nPassed++;
        
    } // end for loop, going through Events
    
    
    
    // Fill counts
    // -----------
    fill_vector(counts, "Generated events \t", nGenerated);
    fill_vector(counts, ">1 lep. kin. cuts\t", nKinematic);
    // fill_vector(counts, ">1 lep. sel. eff.\t", nLepSelect);
    fill_vector(counts, ">1 lep. ID. eff.\t", nLepID);
    fill_vector(counts, ">1 lep. Iso. eff.\t", nLepIso);
    fill_vector(counts, ">1 bjets tagged \t", nbjetSelect);
    fill_vector(counts, "exactly two leptons \t", nDilepton);
    fill_vector(counts, "triggered two leptons \t", nDilepTrig);
    fill_vector(counts, "same sign dileptons \t", nSS2L);
    
    // The following cuts depend on the signal region, so we have to
    //  "dynamically" generate their labels
    
    stringstream nJetComment;
    nJetComment << "at least " << signal_region[iSR].minJets << " jets \t";
    fill_vector(counts, nJetComment.str(), nJets);
    
    stringstream nbJetComment;
    nbJetComment << "at least " << signal_region[iSR].minbJets << " b jets \t";
    fill_vector(counts, nbJetComment.str(), nbJets);
    
    stringstream nMETComment;
    nMETComment << "at least " << signal_region[iSR].minMET << " GeV MET \t";
    fill_vector(counts, nMETComment.str(), nMET);
    
    stringstream HTComment;
    HTComment << "at least " << signal_region[iSR].minHT << " GeV HT \t";
    fill_vector(counts, HTComment.str(), nHT);
    
    stringstream nChargeComment;
    if ( signal_region[iSR].minusminus && !signal_region[iSR].plusplus)
        nChargeComment << "only -- leptons \t";
    else if ( !signal_region[iSR].minusminus && signal_region[iSR].plusplus)
        nChargeComment << "only ++ leptons \t";
    else if ( signal_region[iSR].minusminus && signal_region[iSR].plusplus)
        nChargeComment << "either ++ or -- leptons";
    else nChargeComment << "You fucked up, neither ++ or -- leptons ";
    
    fill_vector(counts, nChargeComment.str(), nCharge);
    
    
    
    
    cout << endl << endl << "STOP: " << pythia.particleData.m0(1000006) << endl;
    cout << "GLUINO: " << pythia.particleData.m0(1000021) << endl;
    cout << "Signal Region " << iSR << endl; 
    cout << "Efficiency: " << double(nPassed) / double(nEvent) << endl;
    
    return double(nPassed) / double(nEvent); 
    
    
} // end void signal_efficiency(...)



double signal_efficiency_b(
    string command_file,                    // Pythia data
    vector< pair<string, int> > &counts,    // intermediate data (for checking)
    int iSR                                 // Signal Region #
    ){
    // For a given parameter space point, outputs the signal efficiency
    // Fills the vector with a list of intermediate counts    
    
    /****************************************************************************
    *   SET UP GENERATION                                                       *
    ****************************************************************************/
    
    Pythia8::Pythia pythia;                 // Declare Pythia object
    Pythia8::Event& event = pythia.event;   // Declare event as a shortcut
    Pythia8::Event& process = pythia.process; 
                                            // only hard process, for b-tag
    
    pythia.readFile(command_file);          // Read in command file

    int nEvent = pythia.mode("Main:numberOfEvents");
    int nAbort = pythia.mode("Main:timesAllowErrors");

    pythia.init();
    
    
    // SIGNAL REGIONS
    vector<signalregion> signal_region;    // as defined in SUS-12-017 Table 1
    fill_signalregions(signal_region);     // fills data from above paper

    
    
    /****************************************************************************
    *   SET UP COUNTERS FOR SANITY CHECK COUNTS                                 *
    ****************************************************************************/
    
    int nGenerated  = 0; // # generated events 
    int nKinematic  = 0; // # events that pass kinematic cuts on leptons
    // int nLepSelect  = 0; // # events that pass lepton selection efficiencies
    int nLepID      = 0; // # events that pass lepton ID efficiencies
    int nLepIso     = 0; // # events that pass lepton isolation efficiencies
    int nbjetSelect = 0; // # events that pass bJet selection efficiencies
    int nDilepton   = 0; // # events that pass dilepton requirement
    int nDilepTrig  = 0; // # events that pass dilep req & trig efficiency
    int nSS2L       = 0; // # events that pass same sign leptons requirement
    int nJets       = 0; // # events with mininum number of jets
    int nbJets      = 0; // # events with minimum number of tagged b jets
    int nMET        = 0; // # events that pass minimum MET requirement
    int nHT         = 0; // # events that pass minimum HT requirement
    int nCharge     = 0; // # events that pass ++ or -- requirement
    int nPassed     = 0; // # events that passed all cuts
    
    
    /****************************************************************************
    *   GENERATE EVENTS                                                         *
    ****************************************************************************/
    
    int iAbort = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) { // event loop
        
        // Quit if too many aborts
        if (!pythia.next()) {                   // if no new event
            if (++iAbort < nAbort) continue;    // if not over abort limit
            cout << " Event generation aborted prematurely, owing to error!\n"; 
            break;
        } // End of 'if no new event'        
        
            
    
        /************************************************************************
        * SET UP EVENT DATA FOR LATER INSPECTION                                *
        ************************************************************************/
        
        vector< pair<int, fastjet::PseudoJet> > preleptons; // from event
        vector< pair<int, fastjet::PseudoJet> > prepartons; // from event
        vector< pair<int, fastjet::PseudoJet> > bpartons;   // from process
        fastjet::PseudoJet METvec (0.0, 0.0, 0.0, 0.0);
        double MET (0.0);
        double HT (0.0);
        
        
        
        /************************************************************************
        * LOOP THROUGH EVENT PARTICLES                                          *
        ************************************************************************/
        
        // LOOP THROUGH TOTAL EVENT PARTICLES
        // ----------------------------------
        // To be distinguished from the loop over the hard
        // event ("process") particles.
        //
        for (int iPart = 0; iPart < event.size(); iPart++){
            
            // Skip things that we can't see
            if (!event[iPart].isFinal()) continue;
            if (!event[iPart].isVisible()) continue;   
            if (abs(event[iPart].eta()) >= 5.0) continue;
            
            fastjet::PseudoJet momentum(event[iPart].px(),
                                        event[iPart].py(),
                                        event[iPart].pz(),
                                        event[iPart].e());
            
            // Fill Missing ET vector                            
            METvec -= momentum;
            
            // Keep track of leptons
            if ((abs(event[iPart].id()) == 11) || 
                (abs(event[iPart].id()) == 13) ) {
                //
                pair<int,fastjet::PseudoJet> cur_lept;
                cur_lept.first = event[iPart].id();     // PDG code
                cur_lept.second = momentum;             // 4-vector
                preleptons.push_back(cur_lept);            // put in list
                continue;
            } // End "if this is an identfiable lepton"  
            
            
            // Anything left is a parton
            pair<int,fastjet::PseudoJet> cur_parton;
            cur_parton.first = event[iPart].id();
            cur_parton.second = momentum;
            prepartons.push_back(cur_parton); 
            HT += momentum.pt();
            
            
        } // End loop through event particles
        
        // LOOP THROUGH HARD EVENT PARTICLES
        // ---------------------------------
        // Use pythia.process to accesses only the hard scattering event
        // We will use this for b-tagging since we have the parton-level
        // b-tagging efficiency. 
        //
        for (int iPart = 0; iPart < process.size(); iPart++){
            
            if (!process[iPart].isFinal()) continue;
            if (!process[iPart].isVisible()) continue;   
            if (abs(process[iPart].eta()) >= 5.0) continue;
            if (!abs(process[iPart].id())==5) continue;     //only bjets
            
            fastjet::PseudoJet momentum(process[iPart].px(),
                                        process[iPart].py(),
                                        process[iPart].pz(),
                                        process[iPart].e());
            
            // bpartons
            pair<int,fastjet::PseudoJet> cur_bparton;
            cur_bparton.first = process[iPart].id();
            cur_bparton.second = momentum;
            bpartons.push_back(cur_bparton); 
    
        } // End loop through process particles
        
    
        // Increment counter
        nGenerated++;        
        
        
        
        /************************************************************************
        * IMPOSE CUTS                                                           *
        ************************************************************************/
        
        // LEPTON KINEMATICS
        // Check that allowed leptons satisfy kinematic cuts 
        // -------------------------------------------------
        
        vector< pair<int, fastjet::PseudoJet> > leptons_kin;
        for(unsigned int iLep = 0; iLep < preleptons.size(); iLep++){
            if (lepton_kinematic_cut(preleptons[iLep]))
                leptons_kin.push_back(preleptons[iLep]);
        } // end for loop over leptons
        
        if (leptons_kin.size() > 1) nKinematic++;
        else continue;
        
        
        
        // PARTON KINEMATICS
        // Check that allowed partons satisfy kinematic cuts
        // -------------------------------------------------
        
        vector< pair<int, fastjet::PseudoJet> > partons;
        for(unsigned int iJet = 0; iJet < prepartons.size(); iJet++){
            if (jet_kinematic_cut(prepartons[iJet]))
                partons.push_back(prepartons[iJet]);
        } // end for loop over partons
                

        
        // SELECTION EFFICIENCIES
        // "selection cuts"
        // (Roll the dice)
        // ----------------------

        // // Instead of using this, we separate into ID * iso efficiency
        // // UNUSED: separate selection eff = ID eff * isolation eff
        //
        // vector< pair<int, fastjet::PseudoJet> > leptons;
        // for(unsigned int iLep = 0; iLep < leptons_kin.size(); iLep++){
        //     if (lepton_selection_cut(leptons_kin[iLep]))
        //         leptons.push_back(leptons_kin[iLep]);
        // } // end for loop over leptons
        // 
        // if (leptons.size() > 1) nLepSelect++;
        // else continue;
        
        vector< pair<int, fastjet::PseudoJet> > leptons_ID;
        for(unsigned int iLep = 0; iLep < leptons_kin.size(); iLep++){
            if (lepton_ID_eff(leptons_kin[iLep]))
                leptons_ID.push_back(leptons_kin[iLep]);
        } // end for loop over leptons
        
        if (leptons_ID.size() > 1) nLepID++;
        else continue;
        
        
        
        vector< pair<int, fastjet::PseudoJet> > leptons;
        for(unsigned int iLep = 0; iLep < leptons_ID.size(); iLep++){
            if (lepton_iso_eff(leptons_ID[iLep], partons))
                leptons.push_back(leptons_ID[iLep]);
        } // end for loop over leptons
        
        if (leptons.size() > 1) nLepIso++;
        else continue;
        
        
        // // bjet tagging at parton level (bpartons)
        //         vector< pair<int, fastjet::PseudoJet> > bJets;
        //         for(unsigned int iJet = 0; iJet < partons.size(); iJet++)
        //             if( ( abs(partons[iJet].first)==5 ) && 
        //                 b_selection_efficiency(partons[iJet]) )
        //                 bJets.push_back(partons[iJet]);
        // 
        //         if (bJets.size() > 1) nbjetSelect++;
        //         else continue;
        
        // bjet tagging at parton level (bpartons)
        // This counts the number of b jets based on
        // the parton level selection efficiency given by CMS
        vector< pair<int, fastjet::PseudoJet> > bJets;
        for(unsigned int iJet = 0; iJet < bpartons.size(); iJet++)
            if( b_selection_efficiency(bpartons[iJet]) )
                bJets.push_back(bpartons[iJet]);

        if (bJets.size() > 1) nbjetSelect++;
        else continue;
        
        
        // Event cuts
        // ----------
        
        // Exactly two leptons
        if (leptons.size() != 2) continue;
        else nDilepton++;
        
        // Trigger efficiency for dilepton
        if (lepton_trig_efficiency(leptons)) continue;
        else nDilepTrig++;
        
        // Same-sign dileptons
        if (leptons[0].first/abs(leptons[0].first) != 
            leptons[1].first/abs(leptons[1].first)) continue;
        else nSS2L++;
        
        
        
        // Signal region cuts: from input
        // ------------------------------
        
        if (partons.size() < signal_region[iSR].minJets) continue;        
        else nJets++;
        
        if (bJets.size() < signal_region[iSR].minbJets) continue;
        else nbJets++;
        
        MET = METvec.pt(); 
        // if (MET < signal_region[iSR].minMET) continue;
        if (!METefficiency(MET,signal_region[iSR].minMET)) continue;
        else nMET++;
        
        // if (HT < signal_region[iSR].minHT) continue;
        if (!HTefficiency(HT,signal_region[iSR].minHT)) continue;
        else nHT++;
        
        bool minmin = (leptons[0].first > 0) && signal_region[iSR].minusminus;
        bool pluplu = (leptons[0].first < 0) && signal_region[iSR].plusplus;
        
        if (!(minmin || pluplu)) continue;
        else nCharge++;
        
        // Made it this far? YOU PASS
        nPassed++;
        
    } // end for loop, going through Events
    
    
    
    // Fill counts
    // -----------
    fill_vector(counts, "Generated events \t", nGenerated);
    fill_vector(counts, ">1 lep. kin. cuts\t", nKinematic);
    // fill_vector(counts, ">1 lep. sel. eff.\t", nLepSelect);
    fill_vector(counts, ">1 lep. ID. eff.\t", nLepID);
    fill_vector(counts, ">1 lep. Iso. eff.\t", nLepIso);
    fill_vector(counts, ">1 bjets tagged \t", nbjetSelect);
    fill_vector(counts, "exactly two leptons \t", nDilepton);
    fill_vector(counts, "triggered two leptons \t", nDilepTrig);
    fill_vector(counts, "same sign dileptons \t", nSS2L);
    
    // The following cuts depend on the signal region, so we have to
    //  "dynamically" generate their labels
    
    stringstream nJetComment;
    nJetComment << "at least " << signal_region[iSR].minJets << " jets \t";
    fill_vector(counts, nJetComment.str(), nJets);
    
    stringstream nbJetComment;
    nbJetComment << "at least " << signal_region[iSR].minbJets << " b jets \t";
    fill_vector(counts, nbJetComment.str(), nbJets);
    
    stringstream nMETComment;
    nMETComment << "at least " << signal_region[iSR].minMET << " GeV MET \t";
    fill_vector(counts, nMETComment.str(), nMET);
    
    stringstream HTComment;
    HTComment << "at least " << signal_region[iSR].minHT << " GeV HT \t";
    fill_vector(counts, HTComment.str(), nHT);
    
    stringstream nChargeComment;
    if ( signal_region[iSR].minusminus && !signal_region[iSR].plusplus)
        nChargeComment << "only -- leptons \t";
    else if ( !signal_region[iSR].minusminus && signal_region[iSR].plusplus)
        nChargeComment << "only ++ leptons \t";
    else if ( signal_region[iSR].minusminus && signal_region[iSR].plusplus)
        nChargeComment << "either ++ or -- leptons";
    else nChargeComment << "You fucked up, neither ++ or -- leptons ";
    
    fill_vector(counts, nChargeComment.str(), nCharge);
    
    
    
    
    cout << endl << endl << "STOP: " << pythia.particleData.m0(1000006) << endl;
    cout << "GLUINO: " << pythia.particleData.m0(1000021) << endl;
    cout << "Signal Region " << iSR << endl; 
    cout << "Efficiency: " << double(nPassed) / double(nEvent) << endl;

    return double(nPassed) / double(nEvent); 
    

    
} // end void signal_efficiency(...)



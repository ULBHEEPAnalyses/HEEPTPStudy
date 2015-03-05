// main41.cc is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Mikhail Kirsanov, Mikhail.Kirsanov@cern.ch, based on main01.cc.
// This program illustrates how HepMC can be interfaced to Pythia8.
// HepMC events are output to the ZPrime_TP.dat file.

// WARNING: typically one needs 25 MB/100 events at the LHC.
// Therefore large event samples may be impractical.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

using namespace Pythia8;

#include <iostream>

int main(int argc, char* argv[]) {
  enum ZprimeTypes { kZprimeSSM, kZprimeI } ;
  
  char* mass_string    = "2500" ;
  char* nEvents_string = "10" ;
  char* model_string   = "0" ;
  char* outputFile     = "output.dat" ;
  if(argc>1) mass_string    = argv[1] ;
  if(argc>2) nEvents_string = argv[2] ;
  if(argc>3) model_string   = argv[3] ;
  if(argc>4) outputFile     = argv[4] ;
  
  int model = kZprimeSSM ;
  int model_int = atoi(model_string) ;
  if(model_int==1) model = kZprimeI ;
  int nEvents = atoi(nEvents_string) ;

  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC::Pythia8ToHepMC ToHepMC;
  
  // Specify file where HepMC events will be stored.
  string model_name = "generic" ;
  if(model==kZprimeI) model_name = "ZprimeI" ;
  char filename[500] ;
  sprintf(filename, "%s", outputFile) ;
  std::cout << filename << std::endl ;
  HepMC::IO_GenEvent ascii_io(filename, std::ios::out);
  
  // Generator. Process selection. LHC initializatio
  Pythia pythia;
  pythia.readString("Beams:eCM = 14000.");
  
  // Turn on process, set mass
  pythia.readString("NewGaugeBoson:ffbar2gmZZprime  = on ");
  pythia.readString("Zprime:gmZmode = 3");
  char massString[50] ;
  sprintf(massString, "32:m0 = %s", mass_string) ;
  pythia.readString(massString); 
  // Switch off all Z/Z' decays and then switch back on those to light leptons.
  pythia.readString("32:onMode = off"); 
  pythia.readString("32:onIfAny = 11");
  
  // Zprime I couplings
  if(model==kZprimeI){
    pythia.readString("Zprime:vd =0.620752");
    pythia.readString("Zprime:ad =-0.620752");
    pythia.readString("Zprime:vu =0");
    pythia.readString("Zprime:au =0.");
    pythia.readString("Zprime:ve =-0.620752");
    pythia.readString("Zprime:ae =-0.620752");
    pythia.readString("Zprime:vnue =-0.620752");
    pythia.readString("Zprime:anue =-0.620752");
  }
  
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  
  pythia.init();

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) continue;
    
    // Construct new empty HepMC event and fill it.
    // Units will be as chosen for HepMC build; but can be changed
    // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    ToHepMC.fill_next_event( pythia, hepmcevt );
    
    // Write the HepMC event to file. Done with it.
    ascii_io << hepmcevt;
    delete hepmcevt;
    
    // End of event loop. Statistics. Histogram.
  }
  pythia.stat();

  // Done.
  return 0;
}

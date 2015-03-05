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
  char const* mass_lower_string    = "2000" ;
  char const* mass_upper_string    = "13000" ;
  char const* nEvents_string = "10" ;
  char const* outputFile     = "output.dat" ;
  
  if(argc>1) outputFile        = argv[1] ;
  if(argc>2) nEvents_string    = argv[2] ;
  if(argc>3) mass_lower_string = argv[3] ;
  if(argc>4) mass_upper_string = argv[4] ;
  
  int nEvents = atoi(nEvents_string) ;

  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC::Pythia8ToHepMC ToHepMC;
  
  char filename[500] ;
  sprintf(filename, "%s", outputFile) ;
  std::cout << filename << std::endl ;
  HepMC::IO_GenEvent ascii_io(filename, std::ios::out);
  
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  
  // Find number of all final charged particles and fill histogram.
  // Turn on process, set mass
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  // Switch off all Z0 decays and then switch back on those to quarks.
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 11");
  
  char mass_lower[500] ;
  char mass_upper[500] ;
  sprintf(mass_lower, "23:mMin = %s", mass_lower_string) ;
  sprintf(mass_upper, "23:mMax = %s", mass_upper_string) ;
  pythia.readString(mass_lower);
  pythia.readString(mass_upper);
  
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

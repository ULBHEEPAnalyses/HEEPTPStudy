#!/bin/bash

#mass=$1
#nEvents=$2
#model=$3

mass=500
nEvents=10
model=1

# Remake ZPrime_TP if necessary
cd ~/TP/CMSSW_7_2_0/src/pythia8200/share/Pythia8/examples
make ZPrime_TP

# Generate samples in pythia
date
cd ~/TP/CMSSW_7_2_0/src/pythia8200/share/Pythia8/examples
make ZPrime_TP
./ZPrime_TP ${mass} ${nEvents} ${model}

# Process samples in Delphes
date
cd ~/TP/CMSSW_7_2_0/src/Delphes_CMS

echo -n "" > run_Delphes.sh
for phase in PhaseI PhaseII
do
if [ "${phase}" = "PhaseI" ]
  then
    card='JetStudies_Phase_I_50PileUp.tcl'
fi
if [ "${phase}" = "PhaseII" ]
  then
    card='JetStudies_Phase_II_50PileUp_conf4.tcl'
fi

output="Delphes_ZP_${phase}_${mass}_ZprimeI.root"
input="../pythia8200/share/Pythia8/examples/ZPrime_TP_${mass}_ZprimeI.dat"
echo "./DelphesHepMC Cards/${card} ${output} ${input}" >> run_Delphes.sh
echo "date" >> run_Delphes.sh
done
source run_Delphes.sh


# Return home
cd ~/TP/CMSSW_7_2_0/src
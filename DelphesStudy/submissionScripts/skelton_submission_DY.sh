#!/bin/bash

date="20150304"

mkdir /localgrid/aidan/DY/${date}

XXX
YYY
nEvents=10000

cd /localgrid/aidan/TP
source cmsset_default.sh
cd /localgrid/aidan/TP/CMSSW_7_2_0/src
eval `scramv1 runtime -sh`

prefix='/localgrid/aidan/DY/'

# Generate samples in pythia
date
cd /localgrid/aidan/TP/pythia8200/share/Pythia8/examples
pythiaFile="${prefix}/${date}/pythia_DY_m${massLower}_${massUpper}_${job}.dat"
rm ${pythiaFile}
echo "./DY_mass ${pythiaFile} ${nEvents} ${massLower} ${massUpper}"
./DY_mass ${pythiaFile} ${nEvents} ${massLower} ${massUpper}


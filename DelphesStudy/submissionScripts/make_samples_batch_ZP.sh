#!/bin/bash

date="20150120_b"

mkdir /localgrid/aidan/TP/${date}

#mass=$1
#$nEvents=$2
#$model=$3

XXX
YYY
nEvents=2500
model=1

cd /localgrid/aidan/TP
source cmsset_default.sh
cd /localgrid/aidan/TP/CMSSW_7_2_0/src
eval `scramv1 runtime -sh`

prefix='/localgrid/aidan/TP'

# Generate samples in pythia
date
cd /localgrid/aidan/TP/pythia8200/share/Pythia8/examples
make ZPrime_TP
pythiaFile="${prefix}/${date}/pythia_ZP_m${mass}_ZPrimeI_${job}.dat"
rm ${pythiaFile}
echo "./ZPrime_TP ${mass} ${nEvents} ${model} ${pythiaFile} ${job}"
./ZPrime_TP ${mass} ${nEvents} ${model} ${pythiaFile}

# Process samples in Delphes
date
cd /localgrid/aidan/TP/Delphes

runDelphes="${prefix}/${date}/run_Delphes_${mass}_${job}.sh"
echo -n "" > ${runDelphes}
for phase in PhaseI PhaseII
do
  if [ "${phase}" = "PhaseI" ]
  then
      card='JetStudies_Phase_I_50PileUp.tcl'
  fi
  if [ "${phase}" = "PhaseII" ]
  then
    card='JetStudies_Phase_II_140PileUp_conf4.tcl'
  fi

  DelphesFile="${prefix}/${date}/Delphes_ZP_${phase}_${mass}_ZprimeI_${job}.root"
  rm ${DelphesFile}
  echo "${prefix}/Delphes/DelphesHepMC ${prefix}/Delphes/Cards/${card} ${DelphesFile} ${pythiaFile}" >> ${runDelphes}
  echo "date" >> ${runDelphes}
done

source ${runDelphes}

#!/bin/bash

date="20150304"

mkdir /localgrid/aidan/TP/${date}
WWW
XXX
YYY
ZZZ

cd /localgrid/aidan/TP
source cmsset_default.sh
cd /localgrid/aidan/TP/CMSSW_7_2_0/src
eval `scramv1 runtime -sh`

prefix='/localgrid/aidan/TP'

# Generate samples in pythia
date
cd /localgrid/aidan/TP/pythia8200/share/Pythia8/examples
#make DY_mass
pythiaFile="${prefix}/${date}/pythia_m${massLower}_${massUpper}_DY_${job}.dat"
rm ${pythiaFile}
echo "./DY_mass ${pythiaFile} ${nEvents}  ${massLower} ${massUpper} ${job}"
./DY_mass ${pythiaFile} ${nEvents} ${massLower} ${massUpper}

# Process samples in Delphes
date
cd /localgrid/aidan/TP/Delphes

runDelphes="${prefix}/${date}/run_Delphes_DY_${massLower}_${massUpper}_${job}.sh"
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

  DelphesFile="${prefix}/${date}/Delphes_ZP_${phase}_${massLower}_${massUpper}_DY_${job}.root"
  rm ${DelphesFile}
  echo "${prefix}/Delphes/DelphesHepMC ${prefix}/Delphes/Cards/${card} ${DelphesFile} ${pythiaFile}" >> ${runDelphes}
  echo "date" >> ${runDelphes}
done

source ${runDelphes}

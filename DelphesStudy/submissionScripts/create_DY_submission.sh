echo "" > submit_all_DY.sh
date=20150304
nEvents=2500

for job in {1..20}
do
  for massRange in {1..12}
  do
  
    massLower=0
    massUpper=13000
  
    if [ ${massRange} -eq  1 ]
    then
      massLower=120
      massUpper=200
    fi
    if [ ${massRange} -eq  2 ]
    then
      massLower=200
      massUpper=400
    fi
    if [ ${massRange} -eq  3 ]
    then
      massLower=400
      massUpper=800
    fi
    if [ ${massRange} -eq 4 ]
    then
      massLower=800
      massUpper=1400
    fi
    if [ ${massRange} -eq  5 ]
    then
      massLower=1400
      massUpper=2300
    fi
    if [ ${massRange} -eq  6 ]
    then
      massLower=2300
      massUpper=3500
    fi
    if [ ${massRange} -eq  7 ]
    then
      massLower=3500
      massUpper=4500
    fi
    if [ ${massRange} -eq  8 ]
    then
      massLower=4500
      massUpper=6000
    fi
    if [ ${massRange} -eq  9 ]
    then
      massLower=6000
      massUpper=7500
    fi
    if [ ${massRange} -eq 10 ]
    then
      massLower=7500
      massUpper=8500
    fi
    if [ ${massRange} -eq 11 ]
    then
      massLower=8500
      massUpper=9500
    fi
    if [ ${massRange} -eq 12 ]
    then
      massLower=9500
      massUpper=13000
    fi
    
    submit=submit_scripts/DY_batch_${massLower}_${massUpper}_${job}.sh
    cat make_samples_batch_DY.sh | sed "s/WWW/job=${job}/" | sed "s/XXX/massLower=${massLower}/" | sed "s/YYY/massUpper=${massUpper}/" | sed "s/ZZZ/nEvents=${nEvents}/" > ${submit}
    echo "qsub -q localgrid@cream02 -o /localgrid/aidan/TP/${date}/${massLower}_${massUpper}_${job}.stdout -e /localgrid/aidan/TP/${date}/${massLower}_${massUpper}_${job}.stderr ${submit}" >> submit_all_DY.sh
    echo 'sleep 2' >> submit_all_DY.sh
  done
done


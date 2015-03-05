echo "" > submit_all.sh
date=20150120_b

for job in {1..40}
do
  #for mass in 500 1000 1500 2000 2500 3000 4000 4500 5000
  for mass in 3500
  do
    submit=submit_scripts/make_samples_batch_${mass}_${job}.sh
    cat make_samples_batch.sh | sed "s/XXX/mass=${mass}/" | sed "s/YYY/job=${job}/" > ${submit}
    echo "qsub -q localgrid@cream02 -o /localgrid/aidan/TP/${date}/${mass}_${job}.stdout -e /localgrid/aidan/TP/${date}/${mass}_${job}.stderr ${submit}" >> submit_all.sh
    echo 'sleep 5' >> submit_all.sh
  done
done

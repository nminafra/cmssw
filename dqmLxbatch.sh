#!/bin/bash
mkdir DQMRun_$1$2
cp src/sendDQMJob.sh DQMRun_$1$2/
cp src/list_of_datGR.py DQMRun_$1$2/
cp src/dqm_run.py DQMRun_$1$2/
cd DQMRun_$1$2
sed -i "s/FIRSTTHREE/${1}/g" list_of_datGR.py
sed -i "s/LASTTHREE/${2}/g" list_of_datGR.py
sed -i "s/FIRSTTHREE/${1}/g" sendDQMJob.sh
sed -i "s/LASTTHREE/${2}/g" sendDQMJob.sh

perl -pe 's/FILELIST/`python list_of_datGR.py`/e' ../dqm_run.py > dqm_run.py
sed -i "s/FIRSTTHREE/${1}/g" dqm_run.py
sed -i "s/LASTTHREE/${2}/g" dqm_run.py

bsub -q 2nd -J DQM$1$2 < sendDQMJob.sh -o log_DQM$1$2.log -e err_DQM$1$2.log
cd ..

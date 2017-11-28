#!/bin/bash
mkdir AnalyzerRun_$1$2
cp /afs/cern.ch/user/n/nminafra/Work/CMSSW/Timing/CMSSW_9_4_0_pre2/src/sendAnalyzerJob.sh AnalyzerRun_$1$2/
cp /afs/cern.ch/user/n/nminafra/Work/CMSSW/Timing/CMSSW_9_4_0_pre2/src/list_of_datGR.py AnalyzerRun_$1$2/
cp /afs/cern.ch/user/n/nminafra/Work/CMSSW/Timing/CMSSW_9_4_0_pre2/src/analyzer_lxbash_run.py AnalyzerRun_$1$2/
cd AnalyzerRun_$1$2
sed -i "s/FIRSTTHREE/${1}/g" list_of_datGR.py
sed -i "s/LASTTHREE/${2}/g" list_of_datGR.py
sed -i "s/FIRSTTHREE/${1}/g" sendAnalyzerJob.sh
sed -i "s/LASTTHREE/${2}/g" sendAnalyzerJob.sh

python list_of_datGR.py

perl -pe 's/FILELIST/`python list_of_datGR.py`/e' ../analyzer_lxbash_run.py > analyzer_lxbash_run.py 
sed -i "s/FIRSTTHREE/${1}/g" analyzer_lxbash_run.py
sed -i "s/LASTTHREE/${2}/g" analyzer_lxbash_run.py

bsub -q 1nd -J Analyzer$1$2 < sendAnalyzerJob.sh -o log_Analyzer$1$2.log -e err_Analyzer$1$2.log
cd ..
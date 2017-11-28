export SCRAM_ARCH=slc6_amd64_gcc530
export BUILD_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/afs/cern.ch/cms
source /afs/cern.ch/cms/cmsset_default.sh
cd /afs/cern.ch/user/n/nminafra/Work/CMSSW/Timing/CMSSW_9_4_0_pre2/src/
eval `scramv1 runtime -sh`
cd -
cmsRun /afs/cern.ch/user/n/nminafra/Work/CMSSW/Timing/CMSSW_9_4_0_pre2/src/AnalyzerRun_FIRSTTHREELASTTHREE/analyzer_lxbash_run.py
mv *.root /afs/cern.ch/user/n/nminafra/Work/public/AnalyzerOut/
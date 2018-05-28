export SCRAM_ARCH=slc6_amd64_gcc530
export BUILD_ARCH=slc5_amd64_gcc462
export VO_CMS_SW_DIR=/afs/cern.ch/cms
source /afs/cern.ch/cms/cmsset_default.sh
cd /afs/cern.ch/work/n/nminafra/CMSSW/UFSD/CMSSW_10_2_0_pre4/src
eval `scramv1 runtime -sh`
cd -
cmsRun /afs/cern.ch/user/n/nminafra/Work/CMSSW/UFSD/CMSSW_10_2_0_pre4/src/DQMRun_FIRSTTHREELASTTHREE/dqm_run.py 2>&1 | grep root
mv DQM_V*.root /afs/cern.ch/user/n/nminafra/Work/public/NinoScanShort/

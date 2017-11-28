## raw data source
import glob
import os
print "*("
filepath='/eos/cms/store/data/Run2017E/ZeroBias/RAW/v1/000/FIRSTTHREE/LASTTHREE/00000/'
onlynames_dat= [os.path.basename(x) for x in glob.glob(filepath+'*.root')]
for i,filename_dat in enumerate(onlynames_dat):
    onlynames_dat[i] ='\'/store/data/Run2017E/ZeroBias/RAW/v1/000/FIRSTTHREE/LASTTHREE/00000/' + filename_dat + '\','
    if i<254000: print onlynames_dat[i]
print ")"
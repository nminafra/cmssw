import os

#Jonathan's Free Running
#path_ = '/store/group/dpg_ctpps/comm_ctpps/HPTDC_calibration/MinidaqRun302086/'

#Global
path_ = '/store/data/Run2017E/ZeroBias/RAW/v1/000/304/447/00000/'
#path_ = '/store/express/Run2017E/ExpressPhysics/FEVT/Express-v1/000/304/292/00000/'
#Minidaq
#path_ = '/store/t0streamer/Minidaq/A/000/303/969/'
file_save = open("./NickAnlzr/MyAnalyzer/python/AutoGenerate_cff.py",'w') 

asps = []
print>>file_save, "readFiles=["
for root, dirs, files in os.walk(r'/eos/cms'+path_):
    for file in files:
        if file.endswith('.root'):
            asps.append(file)
            print>>file_save, "'"+path_ + file + "',"
print>>file_save, "]"
file_save.close()

import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root',',' with the source file you want to use
    fileNames = cms.untracked.vstring(
      *(
        
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/0627D30A-2CA3-E711-8C9B-02163E01A2E6.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/0E3119DB-26A3-E711-8360-02163E01A6AB.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/0E3D0117-29A3-E711-B482-02163E019D8F.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/10527545-2EA3-E711-B309-02163E019B8B.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/1882CBE8-29A3-E711-8A96-02163E01218E.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/1E01BC4A-52A3-E711-8FB2-02163E01A23D.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/1EA7648B-2CA3-E711-908E-02163E01A56E.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/20A01A47-2EA3-E711-9164-02163E01450D.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/2293B910-44A3-E711-BCA6-02163E019D2D.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/32F510C8-2CA3-E711-A102-02163E01418B.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/36CCC1EE-33A3-E711-9CE6-02163E0138FE.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/38B9403A-2FA3-E711-92B2-02163E019E47.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/3A0A3908-2AA3-E711-A171-02163E019E6E.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/4071A4B9-3BA3-E711-979A-02163E01A736.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/42D5ACE1-3BA3-E711-9D45-02163E019C7E.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/4A9E7373-30A3-E711-ACA4-02163E019D2C.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/4AD033C4-2EA3-E711-9110-02163E0146F2.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/4E436866-31A3-E711-AE4A-02163E0137F1.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/4EA7426D-2DA3-E711-99E2-02163E012ABD.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/54743BA4-2AA3-E711-A805-02163E01429D.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/5670C1C8-40A3-E711-8CE2-02163E01A1C1.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/567B6737-30A3-E711-8F79-02163E019D8F.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/5A06C0EA-49A3-E711-AE1F-02163E01234F.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/5CAD057D-37A3-E711-A979-02163E019C4D.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/60D9715D-2EA3-E711-B58B-02163E011E0C.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/60DDD34B-43A3-E711-86EC-02163E012B3A.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/6AF60F71-28A3-E711-BD87-02163E01A47E.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/6AFD1AC7-29A3-E711-8B4A-02163E01464E.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/6C3A9065-2CA3-E711-8E1F-02163E01194C.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/6C68F1D1-26A3-E711-A539-02163E011BB9.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/6CB2961D-32A3-E711-9948-02163E01A26E.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/6EE07093-2AA3-E711-917D-02163E01440F.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/722293A0-32A3-E711-A636-02163E01246F.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/74AABFEB-2BA3-E711-87A9-02163E0144FA.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/7A5D915D-25A3-E711-B085-02163E0146D7.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/8014EFE2-2AA3-E711-8633-02163E01A5F3.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/804AD3F0-39A3-E711-95CA-02163E014627.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/8239B2ED-26A3-E711-A808-02163E014369.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/8A8414E7-2CA3-E711-9CC8-02163E014492.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/8EC601BE-2CA3-E711-8CF4-02163E01464F.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/9091777D-2AA3-E711-BB22-02163E0146D7.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/927231E2-3CA3-E711-B611-02163E01A2D4.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/94AB1FA6-2CA3-E711-A59E-02163E011DBE.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/94FC0177-3CA3-E711-8589-02163E01A468.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/9A1FB695-22A3-E711-B85D-02163E014200.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/9C4D85B8-2DA3-E711-B00D-02163E01A29F.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/9C7C344A-40A3-E711-8D8F-02163E014661.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/A04E8CCA-23A3-E711-83F4-02163E0123DD.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/A25B6C4A-46A3-E711-98A7-02163E014120.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/A2F94D82-34A3-E711-8126-02163E0137F1.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/A43B1DEC-2AA3-E711-8865-02163E0136CC.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/A4914E1B-42A3-E711-8728-02163E013502.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/A4999C0B-42A3-E711-9332-02163E019D8F.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/A4DA5E22-2DA3-E711-A253-02163E0137F1.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/A651EB30-28A3-E711-BF3B-02163E01A2C7.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/A8341CF8-27A3-E711-85D7-02163E019E23.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/AE3F1732-37A3-E711-9223-02163E019E03.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/B0A4263D-29A3-E711-A546-02163E01A208.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/B47697B6-2CA3-E711-B5C5-02163E01A2F5.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/C0AE7631-3EA3-E711-B99E-02163E013543.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/C2768A87-19A3-E711-A042-02163E0144F6.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/C45A08D0-35A3-E711-B22A-02163E01A6C0.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/C65794BE-3EA3-E711-8DAC-02163E01A36C.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/C82A638A-38A3-E711-B94B-02163E0146E6.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/CA83204F-2EA3-E711-A198-02163E01A2D2.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/D0C8F5D2-2FA3-E711-BA5C-02163E01446D.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/D439C2B0-2EA3-E711-9ADD-02163E0139B8.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/DCFAEDC2-2CA3-E711-9BCF-02163E01A4B1.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/E483B867-3BA3-E711-A3BE-02163E01A6C0.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/E662F3E0-32A3-E711-BFE4-02163E0146F2.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/E86CA36F-2EA3-E711-908D-02163E011882.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/EAD580F0-28A3-E711-AE1F-02163E01194C.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/ECEFCFD1-31A3-E711-83CF-02163E01446D.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/F400DE47-3FA3-E711-B70A-02163E01349B.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/F60DF947-29A3-E711-92E3-02163E01A64A.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/F657A057-28A3-E711-AEED-02163E01A36C.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/F677554A-3DA3-E711-8D83-02163E014661.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/F80A077C-6AA3-E711-83F1-02163E0139AF.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/F8713A53-3AA3-E711-82EC-02163E01A384.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/FC83EA78-3AA3-E711-9728-02163E0133F6.root',
        '/store/data/Run2017E/ZeroBias/AOD/PromptReco-v1/000/303/824/00000/FCBE8664-31A3-E711-882E-02163E01349B.root'
     
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/002BA898-8A89-E711-BA50-02163E01A512.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/0650418D-5489-E711-AFFC-02163E01A4C3.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/0842BDF3-7389-E711-A470-02163E012B9B.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/0EB395B7-7489-E711-BA39-02163E019BC5.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/1054EDEA-BE89-E711-B47D-02163E014410.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/12EF7722-9789-E711-A10B-02163E0146E6.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/1C9C5BA7-8B89-E711-8155-02163E012528.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/1ED898A4-5489-E711-9BB6-02163E011C5B.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/2492A22B-7789-E711-BB8C-02163E01A1E4.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/2690F28C-7A89-E711-BCB3-02163E01A5EB.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/2A7C62F7-8289-E711-A8F7-02163E01200E.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/34AA8791-5489-E711-BD32-02163E01A329.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/34D90D77-8789-E711-9420-02163E01A377.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/3A7A5DBE-9489-E711-BAB4-02163E012234.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/3CAE0E90-5489-E711-AF02-02163E019BAD.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/3E5B525B-5589-E711-9D44-02163E011A30.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/3EB0A8BF-5489-E711-9A9B-02163E01345E.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/40DB0BA3-5489-E711-9EDF-02163E014325.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/4869D49E-5489-E711-991A-02163E011E3E.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/580F2761-8989-E711-ADF5-02163E01A7A4.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/58228D8B-5489-E711-8774-02163E01A233.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/5CBF3898-5589-E711-84A3-02163E0143DA.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/681A05C1-8089-E711-8E21-02163E011AA4.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/6A5D0E99-8C89-E711-A907-02163E01241A.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/6AC1AEF8-8789-E711-9971-02163E01A3D1.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/76F1599A-5489-E711-AC3E-02163E01413D.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/7AD71B83-8D89-E711-8966-02163E0142CC.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/845DDAC4-8989-E711-8BF8-02163E0144B9.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/84EDD80B-7389-E711-84D2-02163E01344B.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/88F24956-8B89-E711-8422-02163E01A5FD.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/8A9A368D-5489-E711-B3BD-02163E019E88.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/947E155C-8689-E711-AB50-02163E019B45.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/988C9CF4-7289-E711-8B57-02163E01A41D.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/9E5DA3A8-5489-E711-9194-02163E01A6FF.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/9EECDCA0-7089-E711-9279-02163E01448A.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/A0E4B8A2-9889-E711-8F2B-02163E0119B7.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/A6009093-5489-E711-B448-02163E019D46.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/A813A61A-7589-E711-9B3F-02163E014291.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/AE8AA7EF-6289-E711-956C-02163E01475C.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/C08AB90C-8189-E711-AE33-02163E01A49A.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/CE45DD3F-8989-E711-9095-02163E019C04.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/D2124422-6F89-E711-9C47-02163E019C1A.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/D290F7A1-5489-E711-AE0A-02163E011F68.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/DA59EC97-7089-E711-AB76-02163E01A464.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/DE15ED30-8F89-E711-879C-02163E013621.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/E25C94C9-7789-E711-99B9-02163E019DD0.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/EC1DB989-5489-E711-9CDE-02163E019E69.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/ECBC0340-7889-E711-B37D-02163E014581.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/F6CC8987-5489-E711-B14B-02163E0142D0.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/FAC4B03B-5589-E711-8A77-02163E011A13.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/FAD4948D-6E89-E711-9C3D-02163E01A269.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v3/000/301/627/00000/FADEE23B-7C89-E711-9A50-02163E014663.root',
##    
     )
    ),
    #inputCommands = cms.untracked.vstring(
       #'drop CTPPSPixelCluseredmDetSetVector_ctppsPixelClusters__RECO'
    #)
)
    
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_hlt_relval', '')

# raw-to-digi conversion
process.load("EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff")

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")

#process.ctppsDiamondRecHits.planeInversion = cms.int32(1)               # Put >0 for Runs < 300670     #FIXME
#process.ctppsDiamondRecHits.coarseCorrection = cms.int32(1)               # Put >0 for Runs < ??       #FIXME

# rechits production
process.load('Geometry.VeryForwardGeometry.geometryRP_cfi')
process.load('RecoCTPPS.TotemRPLocal.ctppsDiamondRecHits_cfi')

# local tracks fitter
process.load('RecoCTPPS.TotemRPLocal.ctppsDiamondLocalTracks_cfi')

process.ctppsDiamondLocalTracks.trackingAlgorithmParams.threshold = cms.double(1.9)
process.ctppsDiamondLocalTracks.trackingAlgorithmParams.sigma = cms.double(0.1)
process.ctppsDiamondLocalTracks.trackingAlgorithmParams.resolution = cms.double(0.01) # in mm
#process.ctppsDiamondLocalTracks.trackingAlgorithmParams.pixelEfficiencyFunction = cms.string("(TMath::Erf((x-[0]+0.5*[1])/([2]/4)+2)+1)*TMath::Erfc((x-[0]-0.5*[1])/([2]/4)-2)/4")

# pixel
process.load('RecoCTPPS.PixelLocal.ctppsPixelLocalTracks_cfi')

#process.load("MyAnalyzer.CfiFile_cfi.py")
process.demo = cms.EDAnalyzer('MyAnalyzer',
 tagDigi = cms.InputTag("ctppsDiamondRawToDigi", "TimingDiamond"),
 tagDiamondRecHits = cms.InputTag("ctppsDiamondRecHits"),
 tagDiamondLocalTracks = cms.InputTag("ctppsDiamondLocalTracks"),
 tagLocalTrack = cms.InputTag("totemRPLocalTrackFitter"),
 ctppsPixelLocalTracks = cms.InputTag("ctppsPixelLocalTracks"),
 selectedOOTIndex = cms.int32(0)                                                                        #FIXME
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('303824_th19_sg01_2pl_toll05mm_long.root')                                                   #FIXME
)

process.p = cms.Path(
    #process.ctppsRawToDigi *
    process.recoCTPPS *
    #process.ctppsDiamondRawToDigi *
    process.ctppsDiamondRecHits *
    process.ctppsDiamondLocalTracks *
    process.ctppsPixelLocalTracks *
    process.demo
    )


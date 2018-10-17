import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
	secondaryFileNames = cms.untracked.vstring(),
	# fileNames = cms.untracked.vstring(options.files),
	fileNames = cms.untracked.vstring(
    *(
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/FEE64C09-FC97-E811-804D-FA163E1CED1A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/FE491CAC-1C98-E811-B29B-FA163E3FD732.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/FCE1B388-1798-E811-A722-FA163ECFF8BE.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/FCBE5BC9-0B98-E811-A4F7-FA163E5B858A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/FC8801B7-1298-E811-AA6A-FA163E47622D.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/FC2ED74E-1398-E811-B442-02163E010C02.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/FC192E21-1398-E811-92CD-02163E010C9D.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/FAED9123-1798-E811-AD44-FA163EB85A87.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/FA03C979-1498-E811-B071-A4BF01277792.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/F6C8FAA5-1D98-E811-A38E-FA163E2C305A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/F634C2BA-0C98-E811-9CDD-FA163EA20D66.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/F4F2710B-1498-E811-8EA7-FA163E278EA6.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/F2BE7E9E-1398-E811-922D-FA163E706B48.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/EED458B8-1C98-E811-8F69-FA163E17F3DA.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/EED0B2F2-1898-E811-8005-FA163E0EDCAA.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E8F151FC-1898-E811-940C-FA163ECCD3B2.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E8E67354-1D98-E811-AD7E-FA163EDE923A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E8607CFE-1898-E811-BC52-FA163E23298B.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E68EC88F-1E98-E811-9B44-FA163EA89334.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E4E67EF5-1598-E811-AD5E-FA163E4C1BDF.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E4D6590D-1598-E811-962E-FA163EACFE1F.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E4CA420F-1398-E811-9E73-02163E019F3B.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E48FA5D6-1498-E811-9811-02163E010BCA.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E485B4FE-1498-E811-BC29-FA163EB5453C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E42F8C27-1A98-E811-AD91-FA163EB5453C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E2B339F2-1598-E811-A072-FA163E33FA23.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E2AFD605-1E98-E811-A2CF-FA163E4E87E3.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E271DDFF-1498-E811-97CE-FA163E10365B.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E262CA32-1D98-E811-AED4-FA163E0E994B.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E2281DEC-0E98-E811-9B06-FA163EA2F461.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E0F5E46A-0C98-E811-9CEA-FA163EE4DD06.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E05A177F-0B98-E811-B0E2-FA163EF68D36.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/E003CBB0-1398-E811-A41E-FA163EC92774.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/DE28CC8B-1498-E811-8C60-FA163EB9CE9F.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/DCBB28B1-0B98-E811-B6FB-02163E017FE5.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/DACC6CE2-1D98-E811-AEA4-FA163E22812C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/D8FEC4EA-0F98-E811-BA5A-02163E019FA2.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/D8DDC582-1A98-E811-9543-FA163E120D15.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/D867DECD-1398-E811-8238-FA163E6A7294.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/D6970D99-1398-E811-AB84-FA163EDE5412.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/D4D25469-1398-E811-BE4F-FA163E8F367B.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/D44392F0-0A98-E811-B6FD-FA163E32FF7C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/D439D410-FE97-E811-8244-FA163E03A9B7.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/D2F33068-1998-E811-A389-FA163E2B0B90.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/D2F1A71D-2698-E811-9D62-FA163E984001.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/D2158819-1798-E811-83A7-FA163E607C1F.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/D20AD6CA-1D98-E811-80A1-FA163E2E4F0A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/CE6E72A7-0D98-E811-968B-FA163E6F8F72.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/CE30940D-1498-E811-A7E1-FA163E13EBA8.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/CE2E17BB-3F98-E811-8498-FA163E3DDA6C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/CE23458A-0C98-E811-B929-FA163E74AC8B.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/CCF8BC58-1B98-E811-A81F-FA163E104EBC.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/CCD16447-1098-E811-929C-FA163EBF60A5.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/CC7989C8-0D98-E811-AF5B-02163E016007.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/CAD1FDBC-1298-E811-80C4-FA163E1E616C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/CA395E39-1498-E811-9C3A-FA163EEF5766.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/C8E6481F-1698-E811-BCBE-02163E01A097.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/C6D2BFEE-1298-E811-99E7-FA163EC2F360.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/C6B932AA-1B98-E811-846B-02163E017758.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/C6301386-1098-E811-8D5A-FA163E167E2A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/C4A65DBC-1198-E811-871F-FA163E4A0FB4.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/C45CB49B-FD97-E811-A35C-FA163E2BCBE1.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/C27808F7-1098-E811-9EB2-FA163E6311DA.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/C081D42E-5198-E811-9AED-FA163EF3E48B.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/BEF4E535-1398-E811-848B-FA163E68A5F1.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/BEEF79EF-1098-E811-B3D4-FA163EF87070.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/BED131F5-0798-E811-B993-FA163E8C7780.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/B87FF383-1498-E811-B9FB-FA163E11ED71.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/B850C2B2-1398-E811-AEAC-FA163E0581F0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/B6931413-1698-E811-A659-FA163EFC4E1B.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/B68FBD7E-1498-E811-8DD5-FA163E2FA06D.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/B4EC092A-1298-E811-B60A-FA163E503204.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/B4B08FAE-1198-E811-AEC8-FA163EB2E120.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/B4069CB6-0B98-E811-9CAC-FA163E7C71B8.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/B0EEB902-0E98-E811-A411-02163E017F05.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/AEB7477F-1498-E811-8753-02163E010DFC.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/AEAF16E1-1D98-E811-B29A-FA163E659E42.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/ACA40F94-1598-E811-A83F-FA163E249C0F.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/AC0E1D89-0D98-E811-9486-FA163EAD3FD0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/AA4F530C-1A98-E811-9DE9-FA163E4E87E3.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/A897696F-3C98-E811-BA8C-02163E0152B6.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/A87AEF15-1998-E811-AD15-FA163E5B85EC.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/A821DBCF-0D98-E811-B2CD-FA163E7E98D5.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/A8064FFD-1398-E811-B18C-FA163EE655A1.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/A6CD4A8E-1298-E811-8319-FA163E04457F.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/A6BD92EC-1298-E811-A82F-FA163E617795.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/A64CC1A5-0A98-E811-AA5A-02163E0153BD.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/A4B1D617-0C98-E811-A6E3-FA163E476E6A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/A4279D66-1298-E811-ACA4-FA163E31B17C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/A24DA80B-1198-E811-88C1-FA163EC94D39.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/A23F987F-1A98-E811-90C9-FA163E1AF201.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/A2337734-1A98-E811-8049-FA163E8F473B.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/A0BCDF71-1598-E811-B644-FA163E0B896C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/9EDF3003-0D98-E811-BD6E-FA163E11ED71.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/9E6F38A4-1598-E811-B7F4-02163E019FA2.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/9E1E85B0-1298-E811-8393-02163E010D44.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/9ACC7503-1698-E811-8FA2-FA163E477BED.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/9A949A5B-1798-E811-B7E9-02163E01A140.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/9A769404-1598-E811-9214-FA163E7B2F96.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/96F2F779-1498-E811-A424-FA163EED3BED.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/96846F5B-1098-E811-966B-FA163E01FE18.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/96421EFD-F797-E811-A686-FA163E1E6FE0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/94819451-1898-E811-B7A9-FA163EA1B107.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/9424D6C1-1398-E811-B74B-02163E0153BD.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/92EF51C1-0B98-E811-8008-02163E013029.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/927E8787-0D98-E811-9155-FA163E0EDCAA.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/90C359F8-1898-E811-8457-FA163EBF9BFB.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/90C1FEB9-1D98-E811-B434-FA163ED486E3.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/90207747-2C98-E811-892C-FA163E26831C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/8EF00EEC-1798-E811-BA53-FA163E9673D0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/8E6DDC3F-1598-E811-85C1-FA163EC3F894.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/8E150CCD-0A98-E811-B950-FA163EDC261A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/8CDFEECB-FC97-E811-976F-FA163EBCA265.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/8C7971DB-FD97-E811-BBB7-FA163E3EEE2D.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/8ACE0470-1498-E811-982A-FA163E607C1F.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/8AC7E376-0D98-E811-A60D-FA163E535736.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/8ABCFEA0-1498-E811-97B2-02163E019FEF.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/8AA7F584-1298-E811-9128-FA163EE79415.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/88B7114E-1098-E811-A51B-FA163E490F23.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/88460B84-1098-E811-84B0-02163E01A140.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/86DA7E47-1098-E811-8D2E-FA163E0581F0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/869612A9-1398-E811-A2B3-FA163EE31D2A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/864D47BA-1698-E811-8AAB-FA163E646768.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/8442C62B-1598-E811-A3ED-FA163ECFF8BE.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/8439AF7A-1398-E811-A8B5-FA163E2D3D79.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/842123CC-1A98-E811-A81A-FA163E83935C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/82233500-1698-E811-8BEA-A4BF0114C8F0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/8205F976-1798-E811-B845-FA163E335778.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/80A7EAAB-0F98-E811-A972-FA163E51B1D6.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/80771BFB-1498-E811-9347-FA163E4BB802.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/7EA25E4E-0F98-E811-93CF-FA163EC4736B.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/7E5B1E02-FC97-E811-A7C5-FA163EF74259.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/7CCC1396-0D98-E811-A8D9-FA163E6218E7.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/7ABE1AE4-0E98-E811-880C-FA163E60ED2F.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/7A8F7C1A-1998-E811-AC05-FA163E617795.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/7A4E69A7-1298-E811-B531-FA163E3C94E9.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/789F0F03-FE97-E811-9F68-FA163E692103.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/7886CBAA-1C98-E811-A80A-FA163E1B411E.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/76FF7112-1598-E811-93E7-FA163E476E6A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/749848E2-0E98-E811-B06A-FA163EDE5CFD.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/72B068BC-1298-E811-91F7-FA163EECA815.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/726925B0-1C98-E811-B020-FA163EC7CC24.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/7262D1E4-0698-E811-94DE-02163E017FA3.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/701C2338-1398-E811-839A-FA163E8D01A3.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6EC1EF91-1298-E811-86B0-02163E01A11A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6E9B2031-FF97-E811-A9CE-FA163E35E31D.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6CFE6325-1698-E811-85BE-FA163E7CF385.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6CE2300B-1B98-E811-A012-FA163EA6EFC1.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6CC99427-1398-E811-A941-02163E017FE5.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6C239012-1498-E811-9857-FA163E1DF5C0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6AFB1CA7-1C98-E811-BBBB-A4BF01125A78.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6AA955D0-FD97-E811-A508-FA163E69B278.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6A7846D4-FD97-E811-804B-FA163E6C3064.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6834178A-1798-E811-8868-FA163E2D3D79.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/664896F7-1D98-E811-9506-FA163EA4D9B7.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6629A1E6-0B98-E811-9A61-FA163E60FDAB.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/661CE18F-0C98-E811-B345-FA163E908B47.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/640193F8-1D98-E811-A4C2-FA163E942E01.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/62C0F104-1498-E811-98EF-FA163EDE5CFD.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6273EE59-0F98-E811-AF35-02163E019F3B.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6237E6E9-1898-E811-AF2E-FA163E278EA6.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/6213469D-F797-E811-81F0-02163E017672.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/602A5F3B-0D98-E811-BD36-02163E00BF4D.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/5EF71843-1798-E811-9E9A-FA163E0D54BF.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/5ECD1C48-1098-E811-A139-FA163EA21B5C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/5E3BC978-1A98-E811-8EB3-FA163E5441AB.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/5C8DBE46-1B98-E811-A31E-FA163E878B53.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/5C4AD3D1-0E98-E811-9826-02163E01A142.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/5A659CB7-1E98-E811-9CF6-FA163E9E3073.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/5A3FE7C8-1F98-E811-8785-FA163EEBD149.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/58F7DC24-1398-E811-A733-FA163E884269.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/5432A0EF-1698-E811-8F05-FA163E4842F0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/52F572F6-1798-E811-ABD7-FA163EB9CE9F.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/52E09033-2798-E811-A501-FA163EDDE351.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/52C0B049-2598-E811-BF74-FA163E2F466D.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/520D7388-1398-E811-A1CB-FA163EC92B39.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/50328FCC-1298-E811-8F11-FA163EEE65D6.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/4C6593D2-2998-E811-9FA6-FA163ED01E21.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/4ADE125F-1A98-E811-955D-FA163E1E24E8.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/4AA27830-1598-E811-AA28-FA163EBF9BFB.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/4A899F72-1098-E811-B157-FA163E42FED6.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/4A29456F-0F98-E811-BA3F-FA163E22812C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/4A02113A-1D98-E811-9CC7-FA163E0262D6.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/48612C8F-2298-E811-BC3C-FA163EE4DD06.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/482A092C-1298-E811-9929-FA163EB005BF.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/4665F7B2-0E98-E811-8FA8-02163E010DBD.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/463960EF-1298-E811-9A01-02163E01A029.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/44F051EB-1A98-E811-861E-FA163E81F4D7.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/44300D0D-2198-E811-AAC7-FA163EAF5FB5.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/42F9B346-1098-E811-BDF5-FA163ECF720C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/42CAD63D-1698-E811-B136-02163E019F5A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/42AA3AA5-1A98-E811-B80A-02163E01A142.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/4295D43B-1398-E811-AE9A-FA163E4D6DEF.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/40DD97CF-FD97-E811-90AD-FA163E0D6C88.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/4089DF16-1398-E811-8B0F-FA163EBDCDFA.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/406F7B4D-0F98-E811-A25B-FA163EA81608.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/405F2CE9-0998-E811-8B3E-02163E01A043.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/40030022-1398-E811-8470-FA163E7BAAF4.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/3E316372-0B98-E811-AEF0-FA163E4842F0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/3CD9CF17-0D98-E811-8BAA-FA163E3796B7.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/3C5A0A1B-2698-E811-B6DB-FA163E6BF2DF.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/3C25F574-1B98-E811-A5B6-FA163E8CC774.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/3AD919BB-1C98-E811-9793-A4BF0114CCA8.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/3A7455EB-1098-E811-A16B-02163E01A01E.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/38F9EF32-1A98-E811-BC00-02163E010C02.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/38E30276-1B98-E811-8A30-FA163E6C0084.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/3867C8D1-FD97-E811-AE54-FA163E2D2D27.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/36ABEF0B-1498-E811-A069-FA163ED59DED.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/364E3735-1498-E811-B939-FA163E9240FC.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/361D520E-1598-E811-A0C7-FA163E959A66.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/34D886FC-1598-E811-A479-FA163E78505A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/3028F043-0F98-E811-B210-FA163EC21F46.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/2E14D7A0-1398-E811-B7D0-02163E01A155.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/2C613620-0898-E811-BD22-FA163EC375AA.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/2C01C3DA-FD97-E811-A834-FA163EA3DAC4.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/28E0E383-1798-E811-8D0B-FA163EB7B53C.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/28996427-1798-E811-B4DA-02163E01A029.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/28844551-1798-E811-B061-FA163EC6F3D0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/26CEFB7B-0C98-E811-80B7-FA163E70CD9E.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/2682169A-1098-E811-8D02-FA163EECC2C3.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/24F45077-1098-E811-ADAF-FA163E517EF1.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/22A355CD-1C98-E811-BB5F-02163E019F84.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/221E7D94-1298-E811-8738-FA163ECF8401.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/208B1B95-1298-E811-AD3E-FA163E1B9BC3.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/208AF83C-0C98-E811-9AAA-FA163EC94D39.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/2069D98B-1698-E811-A245-FA163E016D86.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/20581E5E-2E98-E811-ABB4-FA163E47686E.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/1E9A2FAD-1098-E811-810B-FA163E47C8BD.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/1CF04D74-1398-E811-A228-FA163E9673D0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/1CEA34FC-1398-E811-A495-FA163EE4FBEB.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/1C66525E-0998-E811-AC5E-FA163ED0BA63.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/1AE8F871-0698-E811-BE6E-FA163EB2FB89.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/1AA9AE0C-1398-E811-8A95-FA163E42FA7D.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/1AA9370B-1D98-E811-9F1C-FA163EBF8C10.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/18E16FC7-1098-E811-98CD-FA163E104BF4.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/183B74BD-2F98-E811-8381-FA163E7C91C0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/16C61139-0B98-E811-9AD5-02163E019FA2.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/16C05BA9-0B98-E811-A9AD-FA163E7537E4.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/167311FA-0198-E811-A032-FA163E8C7780.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/1625B6E2-0D98-E811-A0F8-FA163EA4AD14.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/161180C5-1A98-E811-9F68-FA163E78505A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/14FE2DC5-1D98-E811-9C18-FA163EB2B58A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/1254E0C9-1798-E811-933E-FA163E10365B.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/124843AA-1C98-E811-A0F7-02163E01A01E.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/120B9090-1798-E811-9B08-FA163E13EBA8.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/10FE4B98-1298-E811-80BC-FA163EAD3FD0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/10D84E2C-1A98-E811-88C4-FA163E721D55.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/10B2AB77-1598-E811-9CE1-FA163ED0947A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/10ACBDA2-2E98-E811-AF02-02163E01767E.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/10A3DCE6-1A98-E811-BAE6-FA163E4041C8.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/10542F43-1A98-E811-A6A1-FA163EEFAAC2.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0ECDF131-1098-E811-A2B5-FA163E2890BA.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0E153D36-0F98-E811-8E64-FA163EE4FBEB.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0E127A58-0D98-E811-808A-02163E013C13.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0CDD7B3F-1A98-E811-B143-FA163EA30F9A.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0C5784B5-1C98-E811-9190-FA163EBDEC47.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0AD85F89-1698-E811-9C10-FA163EE70CD0.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0AAE4B26-1798-E811-8D11-FA163E68A5F1.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0AA33BC5-2498-E811-8A77-FA163EF9E3A4.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0A77D643-1498-E811-BACE-FA163EC94D39.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0A74A7F3-1D98-E811-910A-02163E01A0FF.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0A154C10-1698-E811-8AE7-FA163E48F925.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0873AD01-0D98-E811-926F-FA163ECDE73F.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0848B186-1398-E811-859D-FA163E4C1BDF.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/04C89B6F-1498-E811-A13F-FA163E335778.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/0475A640-1898-E811-8DF3-02163E010E8B.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/02D90789-1B98-E811-A212-FA163EA222CE.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/00E96EC6-1398-E811-BEE3-FA163E4CBD04.root',
'/store/data/Run2018D/EGamma/AOD/PromptReco-v2/000/320/688/00000/006A98E9-0998-E811-9059-02163E017F86.root',
    )),
	skipEvents = cms.untracked.uint32(0)
)

process.demo = cms.EDAnalyzer('VertexDistr',
pfTag = cms.InputTag('particleFlow'),
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('output_long.root')
)

process.p = cms.Path(process.demo)

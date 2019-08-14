import FWCore.ParameterSet.Config as cms

#define the default IDs to produce in VID
_defaultPhoIDModules = [ 
                         'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V1_TrueVtx_cff',
                         'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_OOT_V1_cff',
                       ]

def _setupOOTPhotonPostRECOSequenceMiniAOD(process):
    
    phoSrc = cms.InputTag('slimmedOOTPhotons',processName=cms.InputTag.skipCurrentProcess())
    phoCalibSrc = cms.InputTag('slimmedOOTPhotons',processName=cms.InputTag.skipCurrentProcess())

    process.load('RecoEgamma.EgammaTools.calibratedEgammas_cff')
    process.calibratedPatPhotons.src = phoCalibSrc
    
    energyCorrectionFile = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc"
    process.calibratedPatPhotons.correctionFile = energyCorrectionFile

    process.calibratedPatPhotons.produceCalibratedObjs = False 

    process.egmPhotonIDs.physicsObjectSrc = phoSrc
    process.photonMVAValueMapProducer.src = phoSrc
    process.photonIDValueMapProducer.srcMiniAOD = phoSrc
    process.egmPhotonIsolation.srcToIsolate = phoSrc

    from RecoEgamma.EgammaTools.egammaObjectModificationsInMiniAOD_cff import egamma_modifications
    from RecoEgamma.EgammaTools.photonObjectModifications_tools import makeVIDBitsModifier,makeVIDinPATIDsModifier,makeEnergyScaleAndSmearingSysModifier  
    
    egamma_modifications.append(makeVIDBitsModifier(process,"egmPhotonIDs"))
    egamma_modifications.append(makeVIDinPATIDsModifier(process,"egmPhotonIDs"))
    egamma_modifications.append(makeEnergyScaleAndSmearingSysModifier("calibratedPatPhotons"))

    for pset in egamma_modifications:
        pset.overrideExistingValues = cms.bool(True)
        if hasattr(pset,"photon_config"): pset.photon_config.photonSrc = phoSrc

    process.slimmedOOTPhotons = cms.EDProducer("ModifiedPhotonProducer",
                                               src=phoSrc,
                                               modifierConfig = cms.PSet( modifications = egamma_modifications )
                                               )

    process.ootPhotonScaleSmearTask = cms.Task( process.calibratedPatPhotons, process.slimmedOOTPhotons )

def setupOOTPhotonPostRecoSeq(process):

    from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDPhotonIdProducer,setupAllVIDIdsInModule,DataFormat,setupVIDPhotonSelection

    switchOnVIDPhotonIdProducer(process,DataFormat.MiniAOD)

    for idmod in _defaultPhoIDModules:
        setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

    _setupOOTPhotonPostRECOSequenceMiniAOD(process)
    
    process.ootPhotonScaleSmearSeq = cms.Sequence( process.ootPhotonScaleSmearTask)
    process.ootPhotonPostRecoSeq   = cms.Sequence( process.ootPhotonScaleSmearSeq * process.egmPhotonIDSequence )

    return process

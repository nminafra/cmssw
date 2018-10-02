import FWCore.ParameterSet.Config as cms

#define the default IDs to produce in VID
_defaultPhoIDModules = [ 
                         'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_OOT_V1_cff',
                       ]

def _setupGEDPhotonPostRECOSequenceMiniAOD(process):
    
    phoSrc = cms.InputTag('slimmedPhotons',processName=cms.InputTag.skipCurrentProcess())

    process.egmPhotonIDs.physicsObjectSrc = phoSrc
    process.photonMVAValueMapProducer.src = phoSrc
    process.photonIDValueMapProducer.srcMiniAOD = phoSrc
    process.egmPhotonIsolation.srcToIsolate = phoSrc

    from RecoEgamma.EgammaTools.egammaObjectModificationsInMiniAOD_cff import egamma_modifications
    from RecoEgamma.EgammaTools.photonObjectModifications_tools import makeVIDBitsModifier,makeVIDinPATIDsModifier
    
    egamma_modifications.append(makeVIDBitsModifier(process,"egmPhotonIDs"))
    egamma_modifications.append(makeVIDinPATIDsModifier(process,"egmPhotonIDs"))

    for pset in egamma_modifications:
        #        pset.overrideExistingValues = cms.bool(True)
        if hasattr(pset,"photon_config"): pset.photon_config.photonSrc = phoSrc

    process.slimmedPhotons = cms.EDProducer("ModifiedPhotonProducer",
                                            src=phoSrc,
                                            modifierConfig = cms.PSet( modifications = egamma_modifications )
                                            )

    process.gedPhotonOOTVIDTask = cms.Task( process.slimmedPhotons )

def setupGEDPhotonPostRecoSeq(process):

    from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDPhotonIdProducer,setupAllVIDIdsInModule,DataFormat,setupVIDPhotonSelection

    switchOnVIDPhotonIdProducer(process,DataFormat.MiniAOD)

    for idmod in _defaultPhoIDModules:
        setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

    _setupGEDPhotonPostRECOSequenceMiniAOD(process)
    
    process.gedPhotonOOTVIDSeq   = cms.Sequence( process.gedPhotonOOTVIDTask )
    process.gedPhotonPostRecoSeq = cms.Sequence( process.gedPhotonOOTVIDSeq * process.egmPhotonIDSequence )

    return process

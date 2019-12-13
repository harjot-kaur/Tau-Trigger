import FWCore.ParameterSet.Config as cms

print "Running on MC"


HLTLIST = cms.VPSet(
    #cms.PSet (
    #    HLT = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v"),
    #    path1 = cms.vstring ("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu17LooseIsoPFTau20"),
    #    path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterIsoMu17LooseIsoPFTau20"),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(15)
    #)
    cms.PSet (
        HLT = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09", "hltOverlapFilterSingleIsoMu17LooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterSingleIsoMu17LooseIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    )
)



# filter HLT paths for T&P
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
hltFilter = hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    #TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    #TriggerResultsTag = cms.InputTag("TriggerResults","","RECO"),
    #TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    #TriggerResultsTag = cms.InputTag("TriggerResults","","HLT2"),
    HLTPaths = ['HLT_IsoMu18_v5'],
    #HLTPaths = ['HLT_Mu7p5_L2Mu2_Jpsi_v4'],
    #HLTPaths = ['HLT_IsoMu18_v3'],
    #HLTPaths = ['HLT_IsoMu18_v3'],
    andOr = cms.bool(True), # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(True) #if True: throws exception if a trigger path is invalid
)

## good taus - apply analysis selection
goodTaus = cms.EDFilter("PFTauSelector",
        src = cms.InputTag("hpsPFTauProducer"),
        cut = cms.string(
                'pt > 18 && abs(eta) < 2.5 ' #kinematics
                '&& abs(charge) > 0 && abs(charge) < 2 '
        ),
        discriminators = cms.VPSet(
        cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),selectionCut=cms.double(0.5)),
        cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits"),selectionCut=cms.double(2.5)),
        cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection3"),selectionCut=cms.double(0.5)),
        cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByMVA6VLooseElectronRejection"),selectionCut=cms.double(0.5)),
        ),
        filter = cms.bool(True)
)

from PhysicsTools.JetMCAlgos.TauGenJets_cfi import tauGenJets

trueHadronicTaus = cms.EDFilter(
    "TauGenJetDecayModeSelector",
    src = cms.InputTag("tauGenJets"),
    select = cms.vstring(
        'oneProng0Pi0', 'oneProng1Pi0', 'oneProng2Pi0', 'oneProngOther',
        'threeProng0Pi0', 'threeProng1Pi0', 'threeProngOther', 'rare'),
    filter = cms.bool(False)  
)

genMatchedTaus = cms.EDFilter("genMatchTauFilter_AOD",
        taus = cms.InputTag("goodTaus"),
        genJets = cms.InputTag("trueHadronicTaus")
    )

# Ntuplizer.taus = cms.InputTag("genMatchedTaus")
Ntuplizer_noTagAndProbe = cms.EDAnalyzer("Ntuplizer_noTagAndProbe_AOD",
    treeName = cms.string("TagAndProbe"),
    taus  = cms.InputTag("genMatchedTaus"),
    triggerSet = cms.InputTag("selectedPatTrigger"),
    triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
    #triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
    #triggerResultsLabel = cms.InputTag("TriggerResults", "", "RECO"),
    #triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT2"),
    #triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
    #L1Tau = cms.InputTag("caloStage2Digis", "Tau", "HLT"),
    #L1Tau = cms.InputTag("caloStage2Digis", "Tau", "HLT2"),
    L1Tau = cms.InputTag("caloStage2Digis", "Tau", "RECO"),
    L1EmuTau = cms.InputTag("simCaloStage2Digis", "MP"),
    jetCollection = cms.InputTag("ak4PFJetsCHS"),
    l1tJetCollection = cms.InputTag("caloStage2Digis","Jet"),
    Vertexes = cms.InputTag("offlinePrimaryVertices"),
    triggerList = HLTLIST
)

TAndPseq = cms.Sequence(
    #hltFilter      +
    #goodMuons      +
    tauGenJets +
    trueHadronicTaus +
    goodTaus +
    genMatchedTaus 
)

NtupleSeq = cms.Sequence(
    Ntuplizer_noTagAndProbe
)

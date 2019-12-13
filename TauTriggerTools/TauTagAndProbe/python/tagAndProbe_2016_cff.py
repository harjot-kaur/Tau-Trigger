import FWCore.ParameterSet.Config as cms

print "Running on data"

# filter HLT paths for T&P
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt


HLTLIST_TAG = cms.VPSet(
    #MuTau SingleL1
    cms.PSet (
        HLT = cms.string("HLT_IsoMu27_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(-999)
    ),
)


HLTLIST = cms.VPSet(
    #MuTau 
 cms.PSet (
        HLT = cms.string("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTauJet20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu19LooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterIsoMu19LooseIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09", "hltOverlapFilterSingleIsoMu19LooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterSingleIsoMu19LooseIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTauJet20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu19LooseCombinedIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseCombinedIsoAgainstMuon", "hltOverlapFilterIsoMu19LooseCombinedIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erIsoTau26erL1f0L2f10QL3f19QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu19MediumIsoPFTau32Reg"),
        path2 = cms.vstring ("hltPFTau32TrackPt1MediumIsolationL1HLTMatchedReg", "hltOverlapFilterIsoMu19MediumIsoPFTau32Reg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erIsoTau26erL1f0L2f10QL3f19QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu19MediumCombinedIsoPFTau32Reg"),
        path2 = cms.vstring ("hltPFTau32TrackPt1MediumCombinedIsolationL1HLTMatchedReg", "hltOverlapFilterIsoMu19MediumCombinedIsoPFTau32Reg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erIsoTau26erL1f0L2f10QL3f19QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu19TightCombinedIsoPFTau32Reg"),
        path2 = cms.vstring ("hltPFTau32TrackPt1TightCombinedIsolationL1HLTMatchedReg", "hltOverlapFilterIsoMu19TightCombinedIsoPFTau32Reg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    # the following ones are extra!
    cms.PSet (
        HLT = cms.string("HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v"), 
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu20erIorSingleMu22erL1f0L2f10QL3f21QL3trkIsoFiltered0p09", "hltOverlapFilterSingleIsoMu21LooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterSingleIsoMu21LooseIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
   cms.PSet (
        HLT = cms.string("HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu20erIsoTau26erL1f0L2f10QL3f21QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu21MediumIsoPFTau32Reg"),
        path2 = cms.vstring ("hltPFTau32TrackPt1MediumIsolationL1HLTMatchedReg", "hltOverlapFilterIsoMu21MediumIsoPFTau32Reg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu20erIsoTau26erL1f0L2f10QL3f21QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu21MediumCombinedIsoPFTau32Reg"),
        path2 = cms.vstring ("hltPFTau32TrackPt1MediumCombinedIsolationL1HLTMatchedReg", "hltOverlapFilterIsoMu21MediumCombinedIsoPFTau32Reg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu20erIsoTau26erL1f0L2f10QL3f21QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu21TightCombinedIsoPFTau32Reg"),
        path2 = cms.vstring ("hltPFTau32TrackPt1TightCombinedIsolationL1HLTMatchedReg", "hltOverlapFilterIsoMu21TightCombinedIsoPFTau32Reg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
 )


hltFilter = hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths = ['HLT_IsoMu27_v*'],
    andOr = cms.bool(True), # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(True) #if True: throws exception if a trigger path is invalid
)

## only events where slimmedMuons has exactly 1 muon
muonNumberFilter = cms.EDFilter ("muonNumberFilter",
    src = cms.InputTag("slimmedMuons")
)

## good muons for T&P
goodMuons = cms.EDFilter("PATMuonRefSelector",
        src = cms.InputTag("slimmedMuons"),
        cut = cms.string(
         #       'pt > 5 && abs(eta) < 2.1 ' # kinematics
                'pt > 24 && abs(eta) < 2.1 ' # kinematics
                '&& ( (pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5 * pfIsolationR04().sumPUPt, 0.0)) / pt() ) < 0.1 ' # isolation
                '&& isMediumMuon()' # quality -- medium muon
        ),
        filter = cms.bool(True)
)



## good taus - apply analysis selection
goodTaus = cms.EDFilter("PATTauRefSelector",
        src = cms.InputTag("NewTauIDsEmbedded"),
        cut = cms.string(
        #        'pt > 5 && abs(eta) < 2.1 ' #kinematics
                'pt > 20 && abs(eta) < 2.1 ' #kinematics
                '&& abs(charge) > 0 && abs(charge) < 2 ' #sometimes 2 prongs have charge != 1
                '&& tauID("decayModeFinding") > 0.5 ' # tau ID
                #'&& tauID("byVVLooseIsolationMVArun2v1DBoldDMwLT") > 0.5 '
                '&& tauID("againstMuonTight3") > 0.5 ' # anti Muon tight
                '&& tauID("againstElectronVLooseMVA6") > 0.5 ' # anti-Ele loose
        ),
        filter = cms.bool(True)
)

## b jet veto : no additional b jets in the event (reject tt) -- use in sequence with
bjets = cms.EDFilter("PATJetRefSelector",
        src = cms.InputTag("slimmedJets"),
        cut = cms.string(
                'pt > 20 && abs(eta) < 2.4 ' #kinematics
                '&& bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.8484' # b tag with medium WP
        ),
        #filter = cms.bool(True)
)

TagAndProbe = cms.EDFilter("TauTagAndProbeFilter",
                           taus  = cms.InputTag("goodTaus"),
                           muons = cms.InputTag("goodMuons"),
                           met   = cms.InputTag("slimmedMETs"),
                           useMassCuts = cms.bool(False),
                           electrons = cms.InputTag("slimmedElectrons"),
						   eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wpLoose"),
                           eleVeto = cms.bool(True),
                           bjets = cms.InputTag("bjets")
)



patTriggerUnpacker = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
                                    patTriggerObjectsStandAlone = cms.InputTag("slimmedPatTrigger"),
                                    triggerResults = cms.InputTag('TriggerResults', '', "HLT"),
                                    unpackFilterLabels = cms.bool(True)
                                    )

Ntuplizer = cms.EDAnalyzer("Ntuplizer",
    treeName = cms.string("TagAndProbe"),
    isMC = cms.bool(False),
    genCollection = cms.InputTag(""),
    genPartCollection = cms.InputTag(""),
    muons = cms.InputTag("TagAndProbe"),
    taus = cms.InputTag("TagAndProbe"),
    puInfo = cms.InputTag("slimmedAddPileupInfo"), 
    met   = cms.InputTag("slimmedMETs"),
    triggerList = HLTLIST,
    triggerList_tag = HLTLIST_TAG,
    triggerSet = cms.InputTag("patTriggerUnpacker"),
    triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
    L1Tau = cms.InputTag("caloStage2Digis", "Tau", "RECO"),
    #L1EmuTau = cms.InputTag("simCaloStage2Digis"),
    L1EmuTau = cms.InputTag("simCaloStage2Digis", "MP"),
    Vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
    L2CaloJet_ForIsoPix_Collection = cms.InputTag("hltL2TausForPixelIsolation", "", "TEST"),
    L2CaloJet_ForIsoPix_IsoCollection = cms.InputTag("hltL2TauPixelIsoTagProducer", "", "TEST")   
)

TAndPseq = cms.Sequence(
    hltFilter        +
    muonNumberFilter +
    goodMuons        +
    goodTaus         +
    bjets            +
    TagAndProbe
)

NtupleSeq = cms.Sequence(
    patTriggerUnpacker +
    Ntuplizer
)

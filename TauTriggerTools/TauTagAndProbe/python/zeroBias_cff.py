import FWCore.ParameterSet.Config as cms


#print "Running on data or mc"

HLTLIST = cms.VPSet(
   #Mu-Tau20 (VBF monitoring)
    cms.PSet (
        HLT = cms.string("HLT_IsoMu27_LooseChargedIsoPFTau20_Trk1_eta2p1_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu27LooseChargedIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseChargedIsoAgainstMuon", "hltOverlapFilterIsoMu27LooseChargedIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu27_MediumChargedIsoPFTau20_Trk1_eta2p1_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu27MediumChargedIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackMediumChargedIsoAgainstMuon", "hltOverlapFilterIsoMu27MediumChargedIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu27_TightChargedIsoPFTau20_Trk1_eta2p1_SingleL_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu27TightChargedIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackTightChargedIsoAgainstMuon", "hltOverlapFilterIsoMu27TightChargedIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),

    #Mu-Tau35 (di-tau monitoring)
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24MediumChargedIsoPFTau35MonitoringReg"),
        path2 = cms.vstring ("hltSelectedPFTau35TrackPt1MediumChargedIsolationL1HLTMatchedReg", "hltOverlapFilterIsoMu24MediumChargedIsoPFTau35MonitoringReg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24MediumChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
        path2 = cms.vstring ("hltSelectedPFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsL1HLTMatchedReg", "hltOverlapFilterIsoMu24MediumChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24TightChargedIsoPFTau35MonitoringReg"),
        path2 = cms.vstring ("hltSelectedPFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg", "hltOverlapFilterIsoMu24TightChargedIsoPFTau35MonitoringReg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24TightChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
        path2 = cms.vstring ("hltSelectedPFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsL1HLTMatchedReg", "hltOverlapFilterIsoMu24TightChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    #Mu-Tau50 (Tau+MET monitoring)
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22erIsoTau40erL1f0L2f10QL3f24QL3trkIsoFiltered0p07"),
        path2 = cms.vstring ("hltSelectedPFTau50MediumChargedIsolationL1HLTMatchedMu22IsoTau40"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),

    #Mu-Tau27 (signal path)
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched", "hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27MediumChargedIsolationAgainstMuonL1HLTMatched", "hltOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27TightChargedIsolationAgainstMuonL1HLTMatched", "hltOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27LooseChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched", "hltOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27MediumChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched", "hltOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27TightChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched", "hltOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),


    #ETau CrossL1
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau", "hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30LooseChargedIsolationL1HLTMatched", "hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau", "hltOverlapFilterIsoEle24WPTightGsfMediumIsoPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30MediumChargedIsolationL1HLTMatched", "hltOverlapFilterIsoEle24WPTightGsfMediumIsoPFTau30"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau", "hltOverlapFilterIsoEle24WPTightGsfTightIsoPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30TightChargedIsolationL1HLTMatched", "hltOverlapFilterIsoEle24WPTightGsfTightIsoPFTau30"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau", "hltOverlapFilterIsoEle24WPTightGsfLooseIsoTightOOSCPhotonsPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30LooseChargedIsolationTightOOSCPhotonsL1HLTMatched", "hltOverlapFilterIsoEle24WPTightGsfLooseIsoTightOOSCPhotonsPFTau30"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau", "hltOverlapFilterIsoEle24WPTightGsfMediumIsoTightOOSCPhotonsPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30MediumChargedIsolationTightOOSCPhotonsL1HLTMatched", "hltOverlapFilterIsoEle24WPTightGsfMediumIsoTightOOSCPhotonsPFTau30"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau", "hltOverlapFilterIsoEle24WPTightGsfTightIsoTightOOSCPhotonsPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30TightChargedIsolationTightOOSCPhotonsL1HLTMatched", "hltOverlapFilterIsoEle24WPTightGsfTightIsoTightOOSCPhotonsPFTau30"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
    ),





    #Di-tau
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1TightChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1TightChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1TightChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1TightChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1TIghtChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1TIghtChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),

    #Tau + MET
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v"),
        path1 = cms.vstring ("hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v"),
        path1 = cms.vstring ("hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_v"),
        path1 = cms.vstring ("hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_v"),
        path1 = cms.vstring ("hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_v"),
        path1 = cms.vstring ("hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v"),
        path1 = cms.vstring ("hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),

    #SingleTau
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau100HighPtRelaxedIso_Trk50_eta2p1_1pr_v"),
        path1 = cms.vstring ("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched1Prong"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched1Prong"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau200MediumChargedIsolationL1HLTMatched1Prong"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau220MediumChargedIsolationL1HLTMatched1Prong"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),

    cms.PSet (
        HLT = cms.string("HLT_IsoMu27_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(-1)
    ),

    cms.PSet (
        HLT = cms.string("HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        leg1 = cms.int32(-1),
        leg2 = cms.int32(-1)
    ),

    cms.PSet (
        HLT = cms.string("HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        leg1 = cms.int32(-1),
        leg2 = cms.int32(-1)
    ),

    cms.PSet (
        HLT = cms.string("HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        leg1 = cms.int32(-1),
        leg2 = cms.int32(-1)
    ),



    #Mu+Tau HPS
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltHpsOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltHpsSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched", "hltHpsOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),    

        cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltHpsOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltHpsSelectedPFTau27MediumChargedIsolationAgainstMuonL1HLTMatched", "hltHpsOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),    
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltHpsOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltHpsSelectedPFTau27TightChargedIsolationAgainstMuonL1HLTMatched", "hltHpsOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),    
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltHpsOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltHpsSelectedPFTau27LooseChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched", "hltHpsOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),    
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltHpsOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltHpsSelectedPFTau27MediumChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched", "hltHpsOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),    
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltHpsOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltHpsSelectedPFTau27TightChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched", "hltHpsOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),    

    cms.PSet (
        HLT = cms.string("HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07", "hltHpsOverlapFilterIsoMu27LooseChargedIsoPFTau20"),
        path2 = cms.vstring ("hltHpsPFTau20TrackLooseChargedIsoAgainstMuon", "hltHpsOverlapFilterIsoMu27LooseChargedIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),

    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltHpsOverlapFilterIsoMu24MediumChargedIsoPFTau35MonitoringReg"),
        path2 = cms.vstring ("hltHpsSelectedPFTau35TrackPt1MediumChargedIsolationL1HLTMatchedReg", "hltHpsOverlapFilterIsoMu24MediumChargedIsoPFTau35MonitoringReg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),


    #E+Tau HPS

    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau", "hltHpsOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
        path2 = cms.vstring ("hltHpsSelectedPFTau30LooseChargedIsolationL1HLTMatched", "hltHpsOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau", "hltHpsOverlapFilterIsoEle24WPTightGsfMediumIsoPFTau30"),
        path2 = cms.vstring ("hltHpsSelectedPFTau30MediumChargedIsolationL1HLTMatched", "hltHpsOverlapFilterIsoEle24WPTightGsfMediumIsoPFTau30"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau", "hltHpsOverlapFilterIsoEle24WPTightGsfTightIsoPFTau30"),
        path2 = cms.vstring ("hltHpsSelectedPFTau30TightChargedIsolationL1HLTMatched", "hltHpsOverlapFilterIsoEle24WPTightGsfTightIsoPFTau30"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau", "hltHpsOverlapFilterIsoEle24WPTightGsfLooseIsoTightOOSCPhotonsPFTau30"),
        path2 = cms.vstring ("hltHpsSelectedPFTau30LooseChargedIsolationTightOOSCPhotonsL1HLTMatched", "hltHpsOverlapFilterIsoEle24WPTightGsfLooseIsoTightOOSCPhotonsPFTau30"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau", "hltHpsOverlapFilterIsoEle24WPTightGsfMediumIsoTightOOSCPhotonsPFTau30"),
        path2 = cms.vstring ("hltHpsSelectedPFTau30MediumChargedIsolationTightOOSCPhotonsL1HLTMatched", "hltHpsOverlapFilterIsoEle24WPTightGsfMediumIsoTightOOSCPhotonsPFTau30"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau", "hltHpsOverlapFilterIsoEle24WPTightGsfTightIsoTightOOSCPhotonsPFTau30"),
        path2 = cms.vstring ("hltHpsSelectedPFTau30TightChargedIsolationTightOOSCPhotonsL1HLTMatched", "hltHpsOverlapFilterIsoEle24WPTightGsfTightIsoTightOOSCPhotonsPFTau30"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
    ),

    #Di-tau HPS

    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),

    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau40TrackPt1MediumChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltHpsDoublePFTau40TrackPt1MediumChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau35TrackPt1TightChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltHpsDoublePFTau35TrackPt1TightChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau40TrackPt1TightChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltHpsDoublePFTau40TrackPt1TightChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltHpsDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltHpsDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltHpsDoublePFTau35TrackPt1TIghtChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau40TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltHpsDoublePFTau40TrackPt1TIghtChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),


    #VBF HPS

    cms.PSet (
        HLT = cms.string("HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        leg1 = cms.int32(-1),
        leg2 = cms.int32(-1)
    ),

    cms.PSet (
        HLT = cms.string("HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        leg1 = cms.int32(-1),
        leg2 = cms.int32(-1)
    ),

    cms.PSet (
        HLT = cms.string("HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        leg1 = cms.int32(-1),
        leg2 = cms.int32(-1)
    )



)





ZeroBias = cms.EDAnalyzer("ZeroBias",
    treeName = cms.string("ZeroBias"),
    L1Tau = cms.InputTag("caloStage2Digis", "Tau"),
    L1EmuTau = cms.InputTag("simCaloStage2Digis", "MP"),
    l1tJetCollection = cms.InputTag("caloStage2Digis","Jet"),
    l1tEmuJetCollection = cms.InputTag("simCaloStage2Digis","MP"),
    L1EG = cms.InputTag("caloStage2Digis", "EGamma"),
    L1EmuEG = cms.InputTag("simCaloStage2Digis", "MP"),
    L1Mu = cms.InputTag("hltGtStage2Digis"),
    L1EmuMu = cms.InputTag("simGtStage2Digis"),
    triggerList = HLTLIST,
    #triggerSet = cms.InputTag("slimmedPatTrigger"),
    triggerSet = cms.InputTag("patTriggerUnpacker"),
    triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
    L2CaloJet_L1TauSeeded_Collection = cms.InputTag("hltL2TauJetsL1IsoTauSeeded", "", "MYHLT"),
    L2CaloJet_ForIsoPix_Collection = cms.InputTag("hltL2TausForPixelIsolation", "", "MYHLT"),
    L2CaloJet_ForIsoPix_IsoCollection = cms.InputTag("hltL2TauPixelIsoTagProducer", "", "MYHLT"),                      
    L2CaloJet_IsoPix_Collection = cms.InputTag("hltL2TauJetsIso", "", "MYHLT"),
    PixelTrackCollection = cms.InputTag("hltPixelTracksMergedRegForTau", "", "MYHLT"),
    MergedTrackCollection = cms.InputTag("hltMergedTracksTauReg", "", "MYHLT"),
    PFRegCandCollection = cms.InputTag("hltParticleFlowReg", "", "MYHLT"),
    AK4PFRegJetCollection = cms.InputTag("hltAK4PFJetsReg", "", "MYHLT"),
    PFTauSansRefRegCollection = cms.InputTag("hltPFTausSansRefReg", "", "MYHLT"),
    PFJetRegionCollection = cms.InputTag("hltTauPFJets08RegionReg", "jets", "MYHLT"),
    PFJetChargedHadronAssociation = cms.InputTag("hltTauPFJetsRecoTauChargedHadronsReg", "", "MYHLT"),
    JetPiZeroAssociation = cms.InputTag("hltPFTauPiZerosReg", "", "MYHLT")
)


NtupleZeroBiasSeq = cms.Sequence(
    ZeroBias
)





patTriggerUnpacker = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
                                    patTriggerObjectsStandAlone = cms.InputTag("slimmedPatTrigger"),
                                    triggerResults = cms.InputTag('TriggerResults', '', "HLT"),
                                    unpackFilterLabels = cms.bool(True)
                                    )

patTriggerUnpackerSeq = cms.Sequence(
    patTriggerUnpacker
)

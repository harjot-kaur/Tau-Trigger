/*! Definition of a tuple with variables used as inputs to tau MVA isolation.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#pragma once

#include "TauTriggerTools/Common/interface/SmartTree.h"

#define VAR2(type, name1, name2) VAR(type, name1) VAR(type, name2)
#define VAR3(type, name1, name2, name3) VAR2(type, name1, name2) VAR(type, name3)
#define VAR4(type, name1, name2, name3, name4) VAR3(type, name1, name2, name3) VAR(type, name4)

#define MVA_TUPLE_DATA() \
    /* Event Variables */ \
    VAR(UInt_t, run) /* run number */ \
    VAR(UInt_t, lumi) /* lumi section */ \
    VAR(ULong64_t, evt) /* event number */ \
    /* Basic tau variables */ \
    VAR(Int_t, tau_index) /* index of the tau */ \
    VAR4(Float_t, tau_pt, tau_eta, tau_phi, tau_mass) /* 4-momentum of the tau */ \
    VAR(Int_t, tau_charge) /* tau charge */ \
    VAR(Int_t, lepton_gen_match) /* matching with leptons on the generator level (see Htautau Twiki for details):
                                    Electron = 1, Muon = 2, TauElectron = 3, TauMuon = 4, Tau = 5, NoMatch = 6 */ \
    /* Tau transverse impact paramters.
       See cmssw/RecoTauTag/RecoTau/plugins/PFTauTransverseImpactParameters.cc for details */ \
    VAR3(Float_t, tau_dxy_pca_x, tau_dxy_pca_y, tau_dxy_pca_z) /* The point of closest approach (PCA) of
                                                                  the leadPFChargedHadrCand to the primary vertex */ \
    VAR(Float_t, tau_dxy) /* tau signed transverse impact parameter wrt to the primary vertex */ \
    VAR(Float_t, tau_dxy_error) /* uncertainty of the transverse impact parameter measurement */ \
    VAR(Float_t, tau_ip3d) /* tau signed 3D impact parameter wrt to the primary vertex */ \
    VAR(Float_t, tau_ip3d_error) /* uncertainty of the 3D impact parameter measurement */ \
    VAR(Float_t, tau_dz) /* tau dz of the leadChargedHadrCand wrt to the primary vertex */ \
    VAR(Float_t, tau_dz_error) /* uncertainty of the tau dz measurement */ \
    VAR(Int_t, tau_hasSecondaryVertex) /* tau has the secondary vertex */ \
    VAR3(Float_t, tau_sv_x, tau_sv_y, tau_sv_z) /* position of the secondary vertex */ \
    VAR3(Float_t, tau_flightLength_x, tau_flightLength_y, tau_flightLength_z) /* flight length of the tau */ \
    VAR(Float_t, tau_flightLength_sig) /* significance of the flight length measurement */ \
    VAR(Float_t, tauPt) /*Harjot*/ \
    VAR(Float_t, tauEta) /*Harjot*/ \
    VAR(Float_t, charged_IsoPtSum) /*Harjot*/ \
    VAR(Float_t, neutral_IsoPtSum) /*Harjot*/ \
    VAR(Float_t, pu_CorrPtSum) /*Harjot*/ \
    VAR(Float_t, photon_PtSumOutsideSignalCone) /*Harjot*/ \
    VAR(Float_t, tauDecayMode) /*Harjot*/ \
    VAR(Float_t, n_Photon) /*Harjot*/ \
    VAR(Float_t, ptWeight_DetaStrip) /*Harjot*/ \
    VAR(Float_t, ptWeight_DphiStrip) /*Harjot*/ \
    VAR(Float_t, ptWeighted_DrSignal) /*Harjot*/ \
    VAR(Float_t, ptWeighted_DrIsolation) /*Harjot*/ \
    VAR(Float_t, e_Ratio) /*Harjot*/ \
    VAR(Float_t, taudxy) /*Harjot*/ \
    VAR(Float_t, taudxyabs) /*Harjot*/ \
    VAR(Float_t, taudxysig) /*Harjot*/ \
    VAR(Float_t, tauip3d) /*Harjot*/ \
    VAR(Float_t, tauip3dabs) /*Harjot*/ \
    VAR(Float_t, tauip3dsig) /*Harjot*/ \
    VAR(Float_t, hasSV) /*Harjot*/ \
    VAR(Float_t, decay_DistMag) /*Harjot*/ \
    VAR(Float_t, flightLengthSig) /*Harjot*/ \
    VAR(Float_t, gj_AngleDiff) /*Harjot*/ \

    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(tau_trigger, TauIdVars, MvaTuple, MVA_TUPLE_DATA, "taus")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(tau_trigger, MvaTuple, MVA_TUPLE_DATA)
#undef VAR
#undef VAR2
#undef VAR3
#undef VAR4
#undef MVA_TUPLE_DATA

namespace tau_trigger {

template<typename T>
constexpr T DefaultFillValue() { return std::numeric_limits<T>::lowest(); }
template<>
constexpr float DefaultFillValue<float>() { return -999.; }
template<>
constexpr int DefaultFillValue<int>() { return -999; }

} // namespace tau_trigger

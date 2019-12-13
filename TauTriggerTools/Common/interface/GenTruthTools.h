/*! Tools for working with MC generator truth.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#pragma once

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "AnalysisTypes.h"

namespace analysis {

namespace gen_truth {

struct LeptonMatchResult {
    GenLeptonMatch match{GenLeptonMatch::NoMatch};
    const reco::GenParticle* gen_particle{nullptr};
    std::vector<const reco::GenParticle*> visible_daughters;
    LorentzVectorXYZ visible_daughters_p4;
    int n_chargedParticles;
    int n_neutralParticles;
};

void FindFinalStateDaughters(const reco::GenParticle& particle, std::set<const reco::GenParticle*>& daughters,
                             const std::set<int>& pdg_to_exclude);

LorentzVectorXYZ GetFinalStateMomentum(const reco::GenParticle& particle, std::vector<const reco::GenParticle*>& visible_daughters,
                                       bool excludeInvisible, bool excludeLightLeptons);

LeptonMatchResult LeptonGenMatch(const LorentzVectorM& p4,
    const reco::GenParticleCollection& genParticles);


float GetNumberOfPileUpInteractions(edm::Handle<std::vector<PileupSummaryInfo>>& pu_infos);

} // namespace gen_truth
} // namespace analysis

#ifndef GENMATCHTAUFILTER_AOD_H
#define GENMATCHTAUFILTER_AOD_H

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/InputTag.h>

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include <TLorentzVector.h>

#include <iostream>

using namespace edm;
using namespace std;
 
class genMatchTauFilter_AOD : public edm::EDFilter {

    public:
        genMatchTauFilter_AOD(const edm::ParameterSet &);
        ~genMatchTauFilter_AOD();

private:
  bool filter(edm::Event &, edm::EventSetup const&);
  EDGetTokenT<reco::PFTauCollection>  _tauTag;
  EDGetTokenT<reco::GenJetCollection> _tauGenJetTag;
};

genMatchTauFilter_AOD::genMatchTauFilter_AOD(const edm::ParameterSet & iConfig) :
  _tauTag   (consumes<reco::PFTauCollection> (iConfig.getParameter<InputTag>("taus"))),
  _tauGenJetTag (consumes<reco::GenJetCollection> (iConfig.getParameter<InputTag>("genJets")))
{
    produces <reco::PFTauCollection>  ();
}

genMatchTauFilter_AOD::~genMatchTauFilter_AOD()
{}

bool genMatchTauFilter_AOD::filter(edm::Event & iEvent, edm::EventSetup const& iSetup)
{
    std::unique_ptr<reco::PFTauCollection>  resultTau  ( new  reco::PFTauCollection );
    Handle<reco::PFTauCollection> tauHandle;
    iEvent.getByToken (_tauTag, tauHandle);

    Handle<reco::GenJetCollection> tauGenJetHandle;
    iEvent.getByToken (_tauGenJetTag,tauGenJetHandle);

    int goodTaus = 0;
    int nTaus = 0;

    for(reco::PFTauCollection::const_iterator it=tauHandle->begin();	it!=tauHandle->end(); ++it) 
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiM(it->jetRef()->pt(),it->jetRef()->eta(),it->jetRef()->phi(),it->jetRef()->mass());

	bool matchedToGen = false;

	for(reco::GenJetCollection::const_iterator it2=tauGenJetHandle->begin(); it2!=tauGenJetHandle->end(); ++it2)
	  {
	    TLorentzVector genJet;
	    genJet.SetPtEtaPhiM(it2->pt(),it2->eta(),it2->phi(),it2->mass());
	    if(genJet.DeltaR(jet)<0.5)
	      {
		matchedToGen = true;
		break;
	      }
	  }

	TLorentzVector tau;
	tau.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());
	
	if (it->jetRef()->pt() > 8. && tau.DeltaR(jet)<0.5 && matchedToGen)
	  {
	    ++goodTaus;
	    resultTau->push_back ((*tauHandle)[nTaus]);
	  }
    
	nTaus++;

      }
    if (goodTaus == 0) return false;
    
    iEvent.put(std::move(resultTau));

    return true;
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(genMatchTauFilter_AOD);

#endif

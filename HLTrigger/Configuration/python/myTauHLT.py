import FWCore.ParameterSet.Config as cms

def update(process):
	process.options.wantSummary = cms.untracked.bool(True)
	return process

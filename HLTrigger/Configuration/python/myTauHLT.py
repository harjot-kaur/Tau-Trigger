import FWCore.ParameterSet.Config as cms
import TauTriggerTools.Common

def update(process,summary):
	process.options.wantSummary = cms.untracked.bool(summary)

	return process
#def update(process):
#	process.options.wantSummary = cms.untracked.bool(True)
#	return process

# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: tau_hlt --python_filename tau_hlt.py --customise=HLTrigger/Configuration/myTauHLT.update --step=DIGI,L1,DIGI2RAW,HLT:User --mc --no_exec --conditions auto:run2_mc_GRun --era=Run2_2018 --no_output -n 10 --filein /store/mc/RunIIFall18wmLHEGS/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/GEN-SIM/102X_upgrade2018_realistic_v11-v1/120000/D9393C36-D02E-8C44-948C-D50E889C6690.root
import FWCore.ParameterSet.Config as cms

#Harjot__start
import os
import re
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('wantSummary', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Print run summary at the end of the job.")
options.parseArguments()

#Harjot__end
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process('HLT',Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_User_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoTauTag.Configuration.loadRecoTauTagMVAsFromPrepDB_cfi')
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(30)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIIFall18wmLHEGS/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/GEN-SIM/102X_upgrade2018_realistic_v11-v1/120000/D9393C36-D02E-8C44-948C-D50E889C6690.root'),
    inputCommands = cms.untracked.vstring(
        'keep *', 
        'drop *_genParticles_*_*', 
        'drop *_genParticlesForJets_*_*', 
        'drop *_kt4GenJets_*_*', 
        'drop *_kt6GenJets_*_*', 
        'drop *_iterativeCone5GenJets_*_*', 
        'drop *_ak4GenJets_*_*', 
        'drop *_ak7GenJets_*_*', 
        'drop *_ak8GenJets_*_*', 
        'drop *_ak4GenJetsNoNu_*_*', 
        'drop *_ak8GenJetsNoNu_*_*', 
        'drop *_genCandidatesForMET_*_*', 
        'drop *_genParticlesForMETAllVisible_*_*', 
        'drop *_genMetCalo_*_*', 
        'drop *_genMetCaloAndNonPrompt_*_*', 
        'drop *_genMetTrue_*_*', 
        'drop *_genMetIC5GenJs_*_*'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
#	SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('tau_hlt nevts:1000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.TFileService = cms.Service('TFileService', fileName=cms.string("mvaTauTuple_DY.root"))
# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_GRun', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.myTauHLT
from HLTrigger.Configuration.myTauHLT import update 

#call to customisation function update imported from HLTrigger.Configuration.myTauHLT
process = update(process,options.wantSummary)

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
##print process.dumpPython()

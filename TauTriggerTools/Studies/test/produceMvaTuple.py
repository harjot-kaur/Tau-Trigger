import re
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from TauTriggerTools.Common.ProduceHelpers import *

options = VarParsing('analysis')
options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Text file with a list of the input root files to process.")
options.register('fileNamePrefix', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Prefix to add to input file names.")
options.register('outputTupleFile', 'eventTuple.root', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Event tuple file.")
options.register('skipEvents', -1, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "Number of events to skip")
options.register('eventList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "List of events to process.")
options.register('lumiFile', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "JSON file with lumi mask.")
options.register('period', 'Run2018', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Data taking period")
options.register('isMC', True, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Data or MC")
options.parseArguments()

processName = "TauAnalysis"
process = cms.Process(processName)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.GlobalTag.globaltag = getGlobalTag(options.period, options.isMC)
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring())
process.TFileService = cms.Service('TFileService', fileName=cms.string(options.outputTupleFile))

if len(options.inputFileList) > 0:
    readFileList(process.source.fileNames, options.inputFileList, options.fileNamePrefix)
elif len(options.inputFiles) > 0:
    addFilesToList(process.source.fileNames, options.inputFiles, options.fileNamePrefix)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
if options.maxEvents > 0:
    process.maxEvents.input = cms.untracked.int32(options.maxEvents)
if options.skipEvents > 0:
    process.source.skipEvents = cms.untracked.uint32(options.skipEvents)
if len(options.eventList) > 0:
    process.source.eventsToProcess = cms.untracked.VEventRange(options.eventList.split(','))
if len(options.lumiFile) > 0:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = options.lumiFile).getVLuminosityBlockRange()

process.mvaTupleProducer = cms.EDProducer("TauTriggerMvaTupleProducer",
    isMC            = cms.bool(options.isMC),
    genParticles    = cms.InputTag('prunedGenParticles'),
    taus            = cms.InputTag('slimmedTaus'),
)

process.p = cms.Path(
    process.mvaTupleProducer
)

# Verbosity customization
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = getReportInterval(process.maxEvents.input.value())
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

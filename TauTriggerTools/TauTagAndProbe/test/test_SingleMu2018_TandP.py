import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Config as cms
process = cms.Process("TagAndProbe")
import os

isMC = False
useGenMatch = False
useCustomHLT = False

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

#### handling of cms line options for tier3 submission
#### the following are dummy defaults, so that one can normally use the config changing file list by hand etc.

options = VarParsing.VarParsing ('analysis')
options.register ('skipEvents',
                  -1, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Number of events to skip")
options.register ('JSONfile',
                  "", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "JSON file (empty for no JSON)")
if not isMC:
 	options.outputFile = 'NTuple_SingleMu_Data_2018.root'
else:	
 	options.outputFile = 'NTuple_SingleMu_DYMC_2018.root'

options.inputFiles = []
options.maxEvents  = 1000
options.parseArguments()


def get_cmssw_version():
	"""returns 'CMSSW_X_Y_Z'"""
	return os.environ["CMSSW_RELEASE_BASE"].split('/')[-1]

def get_cmssw_version_number():
	"""returns 'X_Y_Z' (without 'CMSSW_')"""
	return map(int, get_cmssw_version().split("CMSSW_")[1].split("_")[0:3])

def versionToInt(release=9, subversion=4, patch=0):
	return release * 10000 + subversion * 100 + patch

def is_above_cmssw_version( release=10, subversion=2, patch=0):
	split_cmssw_version = get_cmssw_version_number()
	if versionToInt(release, subversion, patch) > versionToInt(split_cmssw_version[0], split_cmssw_version[1], split_cmssw_version[2]):
		return False
	else:
		return True

# START ELECTRON CUT BASED ID SECTION
#
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
#

# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")


#**********************
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
#**********************

process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
# overwrite a default parameter: for miniAOD, the collection name is a slimmed one
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')

from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)

# Define which IDs we want to produce
# Each of these two example IDs contains all four standard 
my_id_modules =[
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',  #Fall17 iso
] 


#Add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


egmMod = 'egmGsfElectronIDs'
mvaMod = 'electronMVAValueMapProducer'
setattr(process,egmMod,process.egmGsfElectronIDs.clone())
setattr(process,mvaMod,process.electronMVAValueMapProducer.clone())
process.electrons = cms.Sequence(getattr(process,mvaMod)*getattr(process,egmMod))


print "The current CMSSW version is", get_cmssw_version()

if is_above_cmssw_version(10, 1, 9):
	print "=== An electron MVA variable helper is added to the sequence for the electron ID ==="
	helpMod = 'electronMVAVariableHelper'
	setattr(process,helpMod,process.electronMVAVariableHelper.clone())
	process.electrons = cms.Sequence(getattr(process,helpMod)*getattr(process,mvaMod)*getattr(process,egmMod))


# =================================
# START TAU MVA BASED ID SECTION   -- NEW WAY --
# =================================

from TauTagAndProbe.TauTagAndProbe.runTauIdMVA import TauIDEmbedder

toKeep = []
toKeep.extend(("2017v1","2017v2","newDM2017v2", "dR0p32017v2"))

na = TauIDEmbedder(process, cms,
    debug=True,
    toKeep = toKeep
)
na.runTauID()


if not isMC:
    from Configuration.AlCa.autoCond import autoCond
    #process.GlobalTag.globaltag = '102X_dataRun2_Sep2018Rereco_v1' # for 2018 ReReco RunABC samples
    process.GlobalTag.globaltag = '102X_dataRun2_Sep2018ABC_v2' # for 2018 ReReco RunABC samples
    #process.GlobalTag.globaltag = '102X_dataRun2_Prompt_v13' # for 2018 PromtReco RunD samples
    process.load('TauTagAndProbe.TauTagAndProbe.tagAndProbe_cff')
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
          #  '/store/data/Run2017E/SingleMuon/MINIAOD/17Nov2017-v1/50000/000DCB8B-2ADD-E711-9100-008CFAF35AC0.root'
         #  '/store/data/Run2017B/SingleMuon/MINIAOD/31Mar2018-v1/100000/001642F1-6638-E811-B4FA-0025905B857A.root'
         #'/store/data/Run2018A/SingleMuon/MINIAOD/06Jun2018-v1/410000/F2A0DDD5-FF83-E811-A183-FA163EFE9CA3.root'
         '/store/data/Run2018D/SingleMuon/MINIAOD/PromptReco-v2/000/321/988/00000/6AB32549-34AF-E811-A47F-FA163EBD19B5.root'
        ),
    )


else:
    process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v18'
    process.load('TauTagAndProbe.TauTagAndProbe.MCanalysis_cff')
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(            
        #   '/store/mc/RunIIFall17MiniAOD/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/12DB4A06-65D7-E711-8DA4-0CC47A78A456.root'
	#'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/20000/00D13F2E-6F44-E811-923E-001E0BED0560.root'
	#' /store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/52E26FAC-370D-0742-96AF-3E5439F7CD9C.root'
	'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/FFDCFC59-4ABE-0646-AABE-BD5D65301169.root',
	'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/FC40474D-EE5A-C845-AEC3-E9B088F15648.root',
	'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/FA870652-56D3-1E45-B1EE-EF8A8FC89506.root',
	'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/FA553171-E5E9-5346-B634-6A64D6DEB64C.root',
	'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/FA553171-E5E9-5346-B634-6A64D6DEB64C.root'
        )
    )


if useCustomHLT:
    process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","MYHLT")
    process.Ntuplizer.triggerSet = cms.InputTag("selectedPatTriggerCustom", "", "MYHLT")
    process.Ntuplizer.triggerResultsLabel = cms.InputTag("TriggerResults", "", "MYHLT")
    process.Ntuplizer.L2CaloJet_ForIsoPix_Collection = cms.InputTag("hltL2TausForPixelIsolation", "", "MYHLT")
    process.Ntuplizer.L2CaloJet_ForIsoPix_IsoCollection = cms.InputTag("hltL2TauPixelIsoTagProducer", "", "MYHLT")


if isMC and useGenMatch:
    process.Ntuplizer.taus = cms.InputTag("genMatchedTaus")


if options.JSONfile:
    print "Using JSON: " , options.JSONfile
    process.source.lumisToProcess = LumiList.LumiList(filename = options.JSONfile).getVLuminosityBlockRange()

if options.inputFiles:
    process.source.fileNames = cms.untracked.vstring(options.inputFiles)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

if options.maxEvents >= -1:
    process.maxEvents.input = cms.untracked.int32(options.maxEvents)
if options.skipEvents >= 0:
    process.source.skipEvents = cms.untracked.uint32(options.skipEvents)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)


process.p = cms.Path(
    process.electrons +
	process.rerunMvaIsolationSequence +
	getattr(process, "NewTauIDsEmbedded") +
    process.TAndPseq +
    process.NtupleSeq
)

if isMC and useGenMatch:
    process.p = cms.Path(
    	process.electrons +
		process.rerunMvaIsolationSequence +
		getattr(process, "NewTauIDsEmbedded") +
    	process.TAndPseq +
	process.genMatchedSeq +
   	process.NtupleSeq
    )

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Adding ntuplizer
process.TFileService=cms.Service('TFileService',fileName=cms.string(options.outputFile))

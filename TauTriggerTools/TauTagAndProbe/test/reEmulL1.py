import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

isMC = True

process = cms.Process("TagAndProbe",eras.Run2_2016)
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')


options = VarParsing.VarParsing ('analysis')
options.register ('secondaryFilesList','',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,  "List of secondary input files")

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
options.outputFile = 'NTuple.root'
options.inputFiles = []
options.maxEvents  = -999

options.parseArguments()

import FWCore.Utilities.FileUtils as FileUtils
listSecondaryFiles = FileUtils.loadListFromFile (options.secondaryFilesList)

if not isMC: # will use 80X
    from Configuration.AlCa.autoCond import autoCond
    process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v8'
    process.load('TauTagAndProbe.TauTagAndProbe.tagAndProbe_cff')
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
            '/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/274/199/00000/7005DB70-4C28-E611-8628-02163E0144DD.root',
        ),

        secondaryFileNames = cms.untracked.vstring(listSecondaryFiles)
        #secondaryFileNames = cms.untracked.vstring('file:AA918EB1-6E64-E611-9BE0-00259074AE54.root')
            #'/store/mc/RunIISpring16MiniAODv2/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/FlatPU20to70HcalNZSRAW_withHLT_80X_mcRun2_asymptotic_v14-v1/50000/302E52FC-8567-E611-B2AA-0CC47A703326.root',
            #),
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/6E300626-5E26-E611-980B-02163E0119A2.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/08F68A47-5D26-E611-B042-02163E012239.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/1011D344-5E26-E611-ABC4-02163E011CF0.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/1E3EFD43-5E26-E611-AE17-02163E0146FF.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/2EA6473B-5E26-E611-B82A-02163E011EAC.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/3CD7BB24-5D26-E611-88A3-02163E014736.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/3CE74444-5E26-E611-8CE6-02163E012545.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/488D1A44-5E26-E611-8057-02163E014285.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/5EECA846-5D26-E611-A99C-02163E01432B.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/7C7FB848-5E26-E611-8A06-02163E014167.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/86E8DC25-5D26-E611-9827-02163E014713.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/A4B22F44-5E26-E611-AFDA-02163E0141F3.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/BA326B3B-5E26-E611-A6E0-02163E0124FA.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/BADC0417-5D26-E611-872F-02163E012A7E.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/C813CA24-5E26-E611-B7A2-02163E011F93.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/E251C323-5E26-E611-825D-02163E011A0F.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/E8FD1844-5E26-E611-B99E-02163E0146CB.root',
#            '/store/data/Run2016B/SingleMuon/RAW/v2/000/274/199/00000/F8B44816-5E26-E611-A87A-02163E011E74.root',
#        ) 
    )

    #process.source.eventsToProcess = cms.untracked.VEventRange('281613:108:12854629')
    #process.source.eventsToProcess = cms.untracked.VEventRange('274199:353:670607108')

else:
    #process.GlobalTag.globaltag = 'auto:run2_mc' #MC 25 ns miniAODv2
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_v14'
    #process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2' #MC 25 ns miniAODv2
    # process.GlobalTag.globaltag = '76X_dataRun2_16Dec2015_v0'
    process.load('TauTagAndProbe.TauTagAndProbe.MCanalysis_cff')
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
            #'/store/mc/RunIISpring16DR80/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/FlatPU20to70HcalNZSRAW_withHLT_80X_mcRun2_asymptotic_v14-v1/40000/AA918EB1-6E64-E611-9BE0-00259074AE54.root'
            '/store/mc/RunIISpring16MiniAODv2/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/FlatPU20to70HcalNZSRAW_withHLT_80X_mcRun2_asymptotic_v14-v1/50000/1A13CB76-9B67-E611-A143-0050560210EC.root'
        ),
       secondaryFileNames = cms.untracked.vstring(listSecondaryFiles)
            #'/store/mc/RunIISummer15GS/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM/MCRUN2_71_V1-v1/00000/08D2C535-5458-E511-B0C0-FA163E83549A.root'
            #'/store/mc/RunIISummer15GS/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM/MCRUN2_71_V1-v1/10000/8A2D3925-4658-E511-80B2-02163E014126.root',
            #'/store/mc/RunIISummer15GS/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM/MCRUN2_71_V1-v1/10000/ECB7FC03-8058-E511-BE9D-02163E0141A2.root',
            #'/store/mc/RunIISummer15GS/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM/MCRUN2_71_V1-v1/60000/4473DFF5-9456-E511-9C4B-002590494C8A.root',
            #'/store/mc/RunIISummer15GS/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM/MCRUN2_71_V1-v1/60000/587A9A92-A956-E511-AC52-0025904B11CC.root',
            #'/store/mc/RunIISummer15GS/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM/MCRUN2_71_V1-v1/60000/7076052C-4A57-E511-B371-00259074AE9A.root',
            #'/store/mc/RunIISummer15GS/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM/MCRUN2_71_V1-v1/60000/7CD6A6C7-9C56-E511-87E5-003048C75840.root',
            #'/store/mc/RunIISummer15GS/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM/MCRUN2_71_V1-v1/60000/7E281C8E-4357-E511-8498-00259074AE80.root',
            #'/store/mc/RunIISummer15GS/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM/MCRUN2_71_V1-v1/60000/98884797-4957-E511-82C6-00259073E504.root',
            #'/store/mc/RunIISummer15GS/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM/MCRUN2_71_V1-v1/60000/BE274E99-4957-E511-88ED-0025907A1A2E.root',

    )

process.schedule = cms.Schedule()

## L1 emulation stuff

if not isMC:
    from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAW 
    process = L1TReEmulFromRAW(process)
else:
    from L1Trigger.Configuration.customiseReEmul import L1TReEmulMCFromRAW
    process = L1TReEmulMCFromRAW(process) 
    from L1Trigger.Configuration.customiseUtils import L1TTurnOffUnpackStage2GtGmtAndCalo 
    process = L1TTurnOffUnpackStage2GtGmtAndCalo(process)

process.load("L1Trigger.L1TCalorimeter.caloStage2Params_2016_v3_2_cfi")

#### handling of cms line options for tier3 submission
#### the following are dummy defaults, so that one can normally use the config changing file list by hand etc.



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

process.p = cms.Path (
    process.TAndPseq  +
    process.RawToDigi +
    process.L1TReEmul +
    process.NtupleSeq
)
process.schedule = cms.Schedule(process.p) # do my sequence pls

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Adding ntuplizer
process.TFileService=cms.Service('TFileService',fileName=cms.string(options.outputFile))

# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TagAndProbe_SingleMu297050'
config.General.workArea = 'DefaultCrab3Area'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test_SingleMu297050.py'

config.section_("Data")
config.Data.inputDataset = '/SingleMuon/Run2017B-PromptReco-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 10000#number of events per jobs
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase = '/store/user/tstreble/TagAndProbeTrees'
config.Data.publication = False
config.Data.outputDatasetTag = 'TagAndProbe_SingleMu297050'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.txt'
# json with 3.99/fb
config.Data.runRange = '297050'

config.section_("Site")
config.Site.storageSite = 'T2_FR_GRIF_LLR'


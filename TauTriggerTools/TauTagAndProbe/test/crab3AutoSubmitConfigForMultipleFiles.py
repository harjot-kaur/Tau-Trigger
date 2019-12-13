# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()
from CRABClient.UserUtilities import getUsernameFromSiteDB
from multiprocessing import Process
from CRABAPI.RawCommand import crabCommand



import datetime
today=datetime.date.today().strftime("%Y-%m-%d")
date=today

def submit(config):
	crabCommand('submit', config = config)
        #except HTTPException as hte:
         #       print "Failed submitting task: %s" % (hte.headers)
        #except ClientException as cle:
         #       print "Failed submitting task: %s" % (cle)


config.section_("General")

config.General.workArea = 'DefaultCrab3Area'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.disableAutomaticOutputCollection = False
	
config.section_("Data")
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic' #EventAwareLumiBased'
#config.Data.unitsPerJob = 180 #number of events per jobs
config.Data.totalUnits = -1  #number of event
config.Data.outLFNDirBase = '/store/user/%s/trigger/TagAndProbeTrees/%s'%(getUsernameFromSiteDB(), date)
config.Data.publication = False
config.Data.allowNonValidInputDataset = True

config.section_("Site")
config.Site.storageSite = 'T2_DE_RWTH'


samples2016New = ["/SingleMuon/Run2016B-17Jul2018_ver2-v1/MINIAOD", "/SingleMuon/Run2016C-17Jul2018-v1/MINIAOD", "/SingleMuon/Run2016D-17Jul2018-v1/MINIAOD", "/SingleMuon/Run2016E-17Jul2018-v1/MINIAOD", "/SingleMuon/Run2016F-17Jul2018-v1/MINIAOD", "/SingleMuon/Run2016G-17Jul2018-v1/MINIAOD", "/SingleMuon/Run2016H-17Jul2018-v1/MINIAOD"]
samples2017 = ["/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD","/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD","/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD","/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD","/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD"]
sampleReReco2018 = ["/SingleMuon/Run2018A-17Sep2018-v2/MINIAOD","/SingleMuon/Run2018B-17Sep2018-v1/MINIAOD","/SingleMuon/Run2018C-17Sep2018-v1/MINIAOD"]
sample2018 = ["/SingleMuon/Run2018D-PromptReco-v2/MINIAOD"]  # different global tags are used for promptreco and rereco, so they need to be run separately

samples2016MC = ["/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM"]
samples2017MC = ["/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM","/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM","/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM"]
samples2018MC = ['/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM']


# ***** Choose which samples do you want to submit!*****
samples = ""
#samples = samples2016MC

if (samples == ""):
        print "===================================================================="
        print " Please choose which samples do you want to submit within the code! "
        print "===================================================================="


for index, sample in enumerate(samples):  

	print "ampleS", sample
	if("Run201" in sample):
		samplenickname = (sample.split("/",1)[-1]).replace("/","_") 
		config.General.requestName = samplenickname[:100] + "_" + date 
		print "here", samplenickname[:100] + "_" + date 
	elif("Jets" in sample):
		samplenickname = (sample.split("/",1)[-1]).replace("/","_") 
		config.General.requestName = samplenickname[:80] + "_" + date 
		print "here", samplenickname[:100] + "_" + date 
	# pSet files are given here depending on the 2016 or 2017 samples
	if "2016" in sample or "Summer16" in sample:
		config.JobType.psetName = 'test_SingleMu2016_TandP.py'
	elif "2017" in sample or "Fall17" in sample:
		config.JobType.psetName = 'test_SingleMu297050_TandP.py'
	elif "2018" in sample or "Autumn18" in sample:
		config.JobType.psetName = 'test_SingleMu2018_TandP.py'
	
	#config.JobType.pyCfgParams = ['outputfilename=%s.root'%(samplenickname)]
	#config.JobType.outputFiles = ['%s.root'%(samplenickname)]
	
	config.Data.inputDataset = sample
	config.Data.outputDatasetTag = 'TagAndProbe_SingleMu_' + sample.split("/")[2]
	
	# lumiMasks for Data samples
	if "Run2016" in sample:
		config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
	elif "Run2017" in sample:
		config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
	elif "Run2018" in sample:
		config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
	
	# RunRange needs to be given if not all run of the given era are not included in the JSON file 
	#if "Run2017F" in sample:
	#	config.Data.runRange = '305040-306462'


	print "========================================================================="
	print "The sample[",index,"]:" , sample, " is being submitted to the crab"
	print "========================================================================="
	
	p = Process(target=submit, args=(config,))
	p.start()
	p.join()




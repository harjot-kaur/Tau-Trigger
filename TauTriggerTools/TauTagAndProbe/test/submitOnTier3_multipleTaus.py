import os

isMC = True
#isMC = False

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

def splitInBlocks (l, n):
    """split the list l in n blocks of equal size"""
    k = len(l) / n
    r = len(l) % n

    i = 0
    blocks = []
    while i < len(l):
        if len(blocks)<r:
            blocks.append(l[i:i+k+1])
            i += k+1
        else:
            blocks.append(l[i:i+k])
            i += k

    return blocks

###########

njobs = 200
#filelist = open("fileListCertification_20Oct16.txt")
#filelist = open("fileListCertification_27Oct16.txt")
#filelist = open("fileListCertification_3Nov16.txt")
#filelist = open("fileList_Run2016B-23Sep2016-v3_141116.txt")
#filelist = open("fileList_Run2016C-23Sep2016-v1_141116.txt")
#filelist = open("fileList_Run2016D-23Sep2016-v1_141116.txt")
#filelist = open("fileList_Run2016E-23Sep2016-v1_141116.txt")
#filelist = open("fileList_Run2016E-23Sep2016-v1_141116.txt")
#filelist = open("fileList_Run2016F-23Sep2016-v1_141116.txt")
#filelist = open("fileList_Run2016G-23Sep2016-v1_141116.txt")
#filelist = open("fileList_Run2016H-PromptReco-v2_141116.txt")
#filelist = open("fileList_Run2016H-PromptReco-v3_141116.txt")
#filelist = open("fileList_Run2016H-PromptReco-v2_141116_Run282092.txt")
filelist = open("VBFHToTauTau_M125_13TeV_powheg_pythia8_RunIISummer16MiniAODv2-FlatPU28to62HcalNZSRAW_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1.txt")

#filelist = open("fileList_MC_RECO.txt")
#folder = "MC_RECO_9x9"
#folder = "MC_RECO_12x12"
#folder = "Certification_20Oct16"
#folder = "Certification_27Oct16"

folder = "MC_MiniAOD_multipleTaus_16_05_17"
#folder = "Run2016H-PromptReco-v2_141116_282092_noBtagVeto"
#folder = "Run2016H-PromptReco-v3_141116"
#folder = "Run2016H-PromptReco-v2_141116"
#folder = "Run2016G-23Sep2016-v1_141116"
#folder = "Run2016F-23Sep2016-v1_141116"
#folder = "Run2016E-23Sep2016-v1_141116"
#folder = "Run2016D-23Sep2016-v1_141116"
#folder = "Run2016C-23Sep2016-v1_141116"
#folder = "Run2016B-23Sep2016-v3_141116"
JSONfile = "/home/llr/cms/davignon/json_NoL1T.txt"
#JSONfile = "/home/llr/cms/davignon/json_DCSONLY.txt"
#/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.txt

###########

os.system ('source /opt/exp_soft/cms/t3/t3setup')

os.system('mkdir ' + folder)
files = [f.strip() for f in filelist]
print "Input has" , len(files) , "files" 
if njobs > len(files) : njobs = len(files)
filelist.close()

fileblocks = splitInBlocks (files, njobs)

for idx, block in enumerate(fileblocks):
    outRootName = folder + '/Ntuple_' + str(idx) + '.root'
    outJobName  = folder + '/job_' + str(idx) + '.sh'
    outListName = folder + "/filelist_" + str(idx) + ".txt"
    outLogName  = os.getcwd() + "/" + folder + "/log_" + str(idx) + ".txt"

    jobfilelist = open(outListName, 'w')
    for f in block: jobfilelist.write(f+"\n")
    jobfilelist.close()

    if not isMC:
        cmsRun = "cmsRun test.py maxEvents=-1 inputFiles_load="+outListName + " outputFile="+outRootName + " JSONfile="+JSONfile + " >& " + outLogName
    else:
        cmsRun = "cmsRun test_noTagAndProbe_multipleTaus.py maxEvents=-1 inputFiles_load="+outListName + " outputFile="+outRootName + " >& " + outLogName        

    skimjob = open (outJobName, 'w')
    skimjob.write ('#!/bin/bash\n')
    skimjob.write ('export X509_USER_PROXY=~/.t3/proxy.cert\n')
    skimjob.write ('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
    skimjob.write ('cd %s\n' % os.getcwd())
    skimjob.write ('export SCRAM_ARCH=slc6_amd64_gcc472\n')
    skimjob.write ('eval `scram r -sh`\n')
    skimjob.write ('cd %s\n'%os.getcwd())
    skimjob.write (cmsRun+'\n')
    skimjob.close ()

    os.system ('chmod u+rwx ' + outJobName)
    command = ('/opt/exp_soft/cms/t3/t3submit_new -long \'' + outJobName +"\'")
    # print command
    os.system (command)

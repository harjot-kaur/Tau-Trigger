import os
import json
from subprocess import Popen, PIPE

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
filelist = open("Data_SingleMu_2016_All_281613_282037.txt")
#filelist = open("Data_SingleMu_2016RunB_PromptRecov2_1Luglio.txt")

folder = "Certification_281613_282037"
#folder = "testSubmitT3TAndP2Luglio"

JSONfile = "/home/llr/cms/davignon/json_DCSONLY.txt"
#JSONfile = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.txt"

###########

os.system ('source /opt/exp_soft/cms/t3/t3setup')

os.system('mkdir ' + folder)
files = [f.strip() for f in filelist]
print "Input has" , len(files) , "files" 
if njobs > len(files) : njobs = len(files)
filelist.close()

fileblocks = splitInBlocks (files, njobs)

from das_client import get_data

for idx, block in enumerate(fileblocks):
    #print idx, block
    secondaryFileName = folder+"/secondaryFileList_split_"+str(idx)+".txt"
    f = open(secondaryFileName, 'w')
    myList = list()

    for currentFile in block:
        #print "finding parent of: ",currentFile
        command_das = "das_client --query=\"parent file="+currentFile+"\" --limit=0"
        #print "    > command is: " , command_das

        pipe = Popen(command_das, shell=True, stdout=PIPE)        
        for line in pipe.stdout:
            line = line.rstrip("\n")
            myList.append(line)
            #print >> f, line

    unique = []
    [unique.append(item) for item in myList if item not in unique]
    #print unique

    for line in unique:
        print >> f, line

    outRootName = folder + '/Ntuple_' + str(idx) + '.root'
    outJobName  = folder + '/job_' + str(idx) + '.sh'
    outListName = folder + "/filelist_" + str(idx) + ".txt"
    outLogName  = os.getcwd() + "/" + folder + "/log_" + str(idx) + ".txt"

    jobfilelist = open(outListName, 'w')
    for f in block: jobfilelist.write(f+"\n")
    jobfilelist.close()

    cmsRun = "cmsRun reEmulL1.py maxEvents=-1 inputFiles_load="+outListName + " secondaryFilesList="+secondaryFileName+" " + " outputFile="+outRootName + " JSONfile="+JSONfile + " >& " + outLogName
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
#    command = ('/opt/exp_soft/cms/t3/t3submit_new -short -q cms \'' + outJobName +"\'")
    print command
    os.system (command)

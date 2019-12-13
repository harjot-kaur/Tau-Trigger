from ROOT import *
import numpy as n

MC = True
DYJets = True   # False means WJet enriched cuts will be used, True means cuts for DYJet enriched samples will be used


if MC:
	saveOnlyOS = True # True; save only OS, False; save both and store weight for bkg sub
	if DYJets: tauGenMatching = True
	if not DYJets: tauGenMatching = False
	excludeLumiSections = False
	print "==> OS events are stored and tau gen matching is applied for MC samples! <=="
else:
	saveOnlyOS = False # True; save only OS, False; save both and store weight for bkg sub
	tauGenMatching = False
	excludeLumiSections = False
	print "==> SS events are stored as weights and applied to suppress the bkg for Data samples! <=="

# the hadd of all the output ntuples
#path = "/eos/user/h/hsert/TriggerStudies/ForkedRepo/Samples_AFS/Samples2016/syncronisedSamples/2017_11_23/"
path = "/eos/user/h/hsert/TriggerStudies/ForkedRepo/Samples2016_Legacy94X/190226/"

#fname = "./NTuple_SingleMu_DYMC_2016.root"
fname = path + "Ntuple_DYJets_RunIISummer16MiniAODv3_94X_mcRun2_ext1ANDext2-v2_190306_PUreweight.root"
#fname =  path +"Ntuple_SingleMuon_Run2016BtoH-17Jul2018_190228.root" 
#fname = path + "Ntuple_SingleMuon_Run2016G-17Jul2018-v1_190226.root"


#pt = [20, 26, 30, 34]
pt = [20, 26, 30, 32, 34]
numberOfHLTTriggers = 10   # the last one for an extra path to be filled with the events passes trigger=3 for 2016G corresponding mediumiso and trigger=4 for 2016H for McombinedIso
numberOfHLTTriggersForFit = 14
	
saveOnlyOS = False # True; save only OS, False: save both and store weight for bkg sub
print fname
#######################################################
fIn = TFile.Open(fname)
tIn = fIn.Get('Ntuplizer/TagAndProbe')
tTriggerNames = fIn.Get("Ntuplizer/triggerNames")
if MC:
	suppressionType = "OStauGenMatched"
else:
	suppressionType = "SSsubtraction"
	
if DYJets:
	outname = fname.replace ('.root', '_' + suppressionType + '_VVLooseWP2017v2_forFit.root')
else:
	outname = fname.replace ('.root', '_' + suppressionType + '_WjetEnriched_MediumWP2017v2_forFit.root')
	
fOut = TFile (outname, 'recreate')
tOut = tIn.CloneTree(0)
tOutNames = tTriggerNames.CloneTree(-1) # copy all

briso   = [n.zeros(1, dtype=int) for x in range (0, len(pt))]
brnoiso = [n.zeros(1, dtype=int) for x in range (0, len(pt))]
bkgSubW = n.zeros(1, dtype=float)
bkgSubANDpuW = n.zeros(1, dtype=float)
#tauPt_35GeV = n.zeros(1, dtype=bool)
#tauPt_27GeV = n.zeros(1, dtype=bool)

hltPathTriggered_OS   = [n.zeros(1, dtype=int) for x in range (0, numberOfHLTTriggersForFit+1)]

for i in range (0, len(pt)):
	name = ("hasL1_" + str(pt[i]))
	tOut.Branch(name, brnoiso[i], name+"/I")
	name += "_iso"
	tOut.Branch(name, briso[i], name+"/I")
	
for i in range (0, numberOfHLTTriggersForFit):
	tTriggerNames.GetEntry(i)
	if(i < numberOfHLTTriggers):
		name = ("hasHLTPath_" + str(i))
	elif(i==numberOfHLTTriggers):
		name = ("hasHLTmutauPath_0") #HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v
	elif(i==numberOfHLTTriggers+1):
		name = ("hasHLTmutauPath_1") # HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v
	elif(i==numberOfHLTTriggers+2):
		name = ("hasHLTetauPath_0and1")# mutau_0and1_and_plusL1Tau26andHLTTau30")
	elif(i==numberOfHLTTriggers+3):
		name = ("hasHLTditauPath_3or4") #HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_v  ===== HLT_3=> HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_v for MC

	tOut.Branch(name, hltPathTriggered_OS[i], name+"/I")

tOut.Branch("bkgSubW", bkgSubW, "bkgSubW/D")
tOut.Branch("bkgSubANDpuW", bkgSubANDpuW, "bkgSubANDpuW/D")


nentries = tIn.GetEntries()
for ev in range (0, nentries):
	tIn.GetEntry(ev)
	if (ev%10000 == 0) : print ev, "/", nentries

	if abs(tIn.tauEta) > 2.1:
		continue
	
	PS_column = tIn.PS_column
	RunNumber = tIn.RunNumber
	
	if saveOnlyOS and not tIn.isOS:
		continue
		
        for i in range (0, len(pt)):
                briso[i][0] = 0
		brnoiso[i][0] = 0

	for i in range (0, numberOfHLTTriggers):
		hltPathTriggered_OS[i][0] = 0
		
	L1iso = True if tIn.l1tIso > 0 else False
	L1pt = tIn.l1tPt
	for i in range(0, len(pt)):
		# print L1pt, pt[i]
		#
		if L1pt > pt[i]:	
			brnoiso[i][0] = 1
        	# print "SUCCESS!! ", brnoiso[i]
        	if L1iso:
        		briso[i][0] = 1
   
	
	triggerBits = tIn.tauTriggerBits
	HLTpt = tIn.hltPt

	for bitIndex in range(0, numberOfHLTTriggers):
		if bitIndex==0:
			if ((triggerBits >> bitIndex) & 1) == 1: # and (HLTpt>=27):  # this is the path for mutau trigger. So no extra requirement is needed like: L1pt and L1iso and HLTpt
				hltPathTriggered_OS[bitIndex][0] = 1
				hltPathTriggered_OS[numberOfHLTTriggers][0] = 1
			else:
				hltPathTriggered_OS[bitIndex][0] = 0
				hltPathTriggered_OS[numberOfHLTTriggers][0] = 0
		elif bitIndex==1:
			if ((triggerBits >> bitIndex) & 1) == 1: # and (HLTpt>=27): # this is the path for mutau trigger. So no extra requirement is needed like: L1pt and L1iso and HLTpt
				hltPathTriggered_OS[bitIndex][0] = 1
				hltPathTriggered_OS[numberOfHLTTriggers+1][0] = 1
			else:
				hltPathTriggered_OS[bitIndex][0] = 0
				hltPathTriggered_OS[numberOfHLTTriggers+1][0] = 0
		else:
			if ((triggerBits >> bitIndex) & 1) == 1: 
				hltPathTriggered_OS[bitIndex][0] = 1
			else:
				hltPathTriggered_OS[bitIndex][0] = 0
		
	if not MC:
		if ((((RunNumber < 276215 and ((triggerBits >> 1) & 1) == 1)) or (RunNumber > 276215 and RunNumber < 278270 and ((triggerBits >> 0) & 1) == 1)) or (RunNumber > 278270 and ((triggerBits >> 0) & 1) == 1 and HLTpt[0] >=30 and L1pt >= 26 and L1iso>0)):
			hltPathTriggered_OS[numberOfHLTTriggers+2][0] = 1  # this is the path for e-tau trigger. HLTpt and L1 cuts are required to have the same threshold on tau. Different mutau paths are used for etau meausurment for different runs
		else:
			hltPathTriggered_OS[numberOfHLTTriggers+2][0] = 0	
       	elif MC:
		if ((((triggerBits >> 1) & 1) == 1)):
			hltPathTriggered_OS[numberOfHLTTriggers+2][0] = 1  # this is the path for e-tau trigger. No HLTpt and L1 cuts are required to have the same threshold on tau. Only one mutau path is used for etau meausurments in MC, since OR also gives the same efficiencies
		else:
			hltPathTriggered_OS[numberOfHLTTriggers+2][0] = 0
			
	if not MC:
		if ((("2016H" in fname) and ((triggerBits >> 4) & 1) == 1) or (("2016H" not in fname) and ((triggerBits >> 3) & 1) == 1 ) and (HLTpt[3]>=35)): #just run once more!!! I modified this for HLT35GeV but did not produce new results!!!
			hltPathTriggered_OS[numberOfHLTTriggers+3][0] = 1
		else:
	       		hltPathTriggered_OS[numberOfHLTTriggers+3][0] = 0
       	               #bitIndex==4:  # should be 4 for 2016H and 3 for 2016G
       	elif MC:
		if (((((triggerBits >> 4) & 1) == 1) or (((triggerBits >> 3) & 1) == 1))  and (HLTpt[3]>=35)):
			hltPathTriggered_OS[numberOfHLTTriggers+3][0] = 1
		else:
			hltPathTriggered_OS[numberOfHLTTriggers+3][0] = 0
	 
        if not "Run" in fname:
		puweight = tIn.puweight
	else:
		puweight = 1
        bkgSubW[0] = 1. if tIn.isOS else -1.


	if(tIn.byVVLooseIsolationMVArun2017v2DBoldDMwLT2017 > 0.5):
        #Mass cuts, mt and mvis for DY Jets
	    if DYJets:
		    if(tIn.mT < 30 and tIn.mVis >40 and tIn.mVis < 80):
			    if(tauGenMatching):  #for tau gen matching
				    if(tIn.tau_genindex > 0):
					    tOut.Fill()
			    else:
				    tOut.Fill()
        #High mT requirement for WJets
	    elif not DYJets:
		    if(tIn.mT > 30):
			    if(tauGenMatching):    #for tau gen matching
				    if(tIn.tau_genindex > 0):
					    tOut.Fill()
			    else:
				    tOut.Fill()	
			
tOutNames.Write()
tOut.Write()
fOut.Close()


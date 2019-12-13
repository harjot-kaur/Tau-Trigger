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
	excludeLumiSections = True
	print "==> SS events are stored as weights and applied to suppress the bkg for Data samples! <=="
	print "==> The additional lumi sections in \"Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt\" compared to \"ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt\" are removed for Data samples! <=="

disabledPScolumns = False  # True; to remove the disabled columns, False; to consider all columns

# the hadd of all the output ntuples
path = "/eos/user/h/hsert/TriggerStudies/ForkedRepo/Samples/12062018/"
#Final files for the second round of SFs
if MC:
	if DYJets:
		fname =  path + "NTuple_DYJetsToLL_12Apr2018_v1Andext1v1_12062018_PU_1000binMC.root"
	if not DYJets:
		fname =  path + "NTuple_0WJets_12Apr2018_12062018_PU_1000binMC.root"
else:
	fname =  path + "NTuple_Data_Run2017BCDEF_31Mar2018_12062018.root"

#Final files for the first round of SFs
#path = "/afs/cern.ch/user/h/hsert/TriggerStudies/ForkedRepo/Samples/2018_01_14/"
#fname =  path + "NTuple_Data2017BCDEF_17Nov2017-v1_14_01_2018.root"
#fname =  path +"NTuple_DYJets_RunIIFall17MiniAOD-RECOSIMstep_94X_mc2017_realistic_nomPlusExt_14_01_2018_PU.root"

#pt = [20, 26, 30, 34]
pt = [20, 26, 30, 32, 34]
numberOfHLTTriggers = 23
numberOfHLTTriggersForFit = 26


#######################################################
fIn = TFile.Open(fname)
tIn = fIn.Get('Ntuplizer/TagAndProbe')
tTriggerNames = fIn.Get("Ntuplizer/triggerNames")
if MC:
	suppressionType = "OStauGenMatched"
else:
	suppressionType = "SSsubtraction"
if DYJets:
	outname = fname.replace ('.root', '_' + suppressionType + 'TightWPold2017v1_forFit2.root')
else:
	outname = fname.replace ('.root', '_' + suppressionType + '_WjetEnriched_MediumWP2017v2_forFit.root')
fOut = TFile (outname, 'recreate')
tOut = tIn.CloneTree(0)
tOutNames = tTriggerNames.CloneTree(-1) # copy all

briso   = [n.zeros(1, dtype=int) for x in range (0, len(pt))]
brnoiso = [n.zeros(1, dtype=int) for x in range (0, len(pt))]
bkgSubW = n.zeros(1, dtype=float)
bkgSubANDpuW = n.zeros(1, dtype=float)

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
		name = ("hasHLTetauPath_13")# _IsoMu20_LooseChargedIsoPFTau27_plusL1Tau26andHLTTau30")
	elif(i==numberOfHLTTriggers+1):
		name = ("hasHLTmutauPath_13") #_IsoMu20_LooseChargedIsoPFTau27_plusL1Tau32")
	elif(i==numberOfHLTTriggers+2):
		name = ("hasHLTditauPath_9or10or11")#_TightTau35TightIDor_MediumTau35TightIDplusHLTTau40or_TightTau35plusHLTTau40_plusL1Tau32")
	tOut.Branch(name, hltPathTriggered_OS[i], name+"/I")

#tOut.Branch("isoHLT", hltPathTriggered_OS[6], name+"/I")

tOut.Branch("bkgSubW", bkgSubW, "bkgSubW/D")
tOut.Branch("bkgSubANDpuW", bkgSubANDpuW, "bkgSubANDpuW/D")

nentries = tIn.GetEntries()
for ev in range (0, nentries):
    tIn.GetEntry(ev)
    if (ev%10000 == 0) : print ev, "/", nentries

    if abs(tIn.tauEta) > 2.1:
        continue

    if saveOnlyOS and not tIn.isOS:
        continue

    for i in range (0, len(pt)):
        briso[i][0] = 0
        brnoiso[i][0] = 0

    for i in range (0, numberOfHLTTriggersForFit):
        hltPathTriggered_OS[i][0] = 0

    L1iso = True if tIn.l1tIso == 1 else False
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
    RunNumber = tIn.RunNumber
    lumi = tIn.lumi

    for bitIndex in range(0, numberOfHLTTriggers):
	import itertools as it
    	if bitIndex in it.chain(range(6, 13), range(19, 23)):   # apply this L1 cut only for di-tau triggers
            if (bitIndex==9 or bitIndex==10):
                if ((triggerBits >> bitIndex) & 1) == 1 and (L1pt>=32) and HLTpt[bitIndex]>=40:
                    hltPathTriggered_OS[bitIndex][0] = 1
                else:
                    hltPathTriggered_OS[bitIndex][0] = 0
            else:
                if ((triggerBits >> bitIndex) & 1) == 1 and (L1pt>=32):
                    hltPathTriggered_OS[bitIndex][0] = 1
                else:
                    hltPathTriggered_OS[bitIndex][0] = 0
        else:
            if ((triggerBits >> bitIndex) & 1) == 1:
            	hltPathTriggered_OS[bitIndex][0] = 1
            else:
                hltPathTriggered_OS[bitIndex][0] = 0
                        
        if(bitIndex==13):	
            if (((triggerBits >> bitIndex) & 1) == 1 and L1pt>=26 and L1iso>0 and HLTpt[bitIndex]>=30):
                hltPathTriggered_OS[numberOfHLTTriggers][0] = 1	  # this is the path for etau trigger. So (L1iso) should be applied here!
            else:
                hltPathTriggered_OS[numberOfHLTTriggers][0] = 0

            if ((triggerBits >> bitIndex) & 1) == 1:
                hltPathTriggered_OS[numberOfHLTTriggers+1][0] = 1 # this is the path for mutau trigger. So no extra requirement is needed like: L1pt and L1iso and HLTpt
            else:
                hltPathTriggered_OS[numberOfHLTTriggers+1][0] = 0       		
        
    if (((((triggerBits >> 9) & 1) == 1 and HLTpt[9]>=40) or (((triggerBits >> 10) & 1) == 1 and HLTpt[10]>=40) or (((triggerBits >> 11) & 1) == 1))  and L1pt>=32):
        hltPathTriggered_OS[numberOfHLTTriggers+2][0] = 1  # this is the path for di-tau trigger. HLTpt cut is required to have the same threshold on tau + L1Pt is needed due to L1 matching differences between MC and Data
    else:
        hltPathTriggered_OS[numberOfHLTTriggers+2][0] = 0
		
    bkgSubW[0] = 1. if tIn.isOS else -1.
    #if (L1pt > 26) and (L1iso) and (HLTpt > 32) and (((triggerBits >> 2) & 1) == 1):
    #    hltPathTriggered_OS[6][0] = 1
    #else:
    #    hltPathTriggered_OS[6][0] = 0

    if not "Data" in fname:
        puweight = tIn.puweight
    else:
        puweight = 1

    bkgSubANDpuW[0] = bkgSubW[0]*puweight

    if(tIn.byTightIsolationMVArun2017v1DBoldDMwLT2017 > 0.5):
        #Mass cuts, mt and mvis for DY Jets
	    if DYJets:
		    if(tIn.mT < 30 and tIn.mVis >40 and tIn.mVis < 80):
			    if(tauGenMatching):  #for tau gen matching
				    if(tIn.tau_genindex > 0):
					    tOut.Fill()
			    elif(excludeLumiSections):   # for removing the extra lumi sections in Final golden JSON file compared to Rereco ones
				    if not ((RunNumber==299480 and lumi==7) or (RunNumber==301397 and lumi==518)):
					    tOut.Fill()
			    else:
				    tOut.Fill()
        #High mT requirement for WJets
	    elif not DYJets:
		    if(tIn.mT > 30):
			    if(tauGenMatching):    #for tau gen matching
				    if(tIn.tau_genindex > 0):
					    tOut.Fill()
			    elif(excludeLumiSections):    # for removing the extra lumi sections in Final golden JSON file compared to Rereco ones
				    if not ((RunNumber==299480 and lumi==7) or (RunNumber==301397 and lumi==518)):
					    tOut.Fill()
			    else:
				    tOut.Fill()
"""
if(disabledPScolumns):   # for removing the disabled PS columns:
    if((RunNumber<305177 and PS_column>=2) or (RunNumber>=305178 and RunNumber<=305387 and PS_column>=2 and PS_column!=10) or (RunNumber>=305388 and PS_column>=3 and PS_column!=11 and PS_column!=12)):
	    tOut.Fill()
"""

tOutNames.Write()
tOut.Write()
fOut.Close()

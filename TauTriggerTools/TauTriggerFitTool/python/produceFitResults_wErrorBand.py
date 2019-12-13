from ROOT import *
from array import array
gROOT.SetBatch(True)
from math import sqrt


files2017 = ("/afs/cern.ch/user/h/hsert/public/Fall17Samples_31MarData_12AprMC/NTuple_Data_Run2017BCDEF_31Mar2018_SSsubtraction_VVLooseWP2017v2.root", "/afs/cern.ch/user/h/hsert/public/Fall17Samples_31MarData_12AprMC/NTuple_DYJetsToLL_12Apr2018_v1Andext1v1_12062018_puWeightsANDtauEScorrectionIncluded_OStauGenMatched_VVLooseWP2017v2.root")

#files = ("/eos/user/h/hsert/TriggerStudies/ForkedRepo/Samples/12062018/NTuple_Data_Run2017BCDEF_31Mar2018_12062018_SSsubtraction_WjetEnriched_VVLooseWP2017v2__forFit.root", "/eos/user/h/hsert/TriggerStudies/ForkedRepo/Samples/12062018/NTuple_0WJets_12Apr2018_12062018_PU_1000binMC_OStauGenMatched_WjetEnriched_VVLooseWP2017v2__forFit.root")

files2018 = "/eos/user/h/hsert/TriggerStudies/ForkedRepo/Samples2018/190110/Ntuple_SingleMuon_Run2018ABCDReReco17SepPromptRecoD_190121_SSsubtraction_MediumWP2017v2_forFit.root", "/eos/user/h/hsert/TriggerStudies/ForkedRepo/Samples2018/190110/Ntuple_SingleMuon_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_190121_PUreweight1000MCbin_OStauGenMatched_MediumWP2017v2_forFit.root"

files = files2017

gStyle.SetFrameLineWidth(1)
gStyle.SetPadBottomMargin(0.13)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)

triggers = ["ditau", "mutau", "etau"]
types = ["DATA", "MC"]
WPs = ["vvlooseTauMVA", "vlooseTauMVA", "looseTauMVA", "mediumTauMVA", "tightTauMVA", "vtightTauMVA", "vvtightTauMVA"]
tauDMs = ["dm0", "dm1", "dm10"]


isDMspesific = True
if "Samples2018" in files:
	outputname =  "tauTriggerEfficiencies2018_firsttrial.root"


outputname = "tauTriggerEfficiencies2017_final_etauESshift.root"
if(isDMspesific): outputname = "tauTriggerEfficiencies2017_final_etauESshift_forDMspesific.root"

hPtDen = [[],[],[]]
hPtNum = [[],[],[]]

hPtDenDM = [[],[],[],[]]
hPtNumDM = [[],[],[],[]]

for ipath, trigger in enumerate(triggers):
	
	if (ipath == 0 ): #ditau
		edges =[]
		for i in range( 20, 50, 2) :
			edges.append( float(i) )
		for i in range( 50, 75, 5 ) :
			edges.append( float(i) )
		for i in range( 75, 100, 25) :
   			edges.append( float(i) )
		for i in range(100, 150, 50 ) :
			edges.append( float(i) )
		for i in range(150, 500, 300 ) :
 	  	 	edges.append( float(i) )
	if (ipath == 1 ): #mutau
		edges =[]
		for i in range( 20, 50, 2) :
			edges.append( float(i) )
		for i in range( 50, 75, 5 ) :
			edges.append( float(i) )
		for i in range( 75, 100, 25) :
   			edges.append( float(i) )
		for i in range(100, 150, 50 ) :
			edges.append( float(i) )
		for i in range(150, 500, 300 ) :
 	  	 	edges.append( float(i) )
	if (ipath == 2): #etau
		edges =[]
		for i in range( 20, 24, 1) :
			edges.append( float(i) )
		for i in range( 24, 30, 3) :
			edges.append( float(i) )
		for i in range( 30, 55, 5 ) :
			edges.append( float(i) )
		for i in range( 55, 70, 15 ) :
			edges.append( float(i) )
		for i in range( 70, 90, 20) :
			edges.append( float(i) )
		for i in range(90, 160, 50 ) :
			edges.append( float(i) )
		for i in range(160, 600, 350 ) :
			edges.append( float(i) )

	print trigger,"trigger", "binning: N=",  len(edges), edges
	
	hPtNum.append([])
	hPtDen.append([])
	# per DM
	hPtNumDM.append([])
	hPtDenDM.append([])

	for index, typ in enumerate(types):
		hPtNum[ipath].append([])
		hPtDen[ipath].append([])
		# per DM
		hPtNumDM[ipath].append([])
		hPtDenDM[ipath].append([])

		for ind, wp in enumerate(WPs):
			histoname = "histo_" + typ + "_" + wp + "_" + trigger
			hPtNumDM[ipath][index].append([])
			hPtDenDM[ipath][index].append([])
			
			hPtDen[ipath][index].append(TH1F (histoname + "_Den", "", len(edges)-1, array('f',edges)))
			hPtNum[ipath][index].append(TH1F (histoname + "_Num", "", len(edges)-1, array('f',edges)))

			# per DM
			for idm, DM in enumerate(tauDMs):
				if (ipath == 0 ): #ditau
					edgesDM =[]
					if (DM =="dm0"):
						for i in range( 20, 30, 2) :
							edgesDM.append( float(i) )
						for i in range( 30, 55, 5) :
							edgesDM.append( float(i) )
						for i in range( 55, 75, 20 ) :
							edgesDM.append( float(i) )
						for i in range( 75, 112, 37) :
							edgesDM.append( float(i) )
						for i in range(112, 500, 300 ) :
							edgesDM.append( float(i) )
					elif(DM =="dm1"):
						for i in range( 20, 30, 2) :
							edgesDM.append( float(i) )
						for i in range( 30, 55, 5) :
							edgesDM.append( float(i) )
						for i in range( 55, 75, 10 ) :
							edgesDM.append( float(i) )
						for i in range( 75, 90, 25) :
							edgesDM.append( float(i) )
						for i in range(90, 112, 22 ) :
							edgesDM.append( float(i) )
						for i in range(112, 500, 300 ) :
							edgesDM.append( float(i) )
					elif(DM=="dm10"):
						if ( typ == "DATA"):
							for i in range(20, 30, 2) :
								edgesDM.append( float(i) )
							for i in range(30, 50, 5) :
								edgesDM.append( float(i) )
							for i in range(50, 70, 10 ) :
								edgesDM.append( float(i) )
							for i in range(70, 85, 15 ) :
								edgesDM.append( float(i) )
							for i in range(85, 115, 30 ) :
								edgesDM.append( float(i) )
							for i in range(115, 500, 300 ) :
								edgesDM.append( float(i) )
						elif (typ == "MC"):
							for i in range(20, 30, 2) :
								edgesDM.append( float(i) )
							for i in range(30, 50, 5) :
								edgesDM.append( float(i) )
							for i in range(50, 58, 8 ) :
								edgesDM.append( float(i) )
							for i in range(58, 78, 20 ) :
								edgesDM.append( float(i) )
							for i in range(78, 120, 42 ) :
								edgesDM.append( float(i) )
							for i in range(120, 500, 300 ) :
								edgesDM.append( float(i) )
				if (ipath == 1 ): #mutau
					edgesDM =[]
					if (DM =="dm10"):
						#if(typ == "MC"):
						for i in range(20, 24, 1) :
							edgesDM.append( float(i) )
						for i in range(24, 30, 2) :
							edgesDM.append( float(i) )
						for i in range(30, 40, 5 ) :
							edgesDM.append( float(i) )
						for i in range(40, 58, 9 ) :
							edgesDM.append( float(i) )
						for i in range(58, 80, 22 ) :
							edgesDM.append( float(i) )
						for i in range(80, 115, 35) :
							edgesDM.append( float(i) )
						for i in range(115,  500, 300 ) :
							edgesDM.append( float(i) )
					elif(DM =="dm1"):
						for i in range( 20, 30, 2) :
							edgesDM.append( float(i) )
						for i in range( 30, 55, 5) :
							edgesDM.append( float(i) )
						for i in range( 55, 65, 10 ) :
							edgesDM.append( float(i) )
						for i in range( 65, 80, 15) :
							edgesDM.append( float(i) )
						for i in range(80, 110, 30 ) :
							edgesDM.append( float(i) )
						for i in range(110, 500, 300 ) :
							edgesDM.append( float(i) )
					else:
						for i in range( 20, 30, 2) :
							edgesDM.append( float(i) )
						for i in range( 30, 55, 5) :
							edgesDM.append( float(i) )
						for i in range( 55, 65, 10 ) :
							edgesDM.append( float(i) )
						for i in range( 65, 80, 15) :
							edgesDM.append( float(i) )
						for i in range(80, 110, 30 ) :
							edgesDM.append( float(i) )
						for i in range(110, 500, 300 ) :
							edgesDM.append( float(i) )
				if (ipath == 2 ): #etau
					edgesDM =[]
					if (DM =="dm10"):
						for i in range( 20, 30, 2) :
							edgesDM.append( float(i) )
						for i in range( 30, 55, 5) :
							edgesDM.append( float(i) )
						for i in range( 55, 75, 10 ) :
							edgesDM.append( float(i) )
						for i in range( 75, 100, 25) :
							edgesDM.append( float(i) )
						for i in range(100, 125, 25 ) :
							edgesDM.append( float(i) )
						for i in range(125, 500, 325 ) :
							edgesDM.append( float(i) )
					elif (DM =="dm1"):
						for i in range( 20, 30, 2) :
							edgesDM.append( float(i) )
						for i in range( 30, 40, 5) :
							edgesDM.append( float(i) )
						for i in range( 40, 55, 15 ) :
							edgesDM.append( float(i) )
						for i in range( 55, 75, 20) :
							edgesDM.append( float(i) )
						for i in range(75, 125, 50 ) :
							edgesDM.append( float(i) )
						for i in range(125, 500, 350 ) :
							edgesDM.append( float(i) )
					elif (DM =="dm0"):
						if(typ == "DATA"):
							for i in range( 20, 22, 1) :
								edgesDM.append( float(i) )
							for i in range( 22, 30, 2) :
								edgesDM.append( float(i) )
							for i in range( 30, 40, 5) :
								edgesDM.append( float(i) )
							for i in range( 40, 60, 10 ) :
								edgesDM.append( float(i) )
							for i in range( 60, 75, 15) :
								edgesDM.append( float(i) )
							for i in range(75, 95, 20) :
								edgesDM.append( float(i) )
							for i in range(95, 500, 300 ) :
								edgesDM.append( float(i) )
						elif(typ == "MC"):
							for i in range( 20, 22, 1) :
								edgesDM.append( float(i) )
							for i in range( 22, 30, 2) :
								edgesDM.append( float(i) )
							for i in range( 30, 50, 5) :
								edgesDM.append( float(i) )
							for i in range( 50, 60, 10 ) :
								edgesDM.append( float(i) )
							for i in range( 60, 80, 20) :
								edgesDM.append( float(i) )
							for i in range(80, 100, 20 ) :
								edgesDM.append( float(i) )
							for i in range(100, 500, 350 ) :
								edgesDM.append( float(i) )

				hPtDenDM[ipath][index][ind].append(TH1F (histoname + "_Den_" + DM, "", len(edgesDM)-1, array('f',edgesDM)))
				hPtNumDM[ipath][index][ind].append(TH1F (histoname + "_Num_" + DM, "", len(edgesDM)-1, array('f',edgesDM)))


# This function is taken from Tyler Ruggles SF tool
# https://github.com/truggles/TauTriggerSFs2017/blob/master/python/helpers.py#L9-L40
# Function to create TH1Fs from TGraphAsymmErrors
# This does not preserve the asymmetric errors, only
# bin width and value and does a rough approximation
# on symmetric errors.
def getTH1FfromTGraphAsymmErrors( asym, name ) :

    # Holding vals for TH1F binning and y-vals
    xSpacing = array( 'd', [] )
    yVals = array( 'd', [] )
    yErrors = array( 'd', [] )

    nVals = asym.GetN()
    x, y = Double(0.), Double(0.)
    xEPlus, xEMin = 0., 0.
    yEPlus, yEMin = 0., 0.

    for n in range( nVals ) :
        asym.GetPoint( n, x, y )
        xEPlus = asym.GetErrorXhigh( n )
        xEMin = asym.GetErrorXlow( n )
        yEPlus = asym.GetErrorYhigh( n )
        yEMin = asym.GetErrorYlow( n )
        xSpacing.append( x-xEMin )
        yVals.append( y )
        # To simplify, take asymm errors and go to approximation
        # of symmetric for TH1
        yErrors.append( sqrt(yEPlus**2 + yEMin**2) )

    # Don't forget to add the high end of last bin
    xSpacing.append( x+xEPlus )
    
    outH = TH1F( name, name, len(xSpacing)-1, xSpacing )
    for bin in range( 1, outH.GetNbinsX()+1 ) :
        outH.SetBinContent( bin, yVals[bin-1] )
        outH.SetBinError( bin, yErrors[bin-1] )
    return outH
    
def getScaleFactor(func):

	SF = TGraphAsymmErrors()
	for i in range(20, 450):
		SF.SetPoint(i, i, (func[0].Eval(i)/func[1].Eval(i)))
		
	SF.Draw("A*")
	SF.GetXaxis().SetLimits(18,600)
	SF.GetYaxis().SetRangeUser(0,1.5)
	SF.GetXaxis().SetMoreLogLabels()
	SF.SetMarkerStyle(20)
	SF.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
	SF.GetYaxis().SetTitle("SF: Data/MC")
	return SF

def getScaleFactorWithPropagatedError(func, histo):

	SF = TGraphErrors()
	SF2 = TGraphAsymmErrors()

	hEffi1 = TH1F("", "", 480, 20, 500 )
	hEffi2 = TH1F("", "", 480, 20, 500 )

	for i in range(0, histo[0].GetNbinsX()):

		sigmaA = histo[0].GetBinError(i)
		sigmaB = histo[1].GetBinError(i)
		funcA = func[0].Eval(histo[0].GetBinCenter(i))
		funcB = func[1].Eval(histo[1].GetBinCenter(i))

		SF.SetPoint( i, histo[0].GetBin(i), (func[0].Eval(histo[0].GetBin(i))/func[1].Eval(histo[0].GetBin(i))))
		SF.SetPointError(i, 0,  sqrt(pow(funcA/funcB,2)*(pow(sigmaA/funcA ,2) + pow(sigmaB/funcB,2))))

		hEffi1.SetBinContent(i, histo[0].GetBinContent(i))
		hEffi1.SetBinError(i, histo[0].GetBinError(i))
		hEffi2.SetBinContent(i, histo[1].GetBinContent(i))
		hEffi2.SetBinError(i, histo[1].GetBinError(i))

	SF.Draw("A*")
	SF.GetXaxis().SetLimits(18,600)
	SF.GetXaxis().SetMoreLogLabels()
	SF.GetYaxis().SetRangeUser(0,1.5)
	SF.SetMarkerStyle(20)
	SF.SetLineColor(kOrange)
	SF.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
	SF.GetYaxis().SetTitle("SF: Data/MC")
	SF2.GetXaxis().SetLimits(18,600)
	SF2.GetXaxis().SetMoreLogLabels()
	SF2.GetYaxis().SetRangeUser(0,1.5)
	SF2.SetMarkerStyle(20)
	SF2.Draw("AP")

	hEffi1.Divide(hEffi2)
	hEffi1.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
	hEffi1.GetYaxis().SetTitle("SF: Data/MC")
	hEffi1.GetXaxis().SetRangeUser(18,600)
	hEffi1.GetXaxis().SetMoreLogLabels()
	hEffi1.GetYaxis().SetRangeUser(0,1.5)
	hEffi1.SetMarkerStyle(20)
	hEffi1.SetLineWidth(6)

	return SF, hEffi1

def getScaleFactorWithShiftedError(func, histo):

	SF1 = TGraphAsymmErrors()
	SF2 = TGraphAsymmErrors()
	SF3 = TGraphAsymmErrors()
	SF4 = TGraphAsymmErrors()
	SF5 = TGraphAsymmErrors()
	SFnew = TGraphAsymmErrors()
	for i in range(0, histo[0].GetNbinsX()):

		funcA = func[0].Eval(histo[0].GetBinCenter(i))
		funcB = func[1].Eval(histo[0].GetBinCenter(i))
		sigmaA = histo[0].GetBinError(i)
		sigmaB = histo[1].GetBinError(i)

		shift0 = funcA / funcB
		shiftPP = (funcA + sigmaA) / (funcB + sigmaB)
		shiftPM = (funcA + sigmaA) / (funcB - sigmaB)
		shiftMP = (funcA - sigmaA) / (funcB + sigmaB)
		shiftMM = (funcA - sigmaA) / (funcB - sigmaB)

		shiftPNom = (funcA + sigmaA) / (funcB)
		shiftMNom = (funcA - sigmaA) / (funcB)
		shiftNomP = (funcA) / (funcB + sigmaB)
		shiftNomM = (funcA) / (funcB - sigmaB)

		SF1.SetPoint( i, histo[0].GetBinCenter(i), shift0)
		SF2.SetPoint( i, histo[0].GetBinCenter(i), shiftPNom)
		SF3.SetPoint( i, histo[0].GetBinCenter(i), shiftMNom)
		SF4.SetPoint( i, histo[0].GetBinCenter(i), shiftNomP)
		SF5.SetPoint( i, histo[0].GetBinCenter(i), shiftNomM)

		errorP = sqrt(pow((shiftPNom - shift0),2) + pow((shiftNomM - shift0),2))
		errorM = sqrt(pow((shift0 - shiftNomP),2) + pow((shift0 - shiftNomM),2))

		SFnew.SetPoint( i, histo[0].GetBinCenter(i), shift0)
		SFnew.SetPointError( i, 0, 0, errorP, errorM)

	SF1.SetMarkerStyle(20)
	SF1.SetLineColor(kBlue)
	SF2.SetMarkerStyle(20)
	SF3.SetMarkerStyle(20)
	SF4.SetMarkerStyle(20)
	SF5.SetMarkerStyle(20)

	SF1.SetMarkerColor(kBlack)
	SF2.SetMarkerColor(kRed)
	SF3.SetMarkerColor(kBlue)
	SF4.SetMarkerColor(kGreen)
	SF5.SetMarkerColor(kViolet)

	SF1.SetMarkerSize(0.5)
	SF2.SetMarkerSize(0.5)
	SF3.SetMarkerSize(0.5)
	SF4.SetMarkerSize(0.5)
	SF5.SetMarkerSize(0.5)

	SF2.SetLineColor(kOrange)

	SFnew.SetMarkerStyle(20)
	SFnew.SetMarkerColor(kBlack)
	SFnew.SetLineColor(kBlue)

	mg = TMultiGraph()
	mg.Add(SF1)
	mg.Add(SF2)
	mg.Add(SF3)
	mg.Add(SF4)
	mg.Add(SF5)
	mg.Draw("AP")
	mg.GetXaxis().SetLimits(18,600)
	mg.GetYaxis().SetRangeUser(0,1.5)
	mg.GetXaxis().SetMoreLogLabels()

	mg.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
	mg.GetYaxis().SetTitle("SF: Data/MC")
	return SFnew

def getConfidenceInterval(g_efficiency, histo, graph, CL): # give CL as 0.68 or 0.95

	TVirtualFitter.GetFitter().GetConfidenceIntervals(histo, CL)

	for i in range(0, g_efficiency.GetN()):
		graph.SetPoint(i, g_efficiency.GetX()[i], 0)
		TVirtualFitter.GetFitter().GetConfidenceIntervals(graph, CL)

	return histo, graph


def createRelativeErrors(fitresult, g_efficiency, f1):

	# confidence interval#
	values = fitresult.GetConfidenceIntervals(0.68, False)
	interval = TGraphErrors(g_efficiency.GetN())
	ratio = TGraphAsymmErrors(g_efficiency.GetN())
	ratio2 = TGraphAsymmErrors(g_efficiency.GetN())
	for i in range(0, g_efficiency.GetN()):
		interval.SetPoint(i, g_efficiency.GetX()[i], f1[index].Eval(g_efficiency.GetX()[i] ))
		interval.SetPointError(i, 0, values[i] )
		ratio.SetPoint(i, g_efficiency.GetX()[i], (f1[index].Eval(g_efficiency.GetX()[i]) - values[i])/f1[index].Eval(g_efficiency.GetX()[i]))
		ratio2.SetPoint(i, g_efficiency.GetX()[i], (f1[index].Eval(g_efficiency.GetX()[i]) + values[i])/f1[index].Eval(g_efficiency.GetX()[i]))

	ratio.GetXaxis().SetTitleSize(0.05)
	ratio.GetYaxis().SetTitleSize(0.05)
	ratio.GetXaxis().SetTitleOffset(1.1)
	ratio.GetYaxis().SetTitleOffset(1.1)
	ratio.SetLineWidth(2)
	ratio2.SetLineWidth(2)
	ratio.Draw("A*")
	ratio2.Draw("*same")
	ratio.SetMarkerStyle(20)
	ratio2.SetMarkerStyle(20)
	ratio.GetYaxis().SetRangeUser(0.8,1.2)
	ratio.SetTitle("")
	ratio.SetTitle(trigger +"Path_" + wp +"_"+ typ)
	ratio.GetYaxis().SetTitle("#frac{fit #pm error}{fit}")
	ratio.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
	ratio2.SetMarkerColor(2)
	ratio.GetXaxis().SetMoreLogLabels()

	return ratio , ratio2

 
for index, filename in enumerate(files):
	print  "filename", filename
	
	file = TFile.Open(filename)
	tree = file.Get('TagAndProbe')
	triggerNamesTree = file.Get("triggerNames")	
	
	print "Populating histograms"

	Nevts = 0
	for iEv in range (0, tree.GetEntries()):
		tree.GetEntry(iEv)
			
		tauPt = tree.tauPt
		HLTPt = tree.hltPt
		tauEta = tree.tauEta
		tauPhi = tree.tauPhi
		tauDM = tree.tauDM
		isOS = tree.isOS
		Nvtx = tree.Nvtx
		
		if ("MC" in filename and "DYJets" in filename):
			puweight = tree.puweight
			tauPt_ESshifted = tree.tauPt_ESshifted
		else:
			puweight=1
			tauPt_ESshifted = tree.tauPt

		vvlooseWP = tree.byVVLooseIsolationMVArun2017v2DBoldDMwLT2017	
		vlooseWP = tree.byVLooseIsolationMVArun2017v2DBoldDMwLT2017	
		looseWP = tree.byLooseIsolationMVArun2017v2DBoldDMwLT2017	
		mediumWP = tree.byMediumIsolationMVArun2017v2DBoldDMwLT2017	
		tightWP = tree.byTightIsolationMVArun2017v2DBoldDMwLT2017
		vtightWP = tree.byVTightIsolationMVArun2017v2DBoldDMwLT2017	
		vvtightWP = tree.byVVTightIsolationMVArun2017v2DBoldDMwLT2017	
		
		if("Run2017B" in filename or "12062018" in filename):
			hasHLTditauPath_9or10or11 = tree.hasHLTditauPath_9or10or11
			hasHLTetauPath_13 = tree.hasHLTetau_Path_13
			hasHLTmutauPath_13 = tree.hasHLTmutauPath_13
		else:
			hasHLTditauPath_4or5or6noHPS = tree.hasHLTditauPath_4or5or6noHPS
			hasHLTetauPath_8noHPS = tree.hasHLTetauPath_8noHPS
			hasHLTmutauPath_8noHPS = tree.hasHLTmutauPath_8noHPS
			hasHLTditauPath_20HPS = tree.hasHLTditauPath_20HPS
			hasHLTditauPath_15HPS = tree.hasHLTPath_15
			hasHLTetauPath_14HPS = tree.hasHLTetauPath_14HPS
			hasHLTmutauPath_14HPS = tree.hasHLTmutauPath_14HPS

		Nevents = tree.EventNumber
		Nevts =Nevts + 1
		
		#bkgSubW = 1. if tree.isOS else -1.
		weight = tree.bkgSubW*puweight
		
		if("Run2018A" in filename or "Autumn18" in filename):
			#print "2018 hlt paths are considering now..."
			HLTHPSpaths18 = [hasHLTditauPath_20HPS, hasHLTmutauPath_14HPS , hasHLTetauPath_14HPS ]
			HLTHPSpaths18MC = [hasHLTditauPath_15HPS, hasHLTmutauPath_14HPS , hasHLTetauPath_14HPS ]
			HLTpaths18 = [hasHLTditauPath_4or5or6noHPS, hasHLTmutauPath_8noHPS, hasHLTetauPath_8noHPS]
		else:
			HLTpaths17 = [hasHLTditauPath_9or10or11, hasHLTmutauPath_13 , hasHLTetauPath_13 ]

		WPoints = [vvlooseWP, vlooseWP, looseWP, mediumWP, tightWP, vtightWP, vvtightWP]
		
		oneProng = False
		oneProngPiZero= False
		threeProng = False

		if (tauDM ==0):
			oneProng = True
		if (tauDM ==1):
			oneProngPiZero = True
		if (tauDM ==10):
			threeProng = True
		DMs = [oneProng, oneProngPiZero, threeProng]

		if(tauEta > 0 and tauEta < 1.479 and tauPhi > 2.8):
			PixelBarrel = True

		elif((tauEta > 0 and tauEta < 1.479  and tauPhi < 2.8) or (tauEta <= 0 and tauEta > -1.479)):
			BarrelRest = True
		elif(abs(tauEta) > 1.479):
			EndCap = True
		else:
			print "is there other eta=phi regions?"

		# Filling the histograms
		for WPind, WP in enumerate(WPoints):
			if((WP > 0) and ("Run2017B" in filename or "12062018" in filename)):
				for ipath, trigger in enumerate(HLTpaths17):
					hPtDen[ipath][index][WPind].Fill(tauPt_ESshifted, weight)
					for idm, DM in enumerate(DMs):
						if(DM == True):
							hPtDenDM[ipath][index][WPind][idm].Fill(tauPt_ESshifted, weight)
				if ( HLTpaths17[0] >0 ):
					hPtNum[0][index][WPind].Fill(tauPt_ESshifted, weight)
					for idm, DM in enumerate(DMs):
						if(DM == True):
							hPtNumDM[0][index][WPind][idm].Fill(tauPt_ESshifted, weight)
				if( HLTpaths17[1] >0 ):
					hPtNum[1][index][WPind].Fill(tauPt_ESshifted, weight)
					for idm, DM in enumerate(DMs):
						if(DM == True):
							hPtNumDM[1][index][WPind][idm].Fill(tauPt_ESshifted, weight)
				if( HLTpaths17[2] >0 ):
					hPtNum[2][index][WPind].Fill(tauPt_ESshifted, weight)
					for idm, DM in enumerate(DMs):
						if(DM == True):
							hPtNumDM[2][index][WPind][idm].Fill(tauPt_ESshifted, weight)

			elif((WP > 0) and ("Run2018A" in filename or "Autumn18" in filename)):
				if("Run2018A" in filename):
					if(tree.RunNumber < 317509):
						for ipath, trigger in enumerate(HLTpaths18):
							hPtDen[ipath][index][WPind].Fill(tauPt, weight)
							for idm, DM in enumerate(DMs):
								if(DM == True):
									hPtDenDM[ipath][index][WPind][idm].Fill(tauPt, weight)
						if ( HLTpaths18[0] >0 ):
							hPtNum[0][index][WPind].Fill(tauPt, weight)
							for idm, DM in enumerate(DMs):
								if(DM == True):
									hPtNumDM[0][index][WPind][idm].Fill(tauPt, weight)
						if( HLTpaths18[1] >0 ):
							hPtNum[1][index][WPind].Fill(tauPt, weight)
							for idm, DM in enumerate(DMs):
								if(DM == True):
									hPtNumDM[1][index][WPind][idm].Fill(tauPt, weight)
						if( HLTpaths18[2] >0 ):
							hPtNum[2][index][WPind].Fill(tauPt, weight)
							for idm, DM in enumerate(DMs):
								if(DM == True):
									hPtNumDM[2][index][WPind][idm].Fill(tauPt, weight)
					else:
						for ipath, trigger in enumerate(HLTHPSpaths18):
							hPtDen[ipath][index][WPind].Fill(tauPt, weight)
							for idm, DM in enumerate(DMs):
								if(DM == True):
									hPtDenDM[ipath][index][WPind][idm].Fill(tauPt, weight)
						if ( HLTHPSpaths18[0] >0 ):
							hPtNum[0][index][WPind].Fill(tauPt, weight)
							for idm, DM in enumerate(DMs):
								if(DM == True):
									hPtNumDM[0][index][WPind][idm].Fill(tauPt, weight)
						if( HLTHPSpaths18[1] >0 ):
							hPtNum[1][index][WPind].Fill(tauPt, weight)
							for idm, DM in enumerate(DMs):
								if(DM == True):
									hPtNumDM[1][index][WPind][idm].Fill(tauPt, weight)
						if( HLTHPSpaths18[2] >0 ):
							hPtNum[2][index][WPind].Fill(tauPt, weight)
							for idm, DM in enumerate(DMs):
								if(DM == True):
									hPtNumDM[2][index][WPind][idm].Fill(tauPt, weight)

				elif((WP > 0) and ("Autumn18" in filename)):

					for ipath, trigger in enumerate(HLTHPSpaths18MC):
						hPtDen[ipath][index][WPind].Fill(tauPt, weight)
						for idm, DM in enumerate(DMs):
								if(DM == True):
									hPtDenDM[ipath][index][WPind][idm].Fill(tauPt, weight)
					if ( HLTHPSpaths18MC[0] >0 ):
						hPtNum[0][index][WPind].Fill(tauPt, weight)
						for idm, DM in enumerate(DMs):
								if(DM == True):
									hPtNumDM[0][index][WPind][idm].Fill(tauPt, weight)
					if( HLTHPSpaths18MC[1] >0 ):
						hPtNum[1][index][WPind].Fill(tauPt, weight)
						for idm, DM in enumerate(DMs):
								if(DM == True):
									hPtNumDM[1][index][WPind][idm].Fill(tauPt, weight)
					if( HLTHPSpaths18MC[2] >0 ):
						hPtNum[2][index][WPind].Fill(tauPt, weight)
						for idm, DM in enumerate(DMs):
								if(DM == True):
									hPtNumDM[2][index][WPind][idm].Fill(tauPt, weight)



file = TFile( "../data/"+outputname, 'recreate')

# check bin contents of the histograms to avoid cases Num > Den, triggered due to SS subtraction (negative weights)
for ipath, trigger in enumerate(triggers):
	for WPind, wp in enumerate(WPs):
		for index, typ in enumerate(types):

			nbins = hPtNum[ipath][index][WPind].GetNbinsX()
			for binid in range(0, nbins + 1):
				if(hPtNum[ipath][index][WPind].GetBinContent(binid) > hPtDen[ipath][index][WPind].GetBinContent(binid)):
					hPtNum[ipath][index][WPind].SetBinContent(binid, hPtDen[ipath][index][WPind].GetBinContent(binid))

			for idm, DM in enumerate(tauDMs):
				for binid in range(0, nbins + 1):
					if(hPtNumDM[ipath][index][WPind][idm].GetBinContent(binid) > hPtDenDM[ipath][index][WPind][idm].GetBinContent(binid)):
						hPtNumDM[ipath][index][WPind][idm].SetBinContent(binid, hPtDenDM[ipath][index][WPind][idm].GetBinContent(binid))


# efficiency calculation after filling the histograms for 3 different triggers for each WPs of DATA and MC
for ipath, trigger in enumerate(triggers):
	g_errBand68 = []
	h_errBand68 = []
	h_errBandDM68 = [[], []]
	g_errBandDM68 = [[], []]
	for WPind, wp in enumerate(WPs):
		f1 =[]

		for idm, DM in enumerate(tauDMs):
			h_errBandDM68.append([])
			g_errBandDM68.append([])
			for index, typ in enumerate(types):
				h_errBandDM68[idm].append(TH1F(histoname+ "_" + DM +"_CL68","histo of 0.68 confidence band", 480, 20, 500))
				g_errBandDM68[idm].append(TGraphErrors())

		for index, typ in enumerate(types):
			f1.append(TF1( 'f1'+typ, '[5] - ROOT::Math::crystalball_cdf(-x, [0], [1], [2], [3])*([4])' ))
			if(index ==0):
				f1[index].SetLineColor( kBlue)
			else:
				f1[index].SetLineColor( kRed)
			f1[index].SetParName( 0, "alpha" )
			f1[index].SetParName( 1, "n" )
			f1[index].SetParName( 2, "simga" )
			f1[index].SetParName( 3, "x0" )
			f1[index].SetParName( 4, "scale" )
			f1[index].SetParName( 5, "y-rise" )

			f2 = [[],[]]
                for idm, DM in enumerate(tauDMs):
                        f2.append([])
                        for index, typ in enumerate(types):
                                f2[idm].append(TF1( 'f2_'+ DM  +"_" + typ, '[5] - ROOT::Math::crystalball_cdf(-x, [0], [1], [2], [3]\
)*([4])' ))
                                if(idm ==0): f2[idm][index].SetLineColor( kBlue )
                                elif(idm ==1): f2[idm][index].SetLineColor( kRed )
                                elif(idm ==2): f2[idm][index].SetLineColor( kGreen+3 )
                                if( isDMspesific):
                                        if index==0: f2[idm][0].SetLineColor( kBlue )
                                        elif index==1: f2[idm][1].SetLineColor( kRed )
                                f2[idm][index].SetParName( 0, "alpha" )
                                f2[idm][index].SetParName( 1, "n" )
                                f2[idm][index].SetParName( 2, "simga" )
                                f2[idm][index].SetParName( 3, "x0" )
                                f2[idm][index].SetParName( 4, "scale" )
                                f2[idm][index].SetParName( 5, "y-rise" )

		for index, typ in enumerate(types):
			# setting the fit parameters for various cases
			if(ipath == 0): # ditau
				f1[index].SetParameter( 0, 0.2)
				f1[index].SetParameter( 1, 5.0 )
				f1[index].SetParameter( 2, 7.0 )
				f1[index].SetParameter( 3, -30.)
				f1[index].SetParameter( 4, 1.0 )
				f1[index].SetParameter( 5, 1.0)

				for idm, DM in enumerate(tauDMs):
					if(idm ==0 or idm == 1):
						f2[idm][index].SetParameter( 0, 0.2)
						f2[idm][index].SetParameter( 1, 5.0 )
						f2[idm][index].SetParameter( 2, 7.0 )
						f2[idm][index].SetParameter( 3, -30.)
						f2[idm][index].SetParameter( 4, 1.0 )
						f2[idm][index].SetParameter( 5, 1.0)
					elif(idm ==2):
						if (wp == "tightTauMVA"):
							if(index == 0):
								f2[idm][index].SetParameter( 0, 2.0)
								f2[idm][index].SetParameter( 1, 10.0)
								f2[idm][index].SetParameter( 2, 0.3)
								f2[idm][index].SetParameter( 3, -30.)
								f2[idm][index].SetParameter( 4, 1.0 )
								f2[idm][index].SetParameter( 5, 1.0)
							elif(index == 1):
								f2[idm][index].SetParameter( 0, 2.0)
								f2[idm][index].SetParameter( 1, 10.0)
								f2[idm][index].SetParameter( 2, 0.5)
								f2[idm][index].SetParameter( 3, -30.)
								f2[idm][index].SetParameter( 4, 1.0 )
								f2[idm][index].SetParameter( 5, 1.0)
						elif (wp == "vtightTauMVA"):
                                                        if(index == 0):
                                                                f2[idm][index].SetParameter( 0, 2.0)
                                                                f2[idm][index].SetParameter( 1, 10.0)
                                                                f2[idm][index].SetParameter( 2, 0.3)
                                                                f2[idm][index].SetParameter( 3, -30.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
                                                        elif(index == 1):
                                                                f2[idm][index].SetParameter( 0, 2.0)
                                                                f2[idm][index].SetParameter( 1, 10.0)
                                                                f2[idm][index].SetParameter( 2, 0.4)
                                                                f2[idm][index].SetParameter( 3, -30.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)

						elif (wp == "vvtightTauMVA"):
							if(index==0):
								f2[idm][index].SetParameter( 0, 2.0)
								f2[idm][index].SetParameter( 1, 10.0)
								f2[idm][index].SetParameter( 2, 0.4)
								f2[idm][index].SetParameter( 3, -30.)
								f2[idm][index].SetParameter( 4, 1.0 )
								f2[idm][index].SetParameter( 5, 1.0 )
							elif(index == 1):
								f2[idm][index].SetParameter( 0, 2.0)
								f2[idm][index].SetParameter( 1, 10.0)
								f2[idm][index].SetParameter( 2, 0.5)
								f2[idm][index].SetParameter( 3, -30.)
								f2[idm][index].SetParameter( 4, 1.0 )
								f2[idm][index].SetParameter( 5, 1.0)
						elif(wp == "vlooseTauMVA"):
							if(index == 0):
								f2[idm][index].SetParameter( 0, 0.8)
								f2[idm][index].SetParameter( 1, 10.0)
								f2[idm][index].SetParameter( 2, 1.0)
								f2[idm][index].SetParameter( 3, -25.)
								f2[idm][index].SetParameter( 4, 1.0 )
								f2[idm][index].SetParameter( 5, 1.0)
							elif(index == 1):
								f2[idm][index].SetParameter( 0, 2.0)
								f2[idm][index].SetParameter( 1, 10.0)
								f2[idm][index].SetParameter( 2, 0.3)
								f2[idm][index].SetParameter( 3, -30.)
								f2[idm][index].SetParameter( 4, 1.0 )
								f2[idm][index].SetParameter( 5, 1.0)
						elif(wp== "mediumTauMVA"):
							if(index == 0):
                                                                f2[idm][index].SetParameter( 0, 0.8)
                                                                f2[idm][index].SetParameter( 1, 10.0)
                                                                f2[idm][index].SetParameter( 2, 1.0)
                                                                f2[idm][index].SetParameter( 3, -25.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
							elif(index == 1):
                                                                f2[idm][index].SetParameter( 0, 0.8)
                                                                f2[idm][index].SetParameter( 1, 10.0)
                                                                f2[idm][index].SetParameter( 2, 0.3)
                                                                f2[idm][index].SetParameter( 3, -30.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
						else:
                                                     if(index == 0):
                                                                f2[idm][index].SetParameter( 0, 0.8)
                                                                f2[idm][index].SetParameter( 1, 10.0)
                                                                f2[idm][index].SetParameter( 2, 1.0)
                                                                f2[idm][index].SetParameter( 3, -25.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
						     elif(index == 1):
                                                                f2[idm][index].SetParameter( 0, 2.0)
                                                                f2[idm][index].SetParameter( 1, 10.0)
                                                                f2[idm][index].SetParameter( 2, 2.0)
                                                                f2[idm][index].SetParameter( 3, -25.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
			if(ipath == 1): # mutau
				f1[index].SetParameter( 0, 0.2)
				f1[index].SetParameter( 1, 5.0 )
				f1[index].SetParameter( 2, 7.0 )
				f1[index].SetParameter( 3, -30.)
				f1[index].SetParameter( 4, 1.0 )
				f1[index].SetParameter( 5, 1.0)
				for idm, DM in enumerate(tauDMs):
					if(idm ==0 or idm == 1):
						f2[idm][index].SetParameter( 0, 0.5)
						f2[idm][index].SetParameter( 1, 5.0 )
						f2[idm][index].SetParameter( 2, 7.0 )
						f2[idm][index].SetParameter( 3, -20.)
						f2[idm][index].SetParameter( 4, 1.0 )
						f2[idm][index].SetParameter( 5, 1.0)
					elif(idm == 2): # dm10
						if(wp =="vvtightTauMVA"):
							if(typ == "MC"):
								f2[idm][index].SetParameter( 0, 0.8)
								f2[idm][index].SetParameter( 1, 10.0 )
								f2[idm][index].SetParameter( 2, 0.4 )
								f2[idm][index].SetParameter( 3, -30.)
								f2[idm][index].SetParameter( 4, 1.0 )
								f2[idm][index].SetParameter( 5, 1.0)
							elif(typ=="DATA"):
								f2[idm][index].SetParameter( 0, 0.8)
								f2[idm][index].SetParameter( 1, 10.0 )
								f2[idm][index].SetParameter( 2, 0.2 )
								f2[idm][index].SetParameter( 3, -30.)
								f2[idm][index].SetParameter( 4, 1.0 )
								f2[idm][index].SetParameter( 5, 1.0)
						elif(wp =="vtightTauMVA"):
                                                        if(typ == "MC"):
                                                                f2[idm][index].SetParameter( 0, 0.8)
                                                                f2[idm][index].SetParameter( 1, 10.0 )
                                                                f2[idm][index].SetParameter( 2, 0.4 )
                                                                f2[idm][index].SetParameter( 3, -30.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
                                                        elif(typ=="DATA"):
                                                                f2[idm][index].SetParameter( 0, 0.8)
                                                                f2[idm][index].SetParameter( 1, 10.0 )
                                                                f2[idm][index].SetParameter( 2, 0.2 )
                                                                f2[idm][index].SetParameter( 3, -25.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
								f2[idm][index].SetParameter( 5, 1.0)
						elif(wp =="tightTauMVA"):
                                                        if(typ == "MC"):
                                                                f2[idm][index].SetParameter( 0, 0.8)
                                                                f2[idm][index].SetParameter( 1, 10.0 )
                                                                f2[idm][index].SetParameter( 2, 0.5 )
                                                                f2[idm][index].SetParameter( 3, -30.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
                                                        elif(typ=="DATA"):
                                                                f2[idm][index].SetParameter( 0, 0.8)
                                                                f2[idm][index].SetParameter( 1, 10.0 )
                                                                f2[idm][index].SetParameter( 2, 0.2 )
                                                                f2[idm][index].SetParameter( 3, -25.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
                                                                f2[idm][index].SetParameter( 5, 1.0)
						elif(wp =="mediumTauMVA"):
                                                        if(typ == "MC"):
                                                                f2[idm][index].SetParameter( 0, 0.8)
                                                                f2[idm][index].SetParameter( 1, 10.0 )
                                                                f2[idm][index].SetParameter( 2, 0.3 )
                                                                f2[idm][index].SetParameter( 3, -30.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
                                                        elif(typ=="DATA"):
                                                                f2[idm][index].SetParameter( 0, 0.8)
                                                                f2[idm][index].SetParameter( 1, 5.0 )
                                                                f2[idm][index].SetParameter( 2, 0.2 )
                                                                f2[idm][index].SetParameter( 3, -30.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
						elif(wp =="looseTauMVA"):
                                                        if(typ == "MC"):
                                                                f2[idm][index].SetParameter( 0, 0.8)
                                                                f2[idm][index].SetParameter( 1, 10.0 )
                                                                f2[idm][index].SetParameter( 2, 0.3 )
                                                                f2[idm][index].SetParameter( 3, -25.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
                                                        elif(typ=="DATA"):
                                                                f2[idm][index].SetParameter( 0, 0.8)
                                                                f2[idm][index].SetParameter( 1, 10.0 )
                                                                f2[idm][index].SetParameter( 2, 0.3 )
                                                                f2[idm][index].SetParameter( 3, -30.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
						elif(wp =="vlooseTauMVA"):
							if(typ == "MC"):
                                                                f2[idm][index].SetParameter( 0, 0.8)
                                                                f2[idm][index].SetParameter( 1, 10.0 )
                                                                f2[idm][index].SetParameter( 2, 0.3 )
                                                                f2[idm][index].SetParameter( 3, -30.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
                                                        elif(typ=="DATA"):
                                                                f2[idm][index].SetParameter( 0, 0.8)
                                                                f2[idm][index].SetParameter( 1, 10.0 )
                                                                f2[idm][index].SetParameter( 2, 0.4 )
                                                                f2[idm][index].SetParameter( 3, -25.)
                                                                f2[idm][index].SetParameter( 4, 1.0 )
                                                                f2[idm][index].SetParameter( 5, 1.0)
						else:
							f2[idm][index].SetParameter( 0, 0.8)
							f2[idm][index].SetParameter( 1, 10.0)
							f2[idm][index].SetParameter( 2, 0.2)
							f2[idm][index].SetParameter( 3, -25.)
							f2[idm][index].SetParameter( 4, 1.0 )
							f2[idm][index].SetParameter( 5, 1.0)

			if(ipath == 2): # etau
				f1[index].SetParameter( 0, 1.0)
				f1[index].SetParameter( 1, 10.0 )
				f1[index].SetParameter( 2, 1.0 )
				f1[index].SetParameter( 3, -30.)
				f1[index].SetParameter( 4, 1.0)
				f1[index].SetParameter( 5, 1.0)
				for idm, DM in enumerate(tauDMs):
					if(idm ==0):
						f2[idm][index].SetParameter( 0, 0.8)
						f2[idm][index].SetParameter( 1, 10.0 )
						f2[idm][index].SetParameter( 2, 5.0 )
						f2[idm][index].SetParameter( 3, -25.)
						f2[idm][index].SetParameter( 4, 1.0 )
						f2[idm][index].SetParameter( 5, 1.0)
					elif(idm == 1):
						f2[idm][index].SetParameter( 0, 0.8)
						f2[idm][index].SetParameter( 1, 10.0 )
						f2[idm][index].SetParameter( 2, 5.0 )
						f2[idm][index].SetParameter( 3, -20.)
						f2[idm][index].SetParameter( 4, 1.0 )
						f2[idm][index].SetParameter( 5, 1.0)
					elif(DM =="dm10"):
						f2[idm][index].SetParameter( 0, 0.8)
						f2[idm][index].SetParameter( 1, 10.0 )
						f2[idm][index].SetParameter( 2, 0.4 )
						f2[idm][index].SetParameter( 3, -25.)
						f2[idm][index].SetParameter( 4, 1.0 )
						f2[idm][index].SetParameter( 5, 1.0)
			
			g_efficiency =TGraphAsymmErrors()
			h_errBand68.append(TH1F(histoname+"_CL68","histo of 0.68 confidence band", 480, 20, 500))
			g_errBand68.append(TGraphErrors())

			g_efficiency.BayesDivide(hPtNum[ipath][index][WPind],hPtDen[ipath][index][WPind]) #,"cl=0.683 b(1,1) mode")
			h_efficiency = getTH1FfromTGraphAsymmErrors(g_efficiency,"histo_" + trigger + "ErrorBand_" + wp +"_"+ typ )

			# write the histograms/graphs into the output ROOT file before the fit
			h_efficiency.Write("histo_"+ trigger +"Efficiency_" + wp +"_"+ typ)
			g_efficiency.Write("graph_"+ trigger +"Efficiency_" + wp +"_"+ typ)
			
			print "Fit is performed for", trigger, "trigger in", wp ,"WP for", typ 
			print "Fit parameters:", f1[index].GetParameter(0), f1[index].GetParameter(1), f1[index].GetParameter(2), f1[index].GetParameter(3), f1[index].GetParameter(4), f1[index].GetParameter(5)    
			
			fit_result = g_efficiency.Fit('f1'+ typ, 'S')
			h_errBand68[index], g_errBand68[index] = getConfidenceInterval(g_efficiency, h_errBand68[index], g_errBand68[index], 0.68)

			# Set the title of the histograms/graphs and their axes
			g_efficiency.SetTitle(trigger +"Path_" + wp +"_"+ typ)
			g_efficiency.GetYaxis().SetTitle("Efficiency")
			g_efficiency.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
			h_efficiency.SetTitle(trigger +"Path_" + wp +"_"+ typ)
			h_efficiency.GetYaxis().SetTitle("Efficiency")
			h_efficiency.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")

			h_errBand68[index].SetTitle(trigger +"Path_" + wp +"_"+ typ)
			h_errBand68[index].GetYaxis().SetTitle("Efficiency")
			h_errBand68[index].GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
			g_errBand68[index].SetTitle(trigger +"Path_" + wp +"_"+ typ)
			g_errBand68[index].GetYaxis().SetTitle("Efficiency")
			g_errBand68[index].GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
			
			# write the histograms/graphs into the output ROOT file after the fit
			g_efficiency.Write("grFit_"+ trigger +"Efficiency_" + wp +"_"+ typ)
			h_errBand68[index].Write("herrorBand_"+ trigger +"Path_" + wp +"_"+ typ)
			g_errBand68[index].Write("gerrorBand_"+ trigger +"Path_" + wp +"_"+ typ)
			fit_result.Write("fitResult_"+ trigger +"Path_" + wp +"_"+ typ)

			#======== Relative error of the fit: "fit +/- error/ fit " ===================
			relativeErrorUP = TGraphAsymmErrors()
			relativeErrorDown = TGraphAsymmErrors()
			relativeErrorUP, relativeErrorDown = createRelativeErrors(fit_result, g_efficiency, f1)

			relativeErrorUP.Write("relativeErrorUp_"+ trigger  +"Path_"  + wp +"_"+ typ)
			relativeErrorDown.Write("relativeErrorDown_"+ trigger +"Path_"   + wp +"_"+ typ)

			g_efficiencyDM0 = TGraphAsymmErrors()
			g_efficiencyDM1 = TGraphAsymmErrors()
			g_efficiencyDM0.BayesDivide(hPtNumDM[ipath][index][WPind][0],hPtDenDM[ipath][index][WPind][0])
			g_efficiencyDM1.BayesDivide(hPtNumDM[ipath][index][WPind][1],hPtDenDM[ipath][index][WPind][1])

			g_efficiencyDM0.Write("testgraph_"+ trigger +"Efficiency_dm0" +"_" + wp +"_"+ typ)
			g_efficiencyDM1.Write("testgraph_"+ trigger +"Efficiency_dm1"+"_" + wp +"_"+ typ)

			# per DM efficiencies
			for idm, DM in enumerate(tauDMs):
				g_efficiencyDM = TGraphAsymmErrors()
				g_efficiencyDM.BayesDivide(hPtNumDM[ipath][index][WPind][idm],hPtDenDM[ipath][index][WPind][idm])
				h_efficiencyDM = getTH1FfromTGraphAsymmErrors(g_efficiencyDM,"histo_" + trigger + "ErrorBand_" + DM +"_" + wp +"_"+ typ )

				fit_result2 = g_efficiencyDM.Fit('f2_'+ DM +"_" + typ, 'S')
				h_errBandDM68[idm][index], g_errBandDM68[idm][index] = getConfidenceInterval(g_efficiencyDM, h_errBandDM68[idm][index], g_errBandDM68[idm][index], 0.68)

				g_efficiencyDM.Write("graph_"+ trigger +"Efficiency_"+ DM +"_" + wp +"_"+ typ)
				h_efficiencyDM.Write("histo_"+ trigger +"Efficiency_"+ DM +"_" + wp +"_"+ typ)
				g_efficiencyDM.Write("grFit_"+ trigger +"Efficiency_"+ DM +"_" + wp +"_"+ typ)
				fit_result2.Write("fitResult_"+ trigger +"Path_" + DM +"_" + wp +"_"+ typ)
				h_errBandDM68[idm][index].Write("herrorBand_"+ trigger +"Path_" + DM +"_" + wp +"_"+ typ)
				g_errBandDM68[idm][index].Write("gerrorBand_"+ trigger +"Path_" + DM +"_" + wp +"_"+ typ)

				#======== Relative error of the fit: "fit +/- error/ fit " ===================
				relativeErrorDMUP = TGraphAsymmErrors()
				relativeErrorDMDown = TGraphAsymmErrors()
				relativeErrorDMUP, relativeErrorDMDown = createRelativeErrors(fit_result2, g_efficiencyDM, f2[idm])

				relativeErrorDMUP.Write("relativeErrorUp_"+ trigger  +"Path_"  + DM + wp +"_"+ typ)
				relativeErrorDMDown.Write("relativeErrorDown_"+ trigger  +"Path_" + "_" + DM + wp +"_"+ typ)

		# Getting Scale Factors
		SF = TGraphAsymmErrors()
		SF = getScaleFactor(f1)
		SF.Write( 'ScaleFactor_'+ trigger + '_' + wp + '_DataMC')
		
		# Getting Scale Factors with Errors
		SF_err = TGraphAsymmErrors()
		SF_errFromDivide = TGraphAsymmErrors()
		Effi2 = TH1F() # ( "" , "", len(edges)-1, array('f',edges))
		SF_err, Effi2 = getScaleFactorWithPropagatedError(f1, h_errBand68 )
		SF_err.Write( 'ScaleFactorWithPropogatedError_'+ trigger + '_' + wp + '_DataMC')
		Effi2.Write( 'ScaleFactorWithErrorFromDivide_'+ trigger + '_' + wp + '_DataMC')

		SF_err2 = TGraphAsymmErrors()
		SF_err2 = getScaleFactorWithShiftedError(f1, h_errBand68 )
		SF_err2.Write( 'ScaleFactorWithShiftedError_'+ trigger + '_' + wp + '_DataMC')

		# Getting Scale Factors per decay mode
		for idm, DM in enumerate(tauDMs):
			SF_dm = TGraphErrors()
			SF_dm = getScaleFactor(f2[idm])
			SF_dm.Write( 'ScaleFactor_'+ trigger + "_" + DM + '_' + wp + '_DataMC')

			# Getting Scale Factors with Errors
			SFdm_err = TGraphAsymmErrors()
			SFdm_errFromDivide = TGraphAsymmErrors()
			Effidm = TH1F()
			SFdm_err, Effidm = getScaleFactorWithPropagatedError(f2[idm], h_errBandDM68[idm] )
			SFdm_err.Write( 'ScaleFactorWithPropogatedError_'+ trigger + "_" + DM + '_' + wp + '_DataMC')
			Effidm.Write( 'ScaleFactorWithErrorFromDivide_'+ trigger + "_" + DM +'_' + wp + '_DataMC')

file.Close()
print "The output ROOT file has been created: ../data/" + outputname


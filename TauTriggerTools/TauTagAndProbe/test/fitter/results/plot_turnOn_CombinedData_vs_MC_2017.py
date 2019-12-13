import ROOT
import TurnOnPlot_DATA as TurnOnPlot

### Edit here ###

# TRIGGERS MUST BE DECLARED
triggers = ["HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v","HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v"]

#triggers = ["ETauTriggerPath_IsoMu20_LooseChargedIsoPFTau27_plusL1Tau26andHLTTau30","MuTauTriggerPath_IsoMu20_LooseChargedIsoPFTau27_plusL1Tau32", "DiTauTriggerPath_TightTau35orMediumTau40TightIDorTightTau40_plusL1Tau32"]
# PLOT TITLES
#plotTitles = ["HLT MediumIsoPFTau32 Data - MC", "HLT MediumIsoPFTau20 Data - MC"]
# ROOT FILE CONTAINING THE DATA
dataFileName = "FitResults/2018_01_14/fitOutput_Data_MuTau2017BCDEF_SFpaths_SSsubtraction_vtightTauMVAWP.root"
# ROOT FILE CONTAINING THE MC
mcFileName = "FitResults/2018_01_14/fitOutput_MC_MuTau2017_DYJetsFall17_nomPlusExt_SFpaths_OStaugenmatchPositive_vtightTauMVAWP.root"

### Do not edit from here ###

#open turn on file
inputFile_Data = ROOT.TFile.Open(dataFileName)
inputFile_MC = ROOT.TFile.Open(mcFileName)

histo_Data = []
histo_MC = []
fit_Data = []
fit_MC = []
turnon_Data = []
turnon_MC = []
plots = []

for trigger in triggers:
    histo_Data.append(inputFile_Data.Get("histo_" + trigger))
    histo_Data[-1].__class__ = ROOT.RooHist
    histo_MC.append(inputFile_MC.Get("histo_" + trigger))
    histo_MC[-1].__class__ = ROOT.RooHist
    fit_Data.append(inputFile_Data.Get("fit_" + trigger))
    fit_Data[-1].__class__ = ROOT.RooCurve
    fit_MC.append(inputFile_MC.Get("fit_" + trigger))
    fit_MC[-1].__class__ = ROOT.RooCurve
    turnon_Data.append(TurnOnPlot.TurnOn(Name="Stage2_Data", Histo=histo_Data[-1], Fit=fit_Data[-1],
                                    MarkerColor=ROOT.kBlue, MarkerStyle=20, LineColor=ROOT.kBlue,LineStyle=1,
                                    Legend="Data")) # 2017, Run[B-F]"))
    turnon_MC.append(TurnOnPlot.TurnOn(Name="Stage2_MC", Histo=histo_MC[-1], Fit=fit_MC[-1],
                                   MarkerColor=ROOT.kRed, MarkerStyle=20, LineColor=ROOT.kRed,LineStyle=1,
                                   Legend="Simulation"))#MC 2017, DYJets"))
    plots.append(TurnOnPlot.TurnOnPlot(TriggerName = trigger + "Data - MC"))
    plots[-1].name = "turnOn_2017_Data_vs_MC_vtightTauMVAWP_" + trigger
    plots[-1].xRange = (20,500)
    #plots[-1].legendPosition = (0.6,0.2,0.9,0.4)
    plots[-1].legendPosition = (0.55,0.3,0.85,0.45)
    plots[-1].addTurnOn(turnon_MC[-1])
    plots[-1].addTurnOn(turnon_Data[-1])
    

	
canvas = []
for plot in plots:
    canvas.append(plot.plot())
	#print "plot:", plot.name
	#canvas.Update()
	#canvas.Print(plots[-1].name , "png")

inputFile_Data.Close()
inputFile_MC.Close()

raw_input()

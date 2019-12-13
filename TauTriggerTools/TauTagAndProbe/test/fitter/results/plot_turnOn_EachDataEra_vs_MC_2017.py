import ROOT
import TurnOnPlot_DATA as TurnOnPlot

### Edit here ###

# TRIGGERS MUST BE DECLARED
triggers = ["HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v", "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v",  "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"]
# PLOT TITLES
#plotTitles = ["HLT MediumIsoPFTau32 Data - MC", "HLT MediumIsoPFTau20 Data - MC"]
# ROOT FILE CONTAINING THE DATA
data17BFileName = "TurnOn_Sync/2017_11_23/fitOutput_Data_MuTau2017B_PRv1v2_23112017_v2.root"
data17CFileName = "TurnOn_Sync/2017_11_23/fitOutput_Data_MuTau2017C_PRv1v2v3_2311201_v2.root"
data17DFileName = "TurnOn_Sync/2017_11_23/fitOutput_Data_MuTau2017D_PRv1_23112017_v2.root"
data17EFileName = "TurnOn_Sync/2017_11_23/fitOutput_Data_MuTau2017E_PRv1_23112017_v2.root"
data17FFileName = "TurnOn_Sync/2017_11_23/fitOutput_Data_MuTau2017F_PRv1_23112017_v2.root"
# ROOT FILE CONTAINING THE MC
mcFileName = "TurnOn_Sync/2017_11_23/fitOutput_MC_MuTau2017v10_ext1_v1v2_pileup_23112017_DitauFit.root"


### Do not edit from here ###

#open turn on file
inputFile_Data2017B = ROOT.TFile.Open(data17BFileName)
inputFile_Data2017C = ROOT.TFile.Open(data17CFileName)
inputFile_Data2017D = ROOT.TFile.Open(data17DFileName)
inputFile_Data2017E = ROOT.TFile.Open(data17EFileName)
inputFile_Data2017F = ROOT.TFile.Open(data17FFileName)
inputFile_MC = ROOT.TFile.Open(mcFileName)

histo_Data2017B = []
histo_Data2017C = []
histo_Data2017D = []
histo_Data2017E = []
histo_Data2017F = []
histo_MC = []
fit_Data2017B = []
fit_Data2017C = []
fit_Data2017D = []
fit_Data2017E = []
fit_Data2017F = []
fit_MC = []
turnon_Data2017B = []
turnon_Data2017C = []
turnon_Data2017D = []
turnon_Data2017E = []
turnon_Data2017F = []
turnon_MC = []
plots = []

for trigger in triggers:
    histo_Data2017B.append(inputFile_Data2017B.Get("histo_" + trigger))
    histo_Data2017B[-1].__class__ = ROOT.RooHist
    histo_Data2017C.append(inputFile_Data2017C.Get("histo_" + trigger))
    histo_Data2017C[-1].__class__ = ROOT.RooHist
    histo_Data2017D.append(inputFile_Data2017D.Get("histo_" + trigger))
    histo_Data2017D[-1].__class__ = ROOT.RooHist
    histo_Data2017E.append(inputFile_Data2017E.Get("histo_" + trigger))
    histo_Data2017E[-1].__class__ = ROOT.RooHist
    histo_Data2017F.append(inputFile_Data2017F.Get("histo_" + trigger))
    histo_Data2017F[-1].__class__ = ROOT.RooHist
    histo_MC.append(inputFile_MC.Get("histo_" + trigger))
    histo_MC[-1].__class__ = ROOT.RooHist
    fit_Data2017B.append(inputFile_Data2017B.Get("fit_" + trigger))
    fit_Data2017B[-1].__class__ = ROOT.RooCurve
    fit_Data2017C.append(inputFile_Data2017C.Get("fit_" + trigger))
    fit_Data2017C[-1].__class__ = ROOT.RooCurve
    fit_Data2017D.append(inputFile_Data2017D.Get("fit_" + trigger))
    fit_Data2017D[-1].__class__ = ROOT.RooCurve
    fit_Data2017E.append(inputFile_Data2017E.Get("fit_" + trigger))
    fit_Data2017E[-1].__class__ = ROOT.RooCurve
    fit_Data2017F.append(inputFile_Data2017F.Get("fit_" + trigger))
    fit_Data2017F[-1].__class__ = ROOT.RooCurve
    fit_MC.append(inputFile_MC.Get("fit_" + trigger))
    fit_MC[-1].__class__ = ROOT.RooCurve
    turnon_Data2017B.append(TurnOnPlot.TurnOn(Name="Stage2_DataB", Histo=histo_Data2017B[-1], Fit=fit_Data2017B[-1],
                                          MarkerColor=ROOT.kBlue, MarkerStyle=20, LineColor=ROOT.kBlue,LineStyle=1,
                                          Legend="2017 Data B"))
    turnon_Data2017C.append(TurnOnPlot.TurnOn(Name="Stage2_DataC", Histo=histo_Data2017C[-1], Fit=fit_Data2017C[-1],
                                          MarkerColor=ROOT.kGreen+2, MarkerStyle=20, LineColor=ROOT.kGreen+2,LineStyle=1,
                                          Legend="2017 Data C"))
    turnon_Data2017D.append(TurnOnPlot.TurnOn(Name="Stage2_DataD", Histo=histo_Data2017D[-1], Fit=fit_Data2017D[-1],
                                          MarkerColor=ROOT.kAzure+8, MarkerStyle=20, LineColor=ROOT.kAzure+8,LineStyle=1,
                                          Legend="2017 Data D"))
    turnon_Data2017E.append(TurnOnPlot.TurnOn(Name="Stage2_DataE", Histo=histo_Data2017E[-1], Fit=fit_Data2017E[-1],
                                          MarkerColor=ROOT.kOrange+3, MarkerStyle=20, LineColor=ROOT.kOrange+3,LineStyle=1,
                                          Legend="2017 Data E"))
    turnon_Data2017F.append(TurnOnPlot.TurnOn(Name="Stage2_DataF", Histo=histo_Data2017F[-1], Fit=fit_Data2017F[-1],
                                          MarkerColor=ROOT.kViolet-2, MarkerStyle=20, LineColor=ROOT.kViolet-2,LineStyle=1,
                                          Legend="2017 Data F"))
    turnon_MC.append(TurnOnPlot.TurnOn(Name="Stage2_MC", Histo=histo_MC[-1], Fit=fit_MC[-1],
                                   MarkerColor=ROOT.kRed, MarkerStyle=20, LineColor=ROOT.kRed,LineStyle=1,
                                   Legend="Simulation"))
    plots.append(TurnOnPlot.TurnOnPlot(TriggerName = trigger + "Data - MC"))
    plots[-1].name = "turnOn_2017_EachDataEra_vs_MC_" + trigger
    plots[-1].xRange = (20,500)
    #plots[-1].legendPosition = (0.6,0.2,0.9,0.4)
    plots[-1].legendPosition = (0.6,0.2,0.9,0.45)
    plots[-1].addTurnOn(turnon_MC[-1])
    plots[-1].addTurnOn(turnon_Data2017B[-1])
    plots[-1].addTurnOn(turnon_Data2017C[-1])
    plots[-1].addTurnOn(turnon_Data2017D[-1])
    plots[-1].addTurnOn(turnon_Data2017E[-1])
    plots[-1].addTurnOn(turnon_Data2017F[-1])

	
canvas = []
for plot in plots:
    canvas.append(plot.plot())
	#print "plot:", plot.name
	#canvas.Update()
	#canvas.Print(plots[-1].name , "png")

inputFile_Data2017B.Close()
inputFile_Data2017C.Close()
inputFile_Data2017D.Close()
inputFile_Data2017E.Close()
inputFile_Data2017F.Close()
inputFile_MC.Close()

raw_input()

import ROOT
import TurnOnPlot_DATA as TurnOnPlot

### Edit here ###

# TRIGGERS MUST BE DECLARED
#triggers = ["HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_v", "HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v", "Pt_26GeV", "Pt_34GeV"]
#triggers = ["Pt_30GeV"]
# PLOT TITLES
#plotTitles = ["HLT MediumIsoPFTau32 Barrel - Endcaps", "HLT MediumIsoPFTau20 Barrel - Endcaps"]
# ROOT FILE CONTAINING THE Barrel
fileName = "FittedTurnOn_Final_MC.root"
fileName20GeV = "FittedTurnOn_L1_20GeV.root"

### Do not edit from here ###

#open turn on file
inputFile= ROOT.TFile.Open(fileName)
inputFile20GeV= ROOT.TFile.Open(fileName20GeV)

histo_SingleL1 = []
histo_WithL1 = []
fit_SingleL1 = []
fit_WithL1 = []
histo_20GeV = []
fit_20GeV = []

turnon_SingleL1 = []
turnon_WithL1 = []
turnon_20GeV = []
plots = []

#for trigger in triggers:
histo_SingleL1.append(inputFile.Get("histo_Stage2_All_vs_HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v"))
histo_SingleL1[-1].__class__ = ROOT.RooHist
histo_WithL1.append(inputFile.Get("histo_Stage2_All_vs_HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v"))
histo_WithL1[-1].__class__ = ROOT.RooHist
histo_20GeV.append(inputFile20GeV.Get("histo_Stage2_All_vs_Pt_20GeV"))
histo_20GeV[-1].__class__ = ROOT.RooHist
fit_SingleL1.append(inputFile.Get("fit_Stage2_All_vs_HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v"))
fit_SingleL1[-1].__class__ = ROOT.RooCurve
fit_WithL1.append(inputFile.Get("fit_Stage2_All_vs_HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v"))
fit_WithL1[-1].__class__ = ROOT.RooCurve
fit_20GeV.append(inputFile20GeV.Get("fit_Stage2_All_vs_Pt_20GeV"))
fit_20GeV[-1].__class__ = ROOT.RooCurve
turnon_SingleL1.append(TurnOnPlot.TurnOn(Name="Stage2_SingleL1", Histo=histo_SingleL1[-1], Fit=fit_SingleL1[-1],
                                MarkerColor=ROOT.kBlue, MarkerStyle=20, LineColor=ROOT.kBlue,LineStyle=1,
                                Legend="No L1T seed"))
turnon_WithL1.append(TurnOnPlot.TurnOn(Name="Stage2_WithL1", Histo=histo_WithL1[-1], Fit=fit_WithL1[-1],
                               MarkerColor=ROOT.kRed, MarkerStyle=20, LineColor=ROOT.kRed,LineStyle=1,
                               Legend="With L1T seed"))
turnon_20GeV.append(TurnOnPlot.TurnOn(Name="Stage2_20GeV", Histo=histo_20GeV[-1], Fit=fit_20GeV[-1],
                               MarkerColor=ROOT.kGreen, MarkerStyle=20, LineColor=ROOT.kGreen,LineStyle=1,
                               Legend="L1T 20 GeV turn-on"))
plots.append(TurnOnPlot.TurnOnPlot(TriggerName = "HLT LooseIsoPFTau20" + "No L1T seed - With L1T seed"))
plots[-1].name = "turnOn_SingleL1_WithL1_HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v"
plots[-1].xRange = (10,109.9)
#plots[-1].legendPosition = (0.6,0.2,0.9,0.4)
plots[-1].legendPosition = (0.6,0.2,0.9,0.4)
plots[-1].addTurnOn(turnon_SingleL1[-1])
plots[-1].addTurnOn(turnon_WithL1[-1])
plots[-1].addTurnOn(turnon_20GeV[-1])

canvas = []
for plot in plots:
    canvas.append(plot.plot())

inputFile.Close()

raw_input()

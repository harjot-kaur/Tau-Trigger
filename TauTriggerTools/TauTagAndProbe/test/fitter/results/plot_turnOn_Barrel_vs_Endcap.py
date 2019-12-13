import ROOT
import TurnOnPlot_DATA as TurnOnPlot

### Edit here ###

# TRIGGERS MUST BE DECLARED
triggers = ["HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_v", "HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v", "Pt_26GeV", "Pt_34GeV"]
#triggers = ["Pt_30GeV"]
# PLOT TITLES
#plotTitles = ["HLT MediumIsoPFTau32 Barrel - Endcaps", "HLT MediumIsoPFTau20 Barrel - Endcaps"]
# ROOT FILE CONTAINING THE Barrel
fileName = "FittedTurnOn_Final_MC.root"

### Do not edit from here ###

#open turn on file
inputFile= ROOT.TFile.Open(fileName)

histo_Barrel = []
histo_Endcaps = []
fit_Barrel = []
fit_Endcaps = []
turnon_Barrel = []
turnon_Endcaps = []
plots = []

for trigger in triggers:
    histo_Barrel.append(inputFile.Get("histo_Stage2_Barrel_vs_" + trigger))
    histo_Barrel[-1].__class__ = ROOT.RooHist
    histo_Endcaps.append(inputFile.Get("histo_Stage2_Endcaps_vs_" + trigger))
    histo_Endcaps[-1].__class__ = ROOT.RooHist
    fit_Barrel.append(inputFile.Get("fit_Stage2_Barrel_vs_" + trigger))
    fit_Barrel[-1].__class__ = ROOT.RooCurve
    fit_Endcaps.append(inputFile.Get("fit_Stage2_Endcaps_vs_" + trigger))
    fit_Endcaps[-1].__class__ = ROOT.RooCurve
    turnon_Barrel.append(TurnOnPlot.TurnOn(Name="Stage2_Barrel", Histo=histo_Barrel[-1], Fit=fit_Barrel[-1],
                                    MarkerColor=ROOT.kBlue, MarkerStyle=20, LineColor=ROOT.kBlue,LineStyle=1,
                                    Legend="Barrel"))
    turnon_Endcaps.append(TurnOnPlot.TurnOn(Name="Stage2_Endcaps", Histo=histo_Endcaps[-1], Fit=fit_Endcaps[-1],
                                   MarkerColor=ROOT.kRed, MarkerStyle=20, LineColor=ROOT.kRed,LineStyle=1,
                                   Legend="Endcaps"))
    plots.append(TurnOnPlot.TurnOnPlot(TriggerName = trigger + "Barrel - Endcaps"))
    plots[-1].name = "turnOn_Barrel_Endcaps_" + trigger
    plots[-1].xRange = (10,109.9)
    #plots[-1].legendPosition = (0.6,0.2,0.9,0.4)
    plots[-1].legendPosition = (0.6,0.2,0.9,0.4)
    plots[-1].addTurnOn(turnon_Barrel[-1])
    plots[-1].addTurnOn(turnon_Endcaps[-1])

canvas = []
for plot in plots:
    canvas.append(plot.plot())

inputFile.Close()

raw_input()

import ROOT
import TurnOnPlot_DATA as TurnOnPlot

### Edit here ###

# TRIGGERS MUST BE DECLARED
triggers = ["HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_v", "HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v", "Pt_26GeV", "Pt_34GeV"]
# PLOT TITLES
#plotTitles = ["HLT MediumIsoPFTau32 Data - MC", "HLT MediumIsoPFTau20 Data - MC"]
# ROOT FILE CONTAINING THE DATA
dataFileName = "FittedTurnOn_Final_Data_0_500GeV.root"
# ROOT FILE CONTAINING THE MC
mcFileName = "FittedTurnOn_Final_MC_0_500GeV.root"

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
    histo_Data.append(inputFile_Data.Get("histo_Stage2_All_vs_" + trigger))
    histo_Data[-1].__class__ = ROOT.RooHist
    histo_MC.append(inputFile_MC.Get("histo_Stage2_All_vs_" + trigger))
    histo_MC[-1].__class__ = ROOT.RooHist
    fit_Data.append(inputFile_Data.Get("fit_Stage2_All_vs_" + trigger))
    fit_Data[-1].__class__ = ROOT.RooCurve
    fit_MC.append(inputFile_MC.Get("fit_Stage2_All_vs_" + trigger))
    fit_MC[-1].__class__ = ROOT.RooCurve
    turnon_Data.append(TurnOnPlot.TurnOn(Name="Stage2_Data", Histo=histo_Data[-1], Fit=fit_Data[-1],
                                    MarkerColor=ROOT.kBlue, MarkerStyle=20, LineColor=ROOT.kBlue,LineStyle=1,
                                    Legend="Data"))
    turnon_MC.append(TurnOnPlot.TurnOn(Name="Stage2_MC", Histo=histo_MC[-1], Fit=fit_MC[-1],
                                   MarkerColor=ROOT.kRed, MarkerStyle=20, LineColor=ROOT.kRed,LineStyle=1,
                                   Legend="Simulation"))
    plots.append(TurnOnPlot.TurnOnPlot(TriggerName = trigger + "Data - MC"))
    plots[-1].name = "turnOn_Data_MC_" + trigger
    plots[-1].xRange = (10,120)
    #plots[-1].legendPosition = (0.6,0.2,0.9,0.4)
    plots[-1].legendPosition = (0.6,0.2,0.9,0.4)
    plots[-1].addTurnOn(turnon_Data[-1])
    plots[-1].addTurnOn(turnon_MC[-1])

canvas = []
for plot in plots:
    canvas.append(plot.plot())


inputFile_Data.Close()
inputFile_MC.Close()

raw_input()

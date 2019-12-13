import ROOT
import TurnOnPlot_DATA as TurnOnPlot

### Edit here ###

# TRIGGERS MUST BE DECLARED
triggers = ["HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v", "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v", "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"]
# PLOT TITLES
#plotTitles = ["HLT MediumIsoPFTau32 Data - MC", "HLT MediumIsoPFTau20 Data - MC"]
# ROOT FILE CONTAINING THE DATA
dataBFileName = "TestTurnOn/fitOutput_Data_MuTau2017B.root"
dataCFileName = "TestTurnOn/fitOutput_Data_MuTau2017C.root"
# ROOT FILE CONTAINING THE MC
mcFileName = "TestTurnOn/fitOutput_DY_BadPixGT.root"

### Do not edit from here ###

#open turn on file
inputFile_DataB = ROOT.TFile.Open(dataBFileName)
inputFile_DataC = ROOT.TFile.Open(dataCFileName)
inputFile_MC = ROOT.TFile.Open(mcFileName)

histo_DataB = []
histo_DataC = []
histo_MC = []
fit_DataB = []
fit_DataC = []
fit_MC = []
turnon_DataB = []
turnon_DataC = []
turnon_MC = []
plots = []

for trigger in triggers:
    histo_DataB.append(inputFile_DataB.Get("histo_" + trigger))
    histo_DataB[-1].__class__ = ROOT.RooHist
    histo_DataC.append(inputFile_DataC.Get("histo_" + trigger))
    histo_DataC[-1].__class__ = ROOT.RooHist
    histo_MC.append(inputFile_MC.Get("histo_" + trigger))
    histo_MC[-1].__class__ = ROOT.RooHist
    fit_DataB.append(inputFile_DataB.Get("fit_" + trigger))
    fit_DataB[-1].__class__ = ROOT.RooCurve
    fit_DataC.append(inputFile_DataC.Get("fit_" + trigger))
    fit_DataC[-1].__class__ = ROOT.RooCurve
    fit_MC.append(inputFile_MC.Get("fit_" + trigger))
    fit_MC[-1].__class__ = ROOT.RooCurve
    turnon_DataB.append(TurnOnPlot.TurnOn(Name="Stage2_DataB", Histo=histo_DataB[-1], Fit=fit_DataB[-1],
                                          MarkerColor=ROOT.kBlue, MarkerStyle=20, LineColor=ROOT.kBlue,LineStyle=1,
                                          Legend="Data up to 18/07/2017 (4.6 fb^{-1})"))
    turnon_DataC.append(TurnOnPlot.TurnOn(Name="Stage2_DataC", Histo=histo_DataC[-1], Fit=fit_DataC[-1],
                                          MarkerColor=ROOT.kRed, MarkerStyle=20, LineColor=ROOT.kRed,LineStyle=1,
                                          Legend="Data from 18/07/2017 (1.2 fb^{-1})"))
    turnon_MC.append(TurnOnPlot.TurnOn(Name="Stage2_MC", Histo=histo_MC[-1], Fit=fit_MC[-1],
                                   MarkerColor=ROOT.kRed, MarkerStyle=20, LineColor=ROOT.kRed,LineStyle=1,
                                   Legend="Simulation"))
    plots.append(TurnOnPlot.TurnOnPlot(TriggerName = trigger + "Data B - Data C - MC"))
    plots[-1].name = "turnOn_Data_2017B_2017C_DY_" + trigger
    plots[-1].xRange = (20,500)
    #if(trigger=="HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v"):
    #    plots[-1].xRange = (20,500)
    #plots[-1].legendPosition = (0.6,0.2,0.9,0.4)
    plots[-1].legendPosition = (0.4,0.2,0.9,0.4)
    #plots[-1].addTurnOn(turnon_MC[-1])
    plots[-1].addTurnOn(turnon_DataB[-1])
    plots[-1].addTurnOn(turnon_DataC[-1])

canvas = []
for plot in plots:
    canvas.append(plot.plot())


inputFile_DataB.Close()
inputFile_DataC.Close()
inputFile_MC.Close()

raw_input()

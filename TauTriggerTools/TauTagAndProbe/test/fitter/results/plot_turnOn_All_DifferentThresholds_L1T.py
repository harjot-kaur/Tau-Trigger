import ROOT
import TurnOnPlot_DATA as TurnOnPlot


plots = []
plots.append(TurnOnPlot.TurnOnPlot(TriggerName="L1 turn-on curves, iso + non-iso"))
plots[-1].name = "turnOn_All"
plots[-1].xRange = (10,109.9)
#plots[-1].legendPosition = (0.6,0.2,0.9,0.4)
plots[-1].legendPosition = (0.6,0.2,0.9,0.4)

#EE
plots.append(TurnOnPlot.TurnOnPlot(TriggerName="L1 turn-on curves, iso"))
plots[-1].name = "turnOn_All_iso"
plots[-1].xRange = (10,109.9)
#plots[-1].legendPosition = (0.6,0.2,0.9,0.4)
plots[-1].legendPosition = (0.6,0.2,0.9,0.4)

#open turn on file
inputFile = ROOT.TFile.Open("FittedTurnOn_Final_MC.root")

histo_26GeV = inputFile.Get("histo_Stage2_All_vs_Pt_26GeV")
histo_26GeV.__class__ = ROOT.RooHist
histo_30GeV = inputFile.Get("histo_Stage2_All_vs_Pt_30GeV")
histo_30GeV.__class__ = ROOT.RooHist
histo_34GeV = inputFile.Get("histo_Stage2_All_vs_Pt_34GeV")
histo_34GeV.__class__ = ROOT.RooHist

fit_26GeV   = inputFile.Get("fit_Stage2_All_vs_Pt_26GeV")
fit_26GeV.__class__ = ROOT.RooCurve
fit_30GeV   = inputFile.Get("fit_Stage2_All_vs_Pt_30GeV")
fit_30GeV.__class__ = ROOT.RooCurve
fit_34GeV   = inputFile.Get("fit_Stage2_All_vs_Pt_34GeV")
fit_34GeV.__class__ = ROOT.RooCurve

turnon_26GeV = TurnOnPlot.TurnOn(Name="turnOn_All_26", Histo=histo_26GeV, Fit=fit_26GeV,
                                    MarkerColor=ROOT.kBlack, MarkerStyle=20, LineColor=ROOT.kBlack,LineStyle=1,
                                    Legend="pt > 26 GeV")

turnon_30GeV = TurnOnPlot.TurnOn(Name="turnOn_All_30", Histo=histo_30GeV, Fit=fit_30GeV,
                                    MarkerColor=ROOT.kRed, MarkerStyle=20, LineColor=ROOT.kRed,LineStyle=1,
                                    Legend="pt > 30 GeV")

turnon_34GeV = TurnOnPlot.TurnOn(Name="turnOn_All_34", Histo=histo_34GeV, Fit=fit_34GeV,
                                    MarkerColor=ROOT.kBlue, MarkerStyle=20, LineColor=ROOT.kBlue,LineStyle=1,
                                    Legend="pt > 34 GeV")


plots[0].addTurnOn(turnon_26GeV)
plots[0].addTurnOn(turnon_30GeV)
plots[0].addTurnOn(turnon_34GeV)


histo_26GeV_iso = inputFile.Get("histo_Stage2_All_vs_Pt_26GeV_iso")
histo_26GeV_iso.__class__ = ROOT.RooHist
histo_30GeV_iso = inputFile.Get("histo_Stage2_All_vs_Pt_30GeV_iso")
histo_30GeV_iso.__class__ = ROOT.RooHist
histo_34GeV_iso = inputFile.Get("histo_Stage2_All_vs_Pt_34GeV_iso")
histo_34GeV_iso.__class__ = ROOT.RooHist

fit_26GeV_iso   = inputFile.Get("fit_Stage2_All_vs_Pt_26GeV_iso")
fit_26GeV_iso.__class__ = ROOT.RooCurve
fit_30GeV_iso   = inputFile.Get("fit_Stage2_All_vs_Pt_30GeV_iso")
fit_30GeV_iso.__class__ = ROOT.RooCurve
fit_34GeV_iso   = inputFile.Get("fit_Stage2_All_vs_Pt_34GeV_iso")
fit_34GeV_iso.__class__ = ROOT.RooCurve

turnon_26GeV_iso = TurnOnPlot.TurnOn(Name="turnOn_All_26_iso", Histo=histo_26GeV_iso, Fit=fit_26GeV_iso,
                                    MarkerColor=ROOT.kBlack, MarkerStyle=20, LineColor=ROOT.kBlack,LineStyle=1,
                                    Legend="pt > 26 GeV")

turnon_30GeV_iso = TurnOnPlot.TurnOn(Name="turnOn_All_30_iso", Histo=histo_30GeV_iso, Fit=fit_30GeV_iso,
                                    MarkerColor=ROOT.kRed, MarkerStyle=20, LineColor=ROOT.kRed,LineStyle=1,
                                    Legend="pt > 30 GeV")

turnon_34GeV_iso = TurnOnPlot.TurnOn(Name="turnOn_All_34_iso", Histo=histo_34GeV_iso, Fit=fit_34GeV_iso,
                                    MarkerColor=ROOT.kBlue, MarkerStyle=20, LineColor=ROOT.kBlue,LineStyle=1,
                                    Legend="pt > 34 GeV")


plots[1].addTurnOn(turnon_26GeV_iso)
plots[1].addTurnOn(turnon_30GeV_iso)
plots[1].addTurnOn(turnon_34GeV_iso)

canvas = []
for plot in plots:
    canvas.append(plot.plot())


inputFile.Close()

raw_input()

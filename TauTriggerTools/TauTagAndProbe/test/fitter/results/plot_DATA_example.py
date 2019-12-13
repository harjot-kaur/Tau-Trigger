import ROOT
import TurnOnPlot_DATA as TurnOnPlot


plots = []
plots.append(TurnOnPlot.TurnOnPlot())
plots[-1].name = "turnon_EB_plot"
plots[-1].xRange = (10,109.9)
#plots[-1].legendPosition = (0.6,0.2,0.9,0.4)
plots[-1].legendPosition = (0.6,0.2,0.9,0.4)

#EE
plots.append(TurnOnPlot.TurnOnPlot())
plots[-1].name = "turnon_EE_plot"
plots[-1].xRange = (10,109.9)
#plots[-1].legendPosition = (0.6,0.2,0.9,0.4)
plots[-1].legendPosition = (0.6,0.2,0.9,0.4)

#open Non-Calib plot
inputFile = ROOT.TFile.Open("./TurnOnDataStage2NoCalibAll/Tau_stage1_stage2_EB_EE_All_Iso.root")
	
histo_EB = inputFile.Get("histo_Stage2_Barrel_vs_Pt")
histo_EB.__class__ = ROOT.RooHist
histo_EE = inputFile.Get("histo_Stage2_Endcaps_vs_Pt")
histo_EE.__class__ = ROOT.RooHist

	
fit_EB   = inputFile.Get("fit_Stage2_Barrel_vs_Pt")
fit_EB.__class__ = ROOT.RooCurve
fit_EE   = inputFile.Get("fit_Stage2_Endcaps_vs_Pt")
fit_EE.__class__ = ROOT.RooCurve

turnon_EB = TurnOnPlot.TurnOn(Name="Stage2_EB_noIso", Histo=histo_EB, Fit=fit_EB,
                                    MarkerColor=ROOT.kBlack, MarkerStyle=20, LineColor=ROOT.kBlack,LineStyle=1,
                                    Legend="Barrel - uncalibrated")
	
turnon_EE = TurnOnPlot.TurnOn(Name="Stage2_EE_noIso", Histo=histo_EE, Fit=fit_EE,
                                   MarkerColor=ROOT.kBlack, MarkerStyle=20, LineColor=ROOT.kBlack,LineStyle=1,
                                   Legend="Endcap - uncalibrated")


plots[0].addTurnOn(turnon_EB)
plots[1].addTurnOn(turnon_EE)




#open Calib plot
inputFile = ROOT.TFile.Open("./TurnOnDataStage2CalibAll_13Gen/Tau_stage1_stage2_EB_EE_All_Iso.root")
	
histo_EB = inputFile.Get("histo_Stage2_Barrel_vs_Pt")
histo_EB.__class__ = ROOT.RooHist
histo_EE = inputFile.Get("histo_Stage2_Endcaps_vs_Pt")
histo_EE.__class__ = ROOT.RooHist

	
fit_EB   = inputFile.Get("fit_Stage2_Barrel_vs_Pt")
fit_EB.__class__ = ROOT.RooCurve
fit_EE   = inputFile.Get("fit_Stage2_Endcaps_vs_Pt")
fit_EE.__class__ = ROOT.RooCurve

turnon_EB = TurnOnPlot.TurnOn(Name="Stage2_EB_noIso", Histo=histo_EB, Fit=fit_EB,
                                    MarkerColor=ROOT.kRed, MarkerStyle=20, LineColor=ROOT.kRed,LineStyle=1,
                                    Legend="Barrel - calibrated")
	
turnon_EE = TurnOnPlot.TurnOn(Name="Stage2_EE_noIso", Histo=histo_EE, Fit=fit_EE,
                                   MarkerColor=ROOT.kRed, MarkerStyle=20, LineColor=ROOT.kRed,LineStyle=1,
                                   Legend="Endcap - calibrated")


plots[0].addTurnOn(turnon_EB)
plots[1].addTurnOn(turnon_EE)


canvas = []
for plot in plots:
    canvas.append(plot.plot())

inputFile.Close()


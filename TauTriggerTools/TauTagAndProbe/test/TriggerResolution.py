#!/usr/bin/python

from ROOT import *
from array import array

gStyle.SetOptStat(1111)

fIn = TFile.Open('NTuple_10Ago_Riccardo.root')
tree = fIn.Get('Ntuplizer/TagAndProbe')

l1tResolution = TH1I("L1TriggerReso1", "L1T resolution", 200, -1, 1)
hltResolution = TH1I("HLTriggerReso1", "HLT resolution", 200, -1, 1)

for iEv in range (0, tree.GetEntries()):
    tree.GetEntry(iEv)
    if tree.foundJet != 1: continue
    ptOffline = tree.tauPhi
    ptL1T = tree.l1tPhi
    ptHLT = tree.hltPhi
    if ptHLT > 0:
        deltaPtHLT = (ptHLT - ptOffline)/ptOffline
        hltResolution.Fill(deltaPtHLT)
    if ptL1T > 0:
        deltaPtL1T = (ptL1T - ptOffline)/ptOffline
        l1tResolution.Fill(deltaPtL1T)

l1tCanvas = TCanvas("l1tCanvas")
l1tResolution.GetXaxis().SetTitle("\\frac{ptL1T - ptOffline}{ptOffline}")
l1tResolution.GetYaxis().SetTitle("Events")
l1tResolution.Draw()
l1tCanvas.Update()
l1tCanvas.Print("l1tResolution.pdf", "pdf")

hltCanvas = TCanvas("hltCanvas")
hltResolution.GetXaxis().SetTitle("\\frac{ptHLT - ptOffline}{ptOffline}")
hltResolution.GetYaxis().SetTitle("Events")
hltResolution.Draw()
hltCanvas.Update()
hltCanvas.Print("hltResolution.pdf", "pdf")

hlt_l1t_Canvas = TCanvas("hlt_l1t_Canvas")
hltResolution.Draw()
l1tResolution.Draw("SAME")
hlt_l1t_Canvas.Update()
hlt_l1t_Canvas.Print("hlt_l1t_Resolution.pdf", "pdf")

raw_input()

#!/usr/bin/python

from ROOT import *
from array import array

gStyle.SetOptStat(111111)

fIn = TFile.Open('NTuple_10Ago_Riccardo.root')
tree = fIn.Get('Ntuplizer/TagAndProbe')

l1tResolution1 = TH1I("L1TriggerReso1", "L1T resolution for offline pt between 20 and 40 GeV", 100, -100, 100)
hltResolution1 = TH1I("HLTriggerReso1", "HLT resolution for offline pt between 20 and 40 GeV", 100, -100, 100)
l1tResolution2 = TH1I("L1TriggerReso2", "L1T resolution for offline pt between 40 and 100 GeV", 100, -100, 100)
hltResolution2 = TH1I("HLTriggerReso2", "HLT resolution for offline pt between 40 and 100 GeV", 100, -100, 100)
#x = resolution
#y = R
l1tResolutionVsDR = TH2I("L1TriggerResoVsDR", "L1T resolution vs DR OS+SS", 100, -100, 100, 500, -0.1, 0.6)
hltResolutionVsDR = TH2I("HLTriggerResoVsDR", "HLT resolution vs DR OS+SS", 100, -100, 100, 50, -0.1, 0.6)

l1tResolutionVsDRSS = TH2I("L1TriggerResoVsDRSS", "L1T resolution vs DR SS", 100, -100, 100, 50, -0.1, 0.6)
hltResolutionVsDRSS = TH2I("HLTriggerResoVsDRSS", "HLT resolution vs DR SS", 100, -100, 100, 50, -0.1, 0.6)

l1tResolutionVsDROS = TH2I("L1TriggerResoVsDROS", "L1T resolution vs DR OS", 100, -100, 100, 50, -0.1, 0.6)
hltResolutionVsDROS = TH2I("HLTriggerResoVsDROS", "HLT resolution vs DR OS", 100, -100, 100, 50, -0.1, 0.6)

for iEv in range (0, tree.GetEntries()):
    tree.GetEntry(iEv)
    if tree.tauTriggerBits > 0:
        ptOffline = tree.tauPt
        ptL1T = tree.l1tPt
        ptHLT = tree.hltPt
        deltaPtL1T = ptOffline - ptL1T
        deltaPtHLT = ptOffline - ptHLT
        deltaRHLT = ((tree.tauEta - tree.hltEta)**2 + (tree.tauPhi - tree.hltPhi)**2)**0.5
        deltaRL1T = ((tree.tauEta - tree.l1tEta)**2 + (tree.tauPhi - tree.l1tPhi)**2)**0.5

        #if tree.isOS == True:
        if (ptOffline > 20)  and (ptOffline < 40):
            l1tResolution1.Fill(deltaPtL1T)
            hltResolution1.Fill(deltaPtHLT)
        if (ptOffline > 40)  and (ptOffline < 100):
            l1tResolution2.Fill(deltaPtL1T)
            hltResolution2.Fill(deltaPtHLT)


        hltResolutionVsDR.Fill(deltaPtHLT, deltaRHLT)

        l1tResolutionVsDR.Fill(deltaPtL1T, deltaRL1T)

        if tree.isOS == True :
            hltResolutionVsDROS.Fill(deltaPtHLT, deltaRHLT)
            l1tResolutionVsDROS.Fill(deltaPtL1T, deltaRL1T)
        else :
            hltResolutionVsDRSS.Fill(deltaPtHLT, deltaRHLT)
            l1tResolutionVsDRSS.Fill(deltaPtL1T, deltaRL1T)

# l1tResolution1.GetXaxis().SetTitle("ptOffline - ptL1T")
# l1tResolution1.Draw()
# c1.Update()
# c1.Print("l1tReso_40_60.pdf", "pdf")
# l1tResolution2.GetXaxis().SetTitle("ptOffline - ptL1T")
# l1tResolution2.Draw()
# c1.Update()
# c1.Print("l1tReso_60_100.pdf", "pdf")
# hltResolution1.GetXaxis().SetTitle("ptOffline - ptHLT")
# hltResolution1.Draw()
# c1.Update()
# c1.Print("hltReso_20_40.pdf", "pdf")
# hltResolution2.GetXaxis().SetTitle("ptOffline - ptHLT")
# hltResolution2.Draw()
# c1.Update()
# c1.Print("hltReso_40_100.pdf", "pdf")
#
hltResolutionVsDR.GetXaxis().SetTitle("ptOffline - ptHLT")
hltResolutionVsDR.GetYaxis().SetTitle("Delta R")
hltResolutionVsDR.Draw("COLZ")
c1.Update()
c1.Print("hltResoVsDR.pdf", "pdf")
#
# l1tResolutionVsDR.GetXaxis().SetTitle("ptOffline - ptL1T")
# l1tResolutionVsDR.GetYaxis().SetTitle("Delta R")
# l1tResolutionVsDR.Draw("COLZ")
# c1.Update()
# c1.Print("l1tResoVsDR.pdf", "pdf")
#
# hltResolutionVsDROS.GetXaxis().SetTitle("ptOffline - ptHLT")
# hltResolutionVsDROS.GetYaxis().SetTitle("Delta R")
# hltResolutionVsDROS.Draw("COLZ")
# c1.Update()
# c1.Print("hltResoVsDROS.pdf", "pdf")
#
# l1tResolutionVsDROS.GetXaxis().SetTitle("ptOffline - ptL1T")
# l1tResolutionVsDROS.GetYaxis().SetTitle("Delta R")
# l1tResolutionVsDROS.Draw("COLZ")
# c1.Update()
# c1.Print("l1tResoVsDROS.pdf", "pdf")
#
# hltResolutionVsDRSS.GetXaxis().SetTitle("ptOffline - ptHLT")
# hltResolutionVsDRSS.GetYaxis().SetTitle("Delta R")
# hltResolutionVsDRSS.Draw("COLZ")
# c1.Update()
# c1.Print("hltResoVsDRSS.pdf", "pdf")
#
# l1tResolutionVsDRSS.GetXaxis().SetTitle("ptOffline - ptL1T")
# l1tResolutionVsDRSS.GetYaxis().SetTitle("Delta R")
# l1tResolutionVsDRSS.Draw("COLZ")
# c1.Update()
# c1.Print("l1tResoVsDRSS.pdf", "pdf")

raw_input()

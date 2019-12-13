#!/usr/bin/python
#
#
# Epoch B: 273150 - 275376
# Epoch C: 275420 - 276283
# Epoch D: 276315 - 276811
#
from ROOT import *
from array import array

gStyle.SetOptStat(111111)

epochBMinRunNumber = 273150
epochBMaxRunNumber = 275376

epochCMinRunNumber = 275420
epochCMaxRunNumber = 276283

epochDMinRunNumber = 276315
epochDMaxRunNumber = 276811

fIn = TFile.Open('NTuple_10Ago_Riccardo.root')

tree = fIn.Get('Ntuplizer/TagAndProbe')

binning = [4, 8, 10, 12, 14, 16, 18, 19, 20, 21, 22, 23, 24, 26, 28, 30, 32, 35, 40, 45, 50, 60, 70, 90, 110, 140]
bins = array('d', binning)

triggerNamesTree = fIn.Get("Ntuplizer/triggerNames")

triggerNamesList = []

l1tCuts = [28, 34, 42]

# hpass = TH1F ("hpass", "hpass", 75, 0, 150)
# htot = TH1F ("htot", "htot", 75, 0, 150)
hPassListHLT_SS = []
hPassListHLT_OS = []
hTotListHLT_SS = []
hTotListHLT_OS = []
hPassListL1T_SS = []
hPassListL1T_OS = []
hTotListL1T_SS = []
hTotListL1T_OS = []
hPassListL1T_Iso_SS = []
hPassListL1T_Iso_OS = []
hTotListL1T_Iso_SS = []
hTotListL1T_Iso_OS = []

turnOnList_HLT = []
turnOnList_L1T = []
turnOnList_L1T_Iso = []

for iTrig in range (0, 6):
    triggerNamesTree.GetEntry(iTrig)
    triggerNamesList.append(triggerNamesTree.triggerNames.Data())

print "Creating histograms"

#Preparing the Histograms
for bitIndex in range(0, len(triggerNamesList)):
    hPassListHLT_OS.append(TH1F("hPassOS_"+triggerNamesList[bitIndex], "hPassOS_"+triggerNamesList[bitIndex], len(binning)-1, bins))
    hTotListHLT_OS.append(TH1F("hTotOS_"+triggerNamesList[bitIndex], "hTotOS_"+triggerNamesList[bitIndex], len(binning)-1, bins))
    hPassListHLT_SS.append(TH1F("hPassSS_"+triggerNamesList[bitIndex], "hPassSS_"+triggerNamesList[bitIndex], len(binning)-1, bins))
    hTotListHLT_SS.append(TH1F("hTotSS_"+triggerNamesList[bitIndex], "hTotSS_"+triggerNamesList[bitIndex], len(binning)-1, bins))
    turnOnList_HLT.append(TGraphAsymmErrors())

hPassTest_OS = TH1F("hPassTestOS", "hPassTestOS", len(binning)-1, bins)
hPassTest_SS = TH1F("hPassTestSS", "hPassTestSS", len(binning)-1, bins)
hTotTest_OS = TH1F("hTotTestOS", "hTotTestOS", len(binning)-1, bins)
hTotTest_SS = TH1F("hTotTestSS", "hTotTestSS", len(binning)-1, bins)
hTurnOnTest = TGraphAsymmErrors()

for cutIndex in range(0, len(l1tCuts)):
    hTotListL1T_OS.append(TH1F("hTotL1IOS_" + str(l1tCuts[cutIndex]), "hTotL1OS_"+str(l1tCuts[cutIndex]), len(binning)-1, bins))
    hPassListL1T_OS.append(TH1F("hPassL1OS_" + str(l1tCuts[cutIndex]), "hPassL1OS_"+str(l1tCuts[cutIndex]), len(binning)-1, bins))
    hTotListL1T_SS.append(TH1F("hTotL1SS_" + str(l1tCuts[cutIndex]), "hTotL1SS_"+str(l1tCuts[cutIndex]), len(binning)-1, bins))
    hPassListL1T_SS.append(TH1F("hPassL1SS_" + str(l1tCuts[cutIndex]), "hPassL1SS_"+str(l1tCuts[cutIndex]), len(binning)-1, bins))
    hTotListL1T_Iso_OS.append(TH1F("hTotL1IsoOS_" + str(l1tCuts[cutIndex]), "hTotL1IsoOS_"+str(l1tCuts[cutIndex]), len(binning)-1, bins))
    hPassListL1T_Iso_OS.append(TH1F("hPassL1IsoOS_" + str(l1tCuts[cutIndex]), "hPassL1IsoOS_"+str(l1tCuts[cutIndex]), len(binning)-1, bins))
    hTotListL1T_Iso_SS.append(TH1F("hTotL1IsoSS_" + str(l1tCuts[cutIndex]), "hTotL1IsoSS_"+str(l1tCuts[cutIndex]), len(binning)-1, bins))
    hPassListL1T_Iso_SS.append(TH1F("hPassL1IsoSS_" + str(l1tCuts[cutIndex]), "hPassL1IsoSS_"+str(l1tCuts[cutIndex]), len(binning)-1, bins))
    turnOnList_L1T.append(TGraphAsymmErrors())
    turnOnList_L1T_Iso.append(TGraphAsymmErrors())

print "Populating histograms"

#Populating the histograms
for iEv in range (0, tree.GetEntries()):
    tree.GetEntry(iEv)

    #if tree.RunNumber < epochDMinRunNumber and tree.RunNumber > epochDMaxRunNumber:
    #    continue

    if abs(tree.tauEta) > 2.1:
        continue

    pt = tree.tauPt

    #HLT Plots
    triggerBits = tree.tauTriggerBits
    for bitIndex in range(0, len(triggerNamesList)):
        if tree.isOS == True:
            hTotListHLT_OS[bitIndex].Fill(pt)
            if ((triggerBits >> bitIndex) & 1) == 1:
                hPassListHLT_OS[bitIndex].Fill(pt)
        else:
            hTotListHLT_SS[bitIndex].Fill(pt)
            if ((triggerBits >> bitIndex) & 1) == 1:
                hPassListHLT_SS[bitIndex].Fill(pt)

    #L1 Plots
    l1tPt = tree.l1tPt
    for cutIndex in range (0, len(l1tCuts)):
        if tree.isOS == True:
            hTotListL1T_OS[cutIndex].Fill(pt)
            hTotListL1T_Iso_OS[cutIndex].Fill(pt)
            if l1tPt > l1tCuts[cutIndex] :
                hPassListL1T_OS[cutIndex].Fill(pt)
                if tree.l1tIso == 1:
                    hPassListL1T_Iso_OS[cutIndex].Fill(pt)
        else :
            hTotListL1T_SS[cutIndex].Fill(pt)
            hTotListL1T_Iso_SS[cutIndex].Fill(pt)
            if l1tPt > l1tCuts[cutIndex] :
                hPassListL1T_SS[cutIndex].Fill(pt)
                if tree.l1tIso == 1:
                    hPassListL1T_Iso_SS[cutIndex].Fill(pt)

    hltPt = tree.hltPt

    if tree.isOS == True:
        hTotTestOS.Fill(pt)
        if ((triggerBits >> 5) & 1) == 1:
            if (hltPt > 35) and (l1tPt > 27.5):
                hPassTestOS.Fill(pt)
    else:
        hTotTestSS.Fill(pt)
        if ((triggerBits >> 5) & 1) == 1:
            if (hltPt > 35) and (l1tPt > 27.5):
                hPassTestSS.Fill(pt)

print "Calculating efficiencies"

#Calculating and saving the efficiencies

c1 = TCanvas ("c1", "c1", 600, 600)
c1.SetGridx()
c1.SetGridy()
fOut = TFile ("turnOn.root", "recreate")

for bitIndex in range(0, len(triggerNamesList)):

    hPassListHLT_OS[bitIndex].Add(hPassListHLT_SS[bitIndex], -1)
    hTotListHLT_OS[bitIndex].Add(hTotListHLT_SS[bitIndex], -1)

    for binIndex in range(1, hPassListHLT_OS[bitIndex].GetNbinsX() - 1):
        if hPassListHLT_OS[bitIndex].GetBinContent(binIndex) > hTotListHLT_OS[bitIndex].GetBinContent(binIndex):
            hPassListHLT_OS[bitIndex].SetBinContent(binIndex, hTotListHLT_OS[bitIndex].GetBinContent(binIndex))

    turnOnList_HLT[bitIndex].Divide(hPassListHLT_OS[bitIndex], hTotListHLT_OS[bitIndex], "cl=0.683 b(1,1) mode")
    turnOnList_HLT[bitIndex].SetMarkerStyle(8)
    turnOnList_HLT[bitIndex].SetMarkerSize(0.8)
    turnOnList_HLT[bitIndex].SetMarkerColor(kRed)
    turnOnList_HLT[bitIndex].GetXaxis().SetTitle("p_t (GeV)")
    turnOnList_HLT[bitIndex].GetYaxis().SetRangeUser(0, 1.05)
    turnOnList_HLT[bitIndex].GetYaxis().SetTitle("Efficiency")
    turnOnList_HLT[bitIndex].SetTitle(triggerNamesList[bitIndex] + " turn-on curve")
    turnOnList_HLT[bitIndex].Draw("AP")
    c1.Update()
    c1.Print("turnOn_" + triggerNamesList[bitIndex] + ".png", "png")
    hTurnOn = hPassListHLT_OS[bitIndex].Clone("hTurnOn_" + triggerNamesList[bitIndex])
    hTurnOn.Divide(hTotListHLT_OS[bitIndex])
    hTurnOn.Write()
    hPassListHLT_OS[bitIndex].Write()
    hTotListHLT_OS[bitIndex].Write()

for cutIndex in range(0, len(l1tCuts)):

    hPassListL1T_OS[cutIndex].Add(hPassListL1T_SS[cutIndex], -1)
    hTotListL1T_OS[cutIndex].Add(hTotListL1T_SS[cutIndex], -1)

    for binIndex in range(1, hPassListL1T_OS[cutIndex].GetNbinsX() - 1):
        if hPassListL1T_OS[cutIndex].GetBinContent(binIndex) > hTotListL1T_OS[cutIndex].GetBinContent(binIndex):
            hPassListL1T_OS[cutIndex].SetBinContent(binIndex, hTotListL1T_OS[cutIndex].GetBinContent(binIndex))
        if hPassListL1T_Iso_OS[cutIndex].GetBinContent(binIndex) > hTotListL1T_Iso_OS[cutIndex].GetBinContent(binIndex):
            hPassListL1T_Iso_OS[cutIndex].SetBinContent(binIndex, hTotListL1T_Iso_OS[cutIndex].GetBinContent(binIndex))
    turnOnList_L1T[cutIndex].Divide(hPassListL1T_OS[cutIndex], hTotListL1T_OS[cutIndex], "cl=0.683 b(1,1) mode")
    turnOnList_L1T[cutIndex].SetMarkerStyle(8)
    turnOnList_L1T[cutIndex].SetMarkerSize(0.8)
    turnOnList_L1T[cutIndex].SetMarkerColor(kRed)
    turnOnList_L1T[cutIndex].GetXaxis().SetTitle("p_t (GeV)");
    turnOnList_L1T[cutIndex].GetYaxis().SetTitle("Efficiency");
    turnOnList_L1T[cutIndex].SetTitle("L1 trigger cut " + str(l1tCuts[cutIndex]) + " turn-on curve")
    turnOnList_L1T[cutIndex].Draw("AP")
    c1.Update()
    c1.Print("turnOnL1_" + str(l1tCuts[cutIndex]) + ".png", "png")
    turnOnList_L1T_Iso[cutIndex].Divide(hPassListL1T_Iso_OS[cutIndex], hTotListL1T_Iso_OS[cutIndex], "cl=0.683 b(1,1) mode")
    turnOnList_L1T_Iso[cutIndex].SetMarkerStyle(8)
    turnOnList_L1T_Iso[cutIndex].SetMarkerSize(0.8)
    turnOnList_L1T_Iso[cutIndex].SetMarkerColor(kRed)
    turnOnList_L1T_Iso[cutIndex].GetXaxis().SetTitle("p_t (GeV)");
    turnOnList_L1T_Iso[cutIndex].GetYaxis().SetTitle("Efficiency");
    turnOnList_L1T_Iso[cutIndex].SetTitle("L1 trigger iso cut " + str(l1tCuts[cutIndex]) + " turn-on curve")
    turnOnList_L1T_Iso[cutIndex].Draw("AP")
    c1.Update()
    c1.Print("turnOnL1Iso_" + str(l1tCuts[cutIndex]) + ".png", "png")
    hTurnOn = hPassListL1T_OS[cutIndex].Clone("hTurnOnL1_" + str(l1tCuts[cutIndex]))
    hTurnOn.Divide(hTotListL1T_OS[cutIndex])
    hTurnOn.Write()
    hTurnOn_Iso = hPassListL1T_Iso_OS[cutIndex].Clone("hTurnOnL1Iso_" + str(l1tCuts[cutIndex]))
    hTurnOn_Iso.Divide(hTotListL1T_Iso_OS[cutIndex])
    hTurnOn_Iso.Write()
    hPassListL1T_Iso_OS[cutIndex].Write()
    hTotListL1T_Iso_OS[cutIndex].Write()
    hPassListL1T_OS[cutIndex].Write()
    hTotListL1T_OS[cutIndex].Write()

hPassTestOS.Add(hPassTestSS, -1)
hTotTestOS.Add(hTotTestSS, -1)

for binIndex in range(1, hPassTestOS.GetNbinsX() - 1):
    if hPassTestOS.GetBinContent(binIndex) > hTotTestOS.GetBinContent(binIndex):
        hPassTestOS.SetBinContent(binIndex, hTotTestOS.GetBinContent(binIndex))

hTurnOnTest.Divide(hPassTestOS, hTotTestOS, "cl=0.683 b(1,1) mode")
hTurnOnTest.SetMarkerStyle(8)
hTurnOnTest.SetMarkerSize(0.8)
hTurnOnTest.SetMarkerColor(kRed)
hTurnOnTest.GetXaxis().SetTitle("p_t (GeV)")
hTurnOnTest.GetYaxis().SetRangeUser(0, 1.05)
hTurnOnTest.GetYaxis().SetTitle("Efficiency")
hTurnOnTest.SetTitle("test turn-on curve")
hTurnOnTest.Draw("AP")
c1.Update()
c1.Print("turnOnTest.png", "png")


raw_input()

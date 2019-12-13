#!/usr/bin/python
#
#
# Epoch B: 273150 - 275376
# Epoch C: 275420 - 276283
# Epoch D: 276315 - 276811
#
from ROOT import *
from array import array

fIn = TFile.Open('NTuple_10Ago_Riccardo.root')

tree = fIn.Get('Ntuplizer/TagAndProbe')

triggerNamesTree = fIn.Get("Ntuplizer/triggerNames")



print "Creating histograms"

hPtOS = TH1F ("hPtOS", "tauPt distribution - OS pairs", 150, 0, 150)
hPtSS = TH1F ("hPtSS", "tauPt distribution - SS pairs", 150, 0, 150)

hEtaOS = TH1F ("hEtaOS", "tauEta distribution - OS pairs", 100, -2.5, +2.5)
hEtaSS = TH1F ("hEtaSS", "tauEta distribution - SS pairs", 100, -2.5, +2.5)

print "Populating histograms"

#Populating the histograms
for iEv in range (0, tree.GetEntries()):
    tree.GetEntry(iEv)

    #if tree.RunNumber < epochDMinRunNumber and tree.RunNumber > epochDMaxRunNumber:
    #    continue

    if abs(tree.tauEta) > 2.1:
        continue

    tauPt = tree.tauPt
    tauEta = tree.tauEta
    isOS = tree.isOS

    if isOS == True:
        hPtOS.Fill(tauPt)
        hEtaOS.Fill(tauEta)

    if isOS == False:
        hPtSS.Fill(tauPt)
        hEtaSS.Fill(tauEta)

#Calculating and saving the efficiencies

c1 = TCanvas ("c1", "c1", 600, 600)
c1.SetGridx()
c1.SetGridy()

hPtOS_SS = hPtOS.Clone("hPtOS_SS")
hPtOS_SS.Add(hPtSS, -1)
hPtOS_SS.SetTitle("tauPt distribution - OS-SS pairs")
hEtaOS_SS = hEtaOS.Clone("hEtaOS_SS")
hEtaOS_SS.Add(hEtaSS, -1)
hEtaOS_SS.SetTitle("tauEta distribution - OS-SS pairs")

hPtOS.Draw()
c1.Update()
c1.Print("hPtOS.png", "png")
hPtSS.Draw()
c1.Update()
c1.Print("hPtSS.png", "png")
hPtOS_SS.Draw()
c1.Update()
c1.Print("hPtOS_SS.png", "png")

hEtaOS.Draw()
c1.Update()
c1.Print("hEtaOS.png", "png")
hEtaSS.Draw()
c1.Update()
c1.Print("hEtaSS.png", "png")
hEtaOS_SS.Draw()
c1.Update()
c1.Print("hEtaOS_SS.png", "png")

raw_input()

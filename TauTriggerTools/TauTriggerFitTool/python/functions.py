from ROOT import *
from array import array
from math import sqrt

class functions:
	def __init__(self,graphEffi, name, dm, itype, func, histo, graph, fitresult, CL):
		self.graphEffi = graphEffi
		self.name = name
		self.func = func
		self.histo = histo
		self.graph = graph
		self.dm = dm
		self.itype = itype
		self.fitresult = fitresult
		self.CL =CL

    # This function is taken from Tyler Ruggles SF tool
	# https://github.com/truggles/TauTriggerSFs2017/blob/master/python/helpers.py#L9-L40
	# Function to create TH1Fs from TGraphAsymmErrors
	# This does not preserve the asymmetric errors, only
	# bin width and value and does a rough approximation
	# on symmetric errors.
	def getTH1FfromTGraphAsymmErrors( self  ) :
                
	    # Holding vals for TH1F binning and y-vals
                xSpacing = array( 'd', [] )
                yVals = array( 'd', [] )
                yErrors = array( 'd', [] )
            
                nVals = self.graphEffi.GetN()
                x, y = Double(0.), Double(0.)
                xEPlus, xEMin = 0., 0.
                yEPlus, yEMin = 0., 0.

                for n in range( nVals ) :
                        self.graphEffi.GetPoint( n, x, y )
                        xEPlus = self.graphEffi.GetErrorXhigh( n )
                        xEMin = self.graphEffi.GetErrorXlow( n )
                        yEPlus = self.graphEffi.GetErrorYhigh( n )
                        yEMin = self.graphEffi.GetErrorYlow( n )
                        xSpacing.append( x-xEMin )
                        yVals.append( y )
                        # To simplify, take asymm errors and go to approximation
                        # of symmetric for TH1
                        yErrors.append( sqrt(yEPlus**2 + yEMin**2) )

                # Don't forget to add the high end of last bin
                xSpacing.append( x+xEPlus )
    
                outH = TH1F( self.name, self.name, len(xSpacing)-1, xSpacing )
                for bin in range( 1, outH.GetNbinsX()+1 ) :
                        outH.SetBinContent( bin, yVals[bin-1] )
                        outH.SetBinError( bin, yErrors[bin-1] )
                return outH
    	
    	
        def getScaleFactorFromFunction(self):

                SF = TGraphAsymmErrors()
                for i in range(20, 450):
                        SF.SetPoint(i, i, (self.func[0].Eval(i)/self.func[1].Eval(i)))
		
                SF.Draw("A*")
                SF.GetXaxis().SetLimits(18,600)
                SF.GetYaxis().SetRangeUser(0,1.5)
                SF.GetXaxis().SetMoreLogLabels()
                SF.SetMarkerStyle(20)
                SF.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
                SF.GetYaxis().SetTitle("SF: Data/MC")
                return SF

        def getScaleFactor(self):

                hEffi1 = TH1F("", "", 480, 20, 500 )
                hEffi2 = TH1F("", "", 480, 20, 500 )

                for i in range(0, self.histo[0].GetNbinsX()):
                        hEffi1.SetBinContent(i, self.histo[0].GetBinContent(i))
                        hEffi2.SetBinContent(i, self.histo[1].GetBinContent(i))

                hEffi1.Divide(hEffi2)
                hEffi1.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
                hEffi1.GetYaxis().SetTitle("SF: Data/MC")
                hEffi1.GetXaxis().SetRangeUser(18,600)
                hEffi1.GetXaxis().SetMoreLogLabels()
                hEffi1.GetYaxis().SetRangeUser(0,1.5)
                hEffi1.SetMarkerStyle(20)
                hEffi1.SetLineWidth(6)

                return  hEffi1
                
        def getScaleFactorError(self):

                hEffi1 = TH1F("", "", 480, 20, 500 )
                hEffi2 = TH1F("", "", 480, 20, 500 )

                for i in range(0, self.histo[0].GetNbinsX()):

                        hEffi1.SetBinContent(i, self.histo[0].GetBinContent(i))
                        hEffi1.SetBinError(i, self.histo[0].GetBinError(i))
                        hEffi2.SetBinContent(i, self.histo[1].GetBinContent(i))
                        hEffi2.SetBinError(i, self.histo[1].GetBinError(i))

                hEffi1.Divide(hEffi2)
                hEffi1.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
                hEffi1.GetYaxis().SetTitle("SF: Data/MC")
                hEffi1.GetXaxis().SetRangeUser(18,600)
                hEffi1.GetXaxis().SetMoreLogLabels()
                hEffi1.GetYaxis().SetRangeUser(0,1.5)
                hEffi1.SetMarkerStyle(20)
                hEffi1.SetLineWidth(6)

                return  hEffi1
		
                SFnew.SetMarkerStyle(20)
                SFnew.SetMarkerColor(kBlack)
                SFnew.SetLineColor(kBlue)

                mg = TMultiGraph()
                mg.Add(SF1)
                mg.Add(SF2)
                mg.Add(SF3)
                mg.Add(SF4)
                mg.Add(SF5)
                mg.Draw("AP")
                mg.GetXaxis().SetLimits(18,600)
                mg.GetYaxis().SetRangeUser(0,1.5)
                mg.GetXaxis().SetMoreLogLabels()

                mg.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
                mg.GetYaxis().SetTitle("SF: Data/MC")
                return SFnew
	
	
        def getConfidenceInterval( self): 

		TVirtualFitter.GetFitter().GetConfidenceIntervals(self.histo, self.CL)
                for i in range(0, self.graphEffi.GetN()):
                        self.graph.SetPoint(i, self.graphEffi.GetX()[i], 0)
			TVirtualFitter.GetFitter().GetConfidenceIntervals(self.graph, self.CL)

                return self.histo, self.graph

		
        def createRelativeErrors(self):

                # confidence interval#
                values = self.fitresult.GetConfidenceIntervals(0.68, False)
		interval = TGraphErrors(self.graphEffi.GetN())
                ratio = TGraphAsymmErrors(self.graphEffi.GetN())
                ratio2 = TGraphAsymmErrors(self.graphEffi.GetN())
                for i in range(0, self.graphEffi.GetN()):
                        interval.SetPoint(i, self.graphEffi.GetX()[i], self.func[self.itype].Eval(self.graphEffi.GetX()[i] ))
                        interval.SetPointError(i, 0, values[i] )
                        ratio.SetPoint(i, self.graphEffi.GetX()[i], (self.func[self.itype].Eval(self.graphEffi.GetX()[i]) - values[i])/self.func[self.itype].Eval(self.graphEffi.GetX()[i]))
                        ratio2.SetPoint(i, self.graphEffi.GetX()[i], (self.func[self.itype].Eval(self.graphEffi.GetX()[i]) + values[i])/self.func[self.itype].Eval(self.graphEffi.GetX()[i]))

                ratio.GetXaxis().SetTitleSize(0.05)
                ratio.GetYaxis().SetTitleSize(0.05)
                ratio.GetXaxis().SetTitleOffset(1.1)
                ratio.GetYaxis().SetTitleOffset(1.1)
                ratio.SetLineWidth(2)
                ratio2.SetLineWidth(2)
                ratio.Draw("A*")
                ratio2.Draw("*same")
                ratio.SetMarkerStyle(20)
                ratio2.SetMarkerStyle(20)
                ratio.GetYaxis().SetRangeUser(0.8,1.2)
                ratio.SetTitle("")
                ratio.GetYaxis().SetTitle("#frac{fit #pm error}{fit}")
                ratio.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
                ratio2.SetMarkerColor(2)
                ratio.GetXaxis().SetMoreLogLabels()

                return ratio , ratio2
                
        

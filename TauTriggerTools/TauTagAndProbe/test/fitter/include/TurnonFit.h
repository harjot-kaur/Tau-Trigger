/**
 *  @file  TurnonFit.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    05/10/2014
 *
 *  @internal
 *     Created :  05/10/2014
 * Last update :  05/10/2014 21:00:46
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef TURNONFIT_H
#define TURNONFIT_H


#include <string>
#include <vector>

#include "TFile.h"

#include "RooPlot.h"
#include "RooHist.h"
#include "RooCurve.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "FuncCB.h"
#include "FuncCB_cdf.h"

class TurnonFit
{

    public:
        TurnonFit(const std::string& name);
        ~TurnonFit();


        void fit();
        void save(TFile* outputFile);

        void setFileName(const std::string& fileName) {m_fileName = fileName;}
        void setTreeName(const std::string& treeName) {m_treeName = treeName;}
        void setXVar(const std::string& xVar, double min, double max) {m_xVar.SetName(xVar.c_str()); m_xVar.setRange(min, max);}
        void setCut(const std::string& cut) {m_cut = cut;}
        void setSelectionVars(const std::vector<std::string>& selectionVars) {m_selectionVars = selectionVars;}
        void setSelection(const std::string& selection) {m_selection = selection;}
        void setWeightVar(const std::string& weightVar) {m_weightVar = weightVar;}
        void setBinning(const std::vector<double>& binning) {m_binning = binning;
	  std::cout<<"in include:"<<std::endl;
	  for(UInt_t iBin = 0 ; iBin < m_binning.size() ; ++iBin) std::cout<<m_binning[iBin]<<std::endl;
	}
        void setCrystalBall(double max, double max0, double max1,
			    double alpha, double alpha0, double alpha1,
			    double n, double n0, double n1,
			    double mean, double mean0, double mean1,
			    double sigma, double sigma0, double sigma1,
			    double mturn, double mturn0,double mturn1,
			    double p, double p0, double p1,
			    double width, double width0, double width1);

		void setCrystalBall(double max, double max0, double max1,
			    double alpha, double alpha0, double alpha1,
			    double n, double n0, double n1,
			    double mean, double mean0, double mean1,
			    double sigma, double sigma0, double sigma1,
			    double yrise, double yrise0, double yrise1
			   );

        void setNCPU(int nCPU) {m_nCPU = nCPU;}
        void setNoFit(bool noFit) {m_noFit = noFit;}

        const std::string& name() const {return m_name;}

        void printParameters();

    private:
        std::string m_name;
        std::string m_fileName;
        std::string m_treeName;
        std::string m_cut;
        std::vector<std::string> m_selectionVars;
        std::string m_weightVar;
        std::string m_selection;
        std::vector<double> m_binning;
        int m_nCPU;
        bool m_noFit;

        RooRealVar m_xVar;
        RooRealVar m_max;
        RooRealVar m_alpha;
        RooRealVar m_n;
        RooRealVar m_mean;
        RooRealVar m_sigma;
        RooRealVar m_mturn;
		RooRealVar m_p;
		RooRealVar m_width;
	    RooRealVar m_yrise;
        
        //FuncCB* m_function;
        FuncCB_cdf* m_function;
        RooPlot* m_plot;
        RooHist* m_histo;
        RooCurve* m_fit;
        RooCurve* m_fitError1Sigma;
        RooCurve* m_fitError2Sigma;
        RooFitResult* m_fitResult;
};


#endif

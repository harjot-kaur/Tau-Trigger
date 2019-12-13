/**
 *  @file  TurnonManager.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    06/10/2014
 *
 *  @internal
 *     Created :  06/10/2014
 * Last update :  06/10/2014 15:23:14
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#include "TurnonManager.h"
#include "Utilities.h"


using namespace std;



/*****************************************************************/
TurnonManager::TurnonManager()
/*****************************************************************/
{
}


/*****************************************************************/
TurnonManager::~TurnonManager()
/*****************************************************************/
{
    m_outputFile->Close();
    for(auto turnon=m_turnonFits.begin(); turnon!=m_turnonFits.end(); turnon++)
    {
        delete *turnon;
    }
}


/*****************************************************************/
bool TurnonManager::readConfig(const std::string& config)
/*****************************************************************/
{
    int status = m_params.ReadFile(config.c_str(),EEnvLevel(0));
    if(status!=0) 
    {
        cout<<"FATAL: Cannot read configuration file "<<config<<"\n"; 
        return false;
    }
    string outputFileName = m_params.GetValue("OutputFile", "");
    m_outputFile = TFile::Open(outputFileName.c_str(), "RECREATE");

    m_nCPU = m_params.GetValue("NCPU", 1);
    m_noFit = m_params.GetValue("NoFit", false);

    int nTurnons = m_params.GetValue("Turnon.N", 1);
    for(int i=0; i<nTurnons; i++)
    {
        stringstream sName, sFile, sTree, sXVar, sCut, sSelectionVars, sSelection, sWeightVar, sBinning, sFitRange, sCBMax, sCBAlpha, sCBN, sCBMean, sCBSigma, sCBMturn, sCBP, sCBWidth, sCBYRise;
        sName         << "Turnon." << i+1 << ".Name";
        sFile         << "Turnon." << i+1 << ".File";
        sTree         << "Turnon." << i+1 << ".Tree";
        sXVar         << "Turnon." << i+1 << ".XVar";
        sCut          << "Turnon." << i+1 << ".Cut";
        sSelectionVars<< "Turnon." << i+1 << ".SelectionVars";
        sSelection    << "Turnon." << i+1 << ".Selection";
        sWeightVar    << "Turnon." << i+1 << ".WeightVar";
        sBinning      << "Turnon." << i+1 << ".Binning";
        sFitRange     << "Turnon." << i+1 << ".FitRange";
        sCBMax        << "Turnon." << i+1 << ".CB.Max";
        sCBAlpha      << "Turnon." << i+1 << ".CB.Alpha";
        sCBN          << "Turnon." << i+1 << ".CB.N";
        sCBMean       << "Turnon." << i+1 << ".CB.Mean";
        sCBSigma      << "Turnon." << i+1 << ".CB.Sigma";
        sCBMturn      << "Turnon." << i+1 << ".CB.Mturn";
        sCBP          << "Turnon." << i+1 << ".CB.P";
        sCBWidth      << "Turnon." << i+1 << ".CB.Width";
        sCBYRise      << "Turnon." << i+1 << ".CB.YRise";

        string name          = m_params.GetValue(sName         .str().c_str(), "dummy");
        string file          = m_params.GetValue(sFile         .str().c_str(), "dummy");
        string tree          = m_params.GetValue(sTree         .str().c_str(), "dummy");
        string xVar          = m_params.GetValue(sXVar         .str().c_str(), "dummy");
        string cut           = m_params.GetValue(sCut          .str().c_str(), "dummy");
        string selectionVars = m_params.GetValue(sSelectionVars.str().c_str(), "");
        string selection     = m_params.GetValue(sSelection    .str().c_str(), "");
        string weightVar     = m_params.GetValue(sWeightVar    .str().c_str(), "");
        string binning       = m_params.GetValue(sBinning      .str().c_str(), "8 10 12 14 16 18 19 20 21 22 24 26 30 35 40 45 50 60 70 100");
	cout<<"binning in TurnonManager = "<<binning<<endl;
        string fitRange      = m_params.GetValue(sFitRange     .str().c_str(), "0. 100.");
        string cbMax         = m_params.GetValue(sCBMax        .str().c_str(), "1. 0.9 1.");
        string cbAlpha       = m_params.GetValue(sCBAlpha      .str().c_str(), "3. 0.01 50.");
        string cbN           = m_params.GetValue(sCBN          .str().c_str(), "10. 1.001 50.");
        string cbMean        = m_params.GetValue(sCBMean       .str().c_str(), "20. 0. 50.");
        string cbSigma       = m_params.GetValue(sCBSigma      .str().c_str(), "2. 0.01 10.");
        string cbMturn       = m_params.GetValue(sCBMturn      .str().c_str(), "20. 10. 50.");
        string cbP           = m_params.GetValue(sCBP          .str().c_str(), "0.8 0.4 1.");
        string cbWidth       = m_params.GetValue(sCBWidth      .str().c_str(), "10. 1. 50.");
        string cbYRise       = m_params.GetValue(sCBYRise      .str().c_str(), "10. 1. 50.");

        vector<double> bins           = Utilities::stringToVector<double>(binning);
	cout<<"bins in TurnonManager = "<<endl;
	for(UInt_t iBin = 0 ; iBin < bins.size() ; ++iBin) cout<<bins[iBin]<<endl;
        vector<double> fitRangeValues = Utilities::stringToVector<double>(fitRange);
        vector<double> cbMaxValues    = Utilities::stringToVector<double>(cbMax);
        vector<double> cbAlphaValues  = Utilities::stringToVector<double>(cbAlpha);
        vector<double> cbNValues      = Utilities::stringToVector<double>(cbN);
        vector<double> cbMeanValues   = Utilities::stringToVector<double>(cbMean);
        vector<double> cbSigmaValues  = Utilities::stringToVector<double>(cbSigma);
        vector<double> cbMturnValues  = Utilities::stringToVector<double>(cbMturn);
        vector<double> cbPValues      = Utilities::stringToVector<double>(cbP);
        vector<double> cbWidthValues  = Utilities::stringToVector<double>(cbWidth);
        vector<double> cbYRiseValues  = Utilities::stringToVector<double>(cbYRise);
        vector<string> selectionVarsList;
        Utilities::tokenize(selectionVars, selectionVarsList); 



        m_turnonFits.push_back(new TurnonFit(name));
        m_turnonFits.back()->setFileName(file);
        m_turnonFits.back()->setTreeName(tree);
        m_turnonFits.back()->setNCPU(m_nCPU);
        m_turnonFits.back()->setNoFit(m_noFit);
        m_turnonFits.back()->setXVar(xVar, fitRangeValues[0], fitRangeValues[1]);
        m_turnonFits.back()->setCut(cut);
        m_turnonFits.back()->setSelectionVars(selectionVarsList);
        m_turnonFits.back()->setSelection(selection);
        m_turnonFits.back()->setWeightVar(weightVar);
        m_turnonFits.back()->setBinning(bins);
      /*
        m_turnonFits.back()->setCrystalBall(cbMaxValues[0],   cbMaxValues[1],   cbMaxValues[2],
                           cbAlphaValues[0], cbAlphaValues[1], cbAlphaValues[2],
                           cbNValues[0],     cbNValues[1],     cbNValues[2],
                           cbMeanValues[0],   cbMeanValues[1],  cbMeanValues[2],
                           cbSigmaValues[0], cbSigmaValues[1], cbSigmaValues[2],
			   cbMturnValues[0],cbMturnValues[1], cbMturnValues[2],
			   cbPValues[0], cbPValues[1], cbPValues[2],
			   cbWidthValues[0], cbWidthValues[1], cbWidthValues[2]
                           );
       */
        m_turnonFits.back()->setCrystalBall(cbMaxValues[0],   cbMaxValues[1],   cbMaxValues[2],
                           cbAlphaValues[0], cbAlphaValues[1], cbAlphaValues[2],
                           cbNValues[0],     cbNValues[1],     cbNValues[2],
                           cbMeanValues[0],   cbMeanValues[1],  cbMeanValues[2],
                           cbSigmaValues[0], cbSigmaValues[1], cbSigmaValues[2],
                           cbYRiseValues[0], cbYRiseValues[1], cbYRiseValues[2]
                          );
    }

    return true;

}


/*****************************************************************/
void TurnonManager::fit()
/*****************************************************************/
{
    for(auto turnon=m_turnonFits.begin(); turnon!=m_turnonFits.end(); turnon++)
    {
        (*turnon)->fit();
        (*turnon)->save(m_outputFile);
    }
}




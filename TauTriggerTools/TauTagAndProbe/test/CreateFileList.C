#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <iostream>
#include <map>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TROOT.h>
#include <sstream>
#include <TBranchElement.h>
#include <fstream>
#include <TROOT.h>
#include <TSystem.h>

using namespace std;

void CreateFileList(string Dataset = "/SingleMuon/Run2016H-PromptReco-v2/MINIAOD", string Outfile = "fileList.txt", string SpaceSeparatedListOfRuns = "277305 277420")
{
  std::vector<string> vectOfRuns;

  if(SpaceSeparatedListOfRuns!="")
    {
      std::stringstream ss(SpaceSeparatedListOfRuns);
      
      string i;
      
      while (ss >> i)
	{
	  vectOfRuns.push_back(i);
	  
	  if (ss.peek() == ' ')
	    ss.ignore();
	}
    }

  TString RmCommand = "rm "+Outfile;
  gSystem->Exec(RmCommand.Data());

  TString ExecQueryRuns = "python ./das_client.py --query=\"run dataset="+Dataset;
  ExecQueryRuns += "\" --limit=0 ";
  cout<<"ExecQueryRuns = "<<ExecQueryRuns<<endl;

  TString Runs = gSystem->GetFromPipe(ExecQueryRuns.Data());
  cout<<"Runs = "<<Runs<<endl;

  // for (int j=0; j< vectOfRuns.size(); j++) std::cout << vectOfRuns.at(j)<<std::endl;
      
  if(vectOfRuns.size()!=0)
    {
      for (unsigned int j=0; j< vectOfRuns.size(); j++)
	{
	  TString ExecQuery = "python ./das_client.py --query=\"file dataset="+Dataset+" run=";
	  ExecQuery += vectOfRuns.at(j);
	  ExecQuery += "\" --limit=0 >> ";
	  ExecQuery += Outfile ;
	  cout<<ExecQuery<<endl;
	  gSystem->Exec(ExecQuery.Data());
	  
	}
    }
  else
    {
      TString ExecQuery = "python ./das_client.py --query=\"file dataset="+Dataset+" ";
      ExecQuery += "\" --limit=0 >> ";
      ExecQuery += Outfile ;
      cout<<ExecQuery<<endl;
      gSystem->Exec(ExecQuery.Data());
    }

  

}

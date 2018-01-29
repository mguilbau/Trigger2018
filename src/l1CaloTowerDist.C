#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TDatime.h"

#include "include/getLinBins.h"
#include "include/L1AnalysisL1CaloTowerDataFormat.h"
#include "include/plotUtilities.h"
#include "include/histDefUtility.h"

int l1CaloTowerDist(const std::string inFileName)
{
  TDatime* date = new TDatime();
  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".root") != std::string::npos){
    outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");
    outFileName = outFileName + "_L1CaloTowerDist_" + std::to_string(date->GetDate()) + ".root";
  }

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  const Int_t nL1PtBins = 50;
  Float_t l1PtLow = 5.;
  Float_t l1PtHi = 105.;
  Double_t l1PtBins[nL1PtBins+1];
  getLinBins(l1PtLow, l1PtHi, nL1PtBins, l1PtBins);

  TH1F* l1CaloTowerPt_p = new TH1F("l1CaloTowerPt_p", ";L1 Calo. Tower p_{T} (Tot);Counts", nL1PtBins, l1PtBins);
  TH1F* l1CaloTowerPt_EM_p = new TH1F("l1CaloTowerPt_EM_p", ";L1 Calo. Tower p_{T} (EM);Counts", nL1PtBins, l1PtBins);
  TH1F* l1CaloTowerPt_HAD_p = new TH1F("l1CaloTowerPt_HAD_p", ";L1 Calo. Tower p_{T} (HAD);Counts", nL1PtBins, l1PtBins);

  const Int_t nL1PtBins2 = 4;
  const Double_t l1PtBins2Low[nL1PtBins2] = {5, 15, 30, 60};
  const Double_t l1PtBins2Hi[nL1PtBins2] = {15, 30, 60, 105};

  TH1F* l1CaloTowerEta_p[nL1PtBins2];
  TH1F* l1CaloTowerPhi_p[nL1PtBins2];
  TH2F* l1CaloTowerEtaPhi_p[nL1PtBins2];

  for(Int_t lI = 0; lI < nL1PtBins2; ++lI){
    std::string ptStr = "Pt" + prettyString(l1PtBins2Low[lI], 0, true) + "to" + prettyString(l1PtBins2Hi[lI], 0, true);
    std::string ptStr2 = prettyString(l1PtBins2Low[lI], 0, false) + " < p_{T} < " + prettyString(l1PtBins2Hi[lI], 0, false);

    l1CaloTowerEta_p[lI] = new TH1F(("l1CaloTowerEta_" + ptStr + "_p").c_str(), (";L1 Calo Tower #eta (Tot, " + ptStr2 + ");Counts").c_str(), 83, -41.5, 41.5);
    l1CaloTowerPhi_p[lI] = new TH1F(("l1CaloTowerPhi_" + ptStr + "_p").c_str(), (";L1 Calo Tower #phi (Tot, " + ptStr2 + ");Counts").c_str(), 74, -0.5, 73.5);
    l1CaloTowerEtaPhi_p[lI] = new TH2F(("l1CaloTowerEtaPhi_" + ptStr + "_p").c_str(), (";L1 Calo Tower #eta (Tot, " + ptStr2 + ");L1 Calo Tower #phi").c_str(), 83, -41.5, 41.5, 74, -0.5, 73.5);

    centerTitles({l1CaloTowerEta_p[lI], l1CaloTowerPhi_p[lI], l1CaloTowerEtaPhi_p[lI]});
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* l1CaloTree_p = (TTree*)inFile_p->Get("l1CaloTowerEmuTree/L1CaloTowerTree");
  L1Analysis::L1AnalysisL1CaloTowerDataFormat *towers_ = new L1Analysis::L1AnalysisL1CaloTowerDataFormat();
  l1CaloTree_p->SetBranchAddress("L1CaloTower", &towers_);

  const Int_t nEntries = l1CaloTree_p->GetEntries();

  std::cout << "Processing.." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;

    l1CaloTree_p->GetEntry(entry);

    for(unsigned int i = 0; i < towers_->iet.size(); ++i){
      if(towers_->iet.at(i) >= l1PtLow && towers_->iet.at(i) < l1PtHi) l1CaloTowerPt_p->Fill(towers_->iet.at(i), 1./nEntries);
      if(towers_->iem.at(i) >= l1PtLow && towers_->iem.at(i) < l1PtHi) l1CaloTowerPt_p->Fill(towers_->iem.at(i), 1./nEntries);
      if(towers_->ihad.at(i) >= l1PtLow && towers_->ihad.at(i) < l1PtHi) l1CaloTowerPt_p->Fill(towers_->ihad.at(i), 1./nEntries);

      for(Int_t lI = 0; lI < nL1PtBins2; ++lI){
	if(towers_->iet.at(i) >= l1PtBins2Low[lI] && towers_->iet.at(i) < l1PtBins2Hi[lI]){
	  l1CaloTowerEta_p[lI]->Fill(towers_->ieta.at(i));
	  l1CaloTowerPhi_p[lI]->Fill(towers_->iphi.at(i));
	  l1CaloTowerEtaPhi_p[lI]->Fill(towers_->ieta.at(i), towers_->iphi.at(i));
	}
      }

    }
  }

  std::cout << "Finished processing..." << std::endl;

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  l1CaloTowerPt_p->Write("", TObject::kOverwrite);
  l1CaloTowerPt_EM_p->Write("", TObject::kOverwrite);
  l1CaloTowerPt_HAD_p->Write("", TObject::kOverwrite);

  delete l1CaloTowerPt_p;
  delete l1CaloTowerPt_EM_p;
  delete l1CaloTowerPt_HAD_p;

  for(Int_t lI = 0; lI < nL1PtBins2; ++lI){
    l1CaloTowerEta_p[lI]->Write("", TObject::kOverwrite);
    l1CaloTowerPhi_p[lI]->Write("", TObject::kOverwrite);
    l1CaloTowerEtaPhi_p[lI]->Write("", TObject::kOverwrite);

    delete l1CaloTowerEta_p[lI];
    delete l1CaloTowerPhi_p[lI];
    delete l1CaloTowerEtaPhi_p[lI];
  }

  outFile_p->Close();
  delete outFile_p;

  delete date;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./l1CaloTowerDist.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += l1CaloTowerDist(argv[1]);
  return retVal;
}

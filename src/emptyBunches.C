//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TH1F.h"

//Local dependencies
#include "include/checkMakeDir.h"
#include "include/histDefUtility.h"
#include "include/L1AnalysisL1CaloTowerDataFormat.h"
#include "include/plotUtilities.h"
#include "include/returnRootFileContentsList.h"

int emptyBunches(const std::string inFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "inFileName \'" << inFileName << "\' is not valid. return 1" << std::endl;
    return 1;
  }
  else if(inFileName.find(".root") == std::string::npos){
    std::cout << "inFileName \'" << inFileName << "\' is not valid. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;
  
  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), 5, "");

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  
  outFileName = "output/" + dateStr + "/" + outFileName + "_EmptyBunches_" + dateStr + ".root";
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* l1IEtSumNotHF_IEta28_h = new TH1F("l1IEtSumNotHF_IEta28_h", ";#Sigma iEt (|ieta| <= 28);#frac{1}{N_{Evt}} Counts", 100, -0.5, 199.5);
  TH1F* l1IEtSumNotHF_IEta24_h = new TH1F("l1IEtSumNotHF_IEta24_h", ";#Sigma iEt (|ieta| <= 24);#frac{1}{N_{Evt}} Counts", 100, -0.5, 199.5);

  centerTitles({l1IEtSumNotHF_IEta28_h, l1IEtSumNotHF_IEta24_h});
  setSumW2({l1IEtSumNotHF_IEta28_h, l1IEtSumNotHF_IEta24_h});

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> fileList = returnRootFileContentsList(inFile_p, "TTree", "Tower");
  bool hasEmuTower = false;
  bool hasNonEmuTower = false;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    if(fileList.at(fI).find("l1CaloTowerEmuTree") != std::string::npos) hasEmuTower = true;
    if(fileList.at(fI).find("l1CaloTowerTree") != std::string::npos) hasNonEmuTower = true;
  }

  TTree* l1CaloTree_p = NULL;
  if(hasEmuTower) l1CaloTree_p = (TTree*)inFile_p->Get("l1CaloTowerEmuTree/L1CaloTowerTree");
  else if(hasNonEmuTower) l1CaloTree_p = (TTree*)inFile_p->Get("l1CaloTowerTree/L1CaloTowerTree");
  else{
    std::cout << "Has no calo tower tree in \'" << inFileName << "\'. return 1" << std::endl;

    inFile_p->Close();
    delete inFile_p;

    delete l1IEtSumNotHF_IEta28_h;
    delete l1IEtSumNotHF_IEta24_h;

    outFile_p->Close();
    delete outFile_p;

    return 1;
  }

  L1Analysis::L1AnalysisL1CaloTowerDataFormat *towers_ = new L1Analysis::L1AnalysisL1CaloTowerDataFormat();
  l1CaloTree_p->SetBranchStatus("*", 0);
  l1CaloTree_p->SetBranchStatus("L1CaloTower", 1);
  l1CaloTree_p->SetBranchStatus("iet", 1);
  l1CaloTree_p->SetBranchStatus("ieta", 1);

  l1CaloTree_p->SetBranchAddress("L1CaloTower", &(towers_));

  const Int_t nEntries = l1CaloTree_p->GetEntries();
  const Int_t nDiv = TMath::Max(1, nEntries/20);

  std::cout << "Processing " << nEntries << "..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;

    l1CaloTree_p->GetEntry(entry);

    Int_t ietSumNotHF_IEta28_ = 0;
    Int_t ietSumNotHF_IEta24_ = 0;

    for(unsigned int i = 0; i < towers_->iet.size(); ++i){
      if(TMath::Abs(towers_->ieta[i]) >= 29) continue;
      ietSumNotHF_IEta28_ += towers_->iet[i];


      if(TMath::Abs(towers_->ieta[i]) >= 25) continue;      
      ietSumNotHF_IEta24_ += towers_->iet[i];
    }

    l1IEtSumNotHF_IEta28_h->Fill(ietSumNotHF_IEta28_);
    l1IEtSumNotHF_IEta24_h->Fill(ietSumNotHF_IEta24_);
  }

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  l1IEtSumNotHF_IEta28_h->Scale(1./(Double_t)nEntries);
  l1IEtSumNotHF_IEta24_h->Scale(1./(Double_t)nEntries);

  l1IEtSumNotHF_IEta28_h->Write("", TObject::kOverwrite);
  delete l1IEtSumNotHF_IEta28_h;

  l1IEtSumNotHF_IEta24_h->Write("", TObject::kOverwrite);
  delete l1IEtSumNotHF_IEta24_h;

  delete outFile_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/emptyBunches.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += emptyBunches(argv[1]);
  return retVal;
}
  

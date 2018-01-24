#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

int l1FiringFraction(const std::string inFileName)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* l1Tree_p = (TTree*)inFile_p->Get("l1UpgradeEmuTree/L1UpgradeTree");
  
  std::vector<float>* jetEt_p=0;
  std::vector<float>* jetEta_p=0;
  std::vector<float>* jetPhi_p=0;
  
  l1Tree_p->SetBranchStatus("*", 0);
  l1Tree_p->SetBranchStatus("jetEt", 1);
  l1Tree_p->SetBranchStatus("jetEta", 1);
  l1Tree_p->SetBranchStatus("jetPhi", 1);

  l1Tree_p->SetBranchAddress("jetEt", &jetEt_p);
  l1Tree_p->SetBranchAddress("jetEta", &jetEta_p);
  l1Tree_p->SetBranchAddress("jetPhi", &jetPhi_p);

  const Int_t nEntries = l1Tree_p->GetEntries();

  const Int_t nJetThresholds = 10;
  Double_t jetThresholds[nJetThresholds] = {4.0, 8.0, 16.0, 24.0, 32.0, 40.0, 48.0, 56.0, 64.0, 72.0};
  Int_t jetCounts[nJetThresholds];
  for(Int_t jI = 0; jI < nJetThresholds; ++jI){
    jetCounts[jI] = 0;
  }

  for(Int_t entry = 0; entry < nEntries; ++entry){
    l1Tree_p->GetEntry(entry);

    Float_t tempLeadingJtPt_ = -999.;

    for(unsigned int jI = 0; jI < jetEt_p->size(); ++jI){
      if(jetEt_p->at(jI) > tempLeadingJtPt_) tempLeadingJtPt_ = jetEt_p->at(jI);
    }
    
    for(Int_t jI = 0; jI < nJetThresholds; ++jI){
      if(tempLeadingJtPt_ > jetThresholds[jI]) jetCounts[jI] += 1;
      else break;
    }
  }

  for(Int_t jI = 0; jI < nJetThresholds; ++jI){
    std::cout << "Thresh " << jetThresholds[jI] << ": " << jetCounts[jI] << "/" << nEntries << std::endl;
  }

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./l1FiringFraction.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += l1FiringFraction(argv[1]);
  return retVal;
}

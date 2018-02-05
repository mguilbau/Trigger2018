#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "include/L1AnalysisL1UpgradeDataFormat.h"
#include "include/mntToXRootdFileString.h"

int l1FiringFraction(const std::string inFileName)
{
  TFile* inFile_p = NULL;
  inFile_p = TFile::Open(mntToXRootdFileString(inFileName).c_str(), "READ");
  
  if(inFile_p == NULL){
    std::cout << "WARNING: File \'" << inFileName << "\' returns NULL. return 1" << std::endl;
    return 1;
  }
  else if(inFile_p->IsZombie()){
    std::cout << "WARNING: File \'" << inFileName << "\' returns IsZombie. return 1" << std::endl;
    return 1;
  }

  TTree* l1Tree_p = (TTree*)inFile_p->Get("l1UpgradeEmuTree/L1UpgradeTree");

  L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade = new L1Analysis::L1AnalysisL1UpgradeDataFormat();

  l1Tree_p->SetBranchStatus("*", 0);
  l1Tree_p->SetBranchStatus("L1Upgrade", 1);
  l1Tree_p->SetBranchStatus("jetEt", 1);

  l1Tree_p->SetBranchAddress("L1Upgrade", &(upgrade));

  const Int_t nEntries = l1Tree_p->GetEntries();

  const Int_t nJetThresholds = 10;
  Double_t jetThresholds[nJetThresholds] = {4.0, 8.0, 16.0, 24.0, 32.0, 40.0, 48.0, 56.0, 64.0, 72.0};
  Int_t jetCounts[nJetThresholds];
  for(Int_t jI = 0; jI < nJetThresholds; ++jI){
    jetCounts[jI] = 0;
  }


  std::cout << "Processing..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    l1Tree_p->GetEntry(entry);

    Float_t tempLeadingJtPt_ = -999.;

    for(unsigned int jI = 0; jI < upgrade->jetEt.size(); ++jI){
      if(upgrade->jetEt.at(jI) > tempLeadingJtPt_) tempLeadingJtPt_ = upgrade->jetEt.at(jI);
    }
    
    for(Int_t jI = 0; jI < nJetThresholds; ++jI){
      if(tempLeadingJtPt_ > jetThresholds[jI]) jetCounts[jI] += 1;
      else break;
    }
  }

  for(Int_t jI = 0; jI < nJetThresholds; ++jI){
    std::cout << "Thresh " << jetThresholds[jI] << ": " << jetCounts[jI] << "/" << nEntries << " = " << (Double_t(jetCounts[jI]))/(Double_t(nEntries))<< std::endl;
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

//cpp dependencies
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"

//local dependencies
#include "include/checkMakeDir.h"
#include "include/L1AnalysisL1UpgradeDataFormat.h"

int l1ToTTree(const std::string inFileName)
{
  if(!checkFile(inFileName) || inFileName.find(".txt") == std::string::npos){
    std::cout << "inFileName \'" << inFileName << "\' is not valid. return 1" << std::endl;
    return 1;
  }

  std::vector<std::string> fileList;
  std::ifstream file(inFileName.c_str());
  std::string tempStr;
  while(std::getline(file, tempStr)){
    if(checkFile(tempStr) && tempStr.find(".root") != std::string::npos) fileList.push_back(tempStr);
  }
  file.close();
  
  if(fileList.size() == 0){
    std::cout << "inFileName \'" << inFileName << "\' contains no root files. return 1" << std::endl;
    return 1;
  }
  
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  outFileName.replace(outFileName.find(".txt"),4,"");
  outFileName = "output/" + dateStr + "/" + outFileName + "_" + dateStr + ".root";

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  Int_t nL1Jt_;
  std::vector<short>* l1JtIEt_p = new std::vector<short>;
  std::vector<short>* l1JtIEta_p = new std::vector<short>;
  std::vector<short>* l1JtIPhi_p = new std::vector<short>;

  TTree* l1BasicTree_p = new TTree("l1BasicTree", "");
  l1BasicTree_p->Branch("nL1Jt", &nL1Jt_, "nL1Jt/I");
  l1BasicTree_p->Branch("l1JtIEt", &l1JtIEt_p);
  l1BasicTree_p->Branch("l1JtIEta", &l1JtIEta_p);
  l1BasicTree_p->Branch("l1JtIPhi", &l1JtIPhi_p);

  
  std::cout << "Processing " << fileList.size() << " files..." << std::endl;
  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << " File " << fI << "/" << fileList.size() << std::endl;
    TFile* inFile_p = new TFile(fileList.at(fI).c_str(), "READ");
    L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
    TTree* l1Tree_p = (TTree*)inFile_p->Get("l1UpgradeEmuTree/L1UpgradeTree");

    l1Tree_p->SetBranchStatus("*", 0);
    l1Tree_p->SetBranchStatus("L1Upgrade", 1);
    l1Tree_p->SetBranchStatus("jetIEt", 1);
    l1Tree_p->SetBranchStatus("jetIEta", 1);
    l1Tree_p->SetBranchStatus("jetIPhi", 1);

    l1Tree_p->SetBranchAddress("L1Upgrade", &(upgrade));



    const Int_t nEntries = l1Tree_p->GetEntries();
    std::cout << " Processing " << nEntries << " entries..." << std::endl;
    Int_t printInt = TMath::Max(1, nEntries/20);

    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(entry%printInt == 0) std::cout << "  Entry " << entry << "/" << nEntries << std::endl;
      l1Tree_p->GetEntry(entry);

      nL1Jt_ = 0;
      l1JtIEt_p->clear();
      l1JtIEta_p->clear();
      l1JtIPhi_p->clear();

      for(unsigned int jI = 0; jI < upgrade->jetIEt.size(); ++jI){
	l1JtIEt_p->push_back(upgrade->jetIEt.at(jI));
	l1JtIPhi_p->push_back(upgrade->jetIPhi.at(jI));
	l1JtIEta_p->push_back(upgrade->jetIEta.at(jI));
	++nL1Jt_;
      }
      
      l1BasicTree_p->Fill();
    }

    delete upgrade;

    inFile_p->Close();
    delete inFile_p;
  }


  outFile_p->cd();
  
  l1BasicTree_p->Write("", TObject::kOverwrite);
  delete l1BasicTree_p;

  outFile_p->Close();
  delete outFile_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/l1ToTTree.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += l1ToTTree(argv[1]);
  return retVal;
}

//cpp dependencies
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TDatime.h"

//Local dependencies
#include "include/checkMakeDir.h"
#include "include/returnRootFileContentsList.h"

int hltFiringFraction(const std::string inFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "Warning: Given inFileName \'" << inFileName << "\' is not a valid file. return 1" << std::endl;
    return 1;
  }
  else if(inFileName.find(".root") == std::string::npos){
    std::cout << "Warning: Given inFileName \'" << inFileName << "\' is not a valid .root file. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");

  //Pick TTree
  std::vector<std::string> listOfHLTTrees = returnRootFileContentsList(inFile_p, "TTree", "hlt");
  if(listOfHLTTrees.size() == 0){
    std::cout << "Given inFileName \'" << inFileName << "\' contains no TTree w/ \'hlt\' in name. return 1" << std::endl;
    return 1;
  }
  else if(listOfHLTTrees.size() > 1){
    std::cout << "Warning: Given inFileName \'" << inFileName << "\' contains multiple TTree w/ \'hlt\' in name. Picking first, \'" << listOfHLTTrees.at(0) << "\'. Please check file that this is right choice" << std::endl;
  }
  TTree* hltTree_p = (TTree*)inFile_p->Get(listOfHLTTrees.at(0).c_str());

  //Pick Branches
  std::vector<std::string> finalListOfBranches;
  TObjArray* initListOfBranches = hltTree_p->GetListOfBranches();
  for(Int_t bI = 0; bI < initListOfBranches->GetEntries(); ++bI){
    std::string branchName = initListOfBranches->At(bI)->GetName();

    if(branchName.find("Prescl") != std::string::npos) continue;
    if(branchName.size() < 4) continue;
    if(branchName.substr(0,4).find("HLT_") == std::string::npos) continue;

    finalListOfBranches.push_back(branchName);
  }

  const Int_t nBranches = finalListOfBranches.size();
  Int_t branchVal_[nBranches];
  Int_t branchFires_[nBranches];
  
  hltTree_p->SetBranchStatus("*", 0);

  for(Int_t bI = 0; bI < nBranches; ++bI){
    hltTree_p->SetBranchStatus(finalListOfBranches.at(bI).c_str(), 1);

    hltTree_p->SetBranchAddress(finalListOfBranches.at(bI).c_str(), &(branchVal_[bI]));

    branchFires_[bI] = 0;
  }

  const Int_t nEntries = hltTree_p->GetEntries();
  const Int_t nDiv = TMath::Max(1, nEntries/20);

  std::cout << "Processing " << nEntries << "..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    hltTree_p->GetEntry(entry);

    for(Int_t bI = 0; bI < nBranches; ++bI){
      if(branchVal_[bI] == 1) ++(branchFires_[bI]);
    }
  }

  inFile_p->Close();
  delete inFile_p;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);

  std::string outFileName = inFileName.substr(0, inFileName.find(".root"));
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  outFileName = "output/" + dateStr + "/" + outFileName + "_FiringFrac_" + dateStr + ".csv";
  
  std::ofstream outFile(outFileName.c_str());

  outFile << "TriggerName,Fires (Out of " << nEntries << "),Fraction,Rate at 40kHz (Hz)," << std::endl;
  for(Int_t bI = 0; bI < nBranches; ++bI){
    Double_t fraction = branchFires_[bI];
    fraction /= (Double_t)nEntries;

    outFile << finalListOfBranches.at(bI) << "," << branchFires_[bI] << "," << fraction << "," << fraction*40000. << "," << std::endl;
  }

  outFile.close();

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/hltFiringFraction.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += hltFiringFraction(argv[1]);
  return retVal;
}

//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"

//Local dependencies
#include "include/checkMakeDir.h"
#include "include/returnRootFileContentsList.h"

int listHLT(const std::string inFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "Warning: Given inFileName \'" << inFileName << "\' is not a valid file. return 1" << std::endl;
    return 1;
  }
  else if(inFileName.find(".root") == std::string::npos){
    std::cout << "Warning: Given inFileName \'" << inFileName << "\' is not a valid .root file. return 1" << std::endl;
    return 1;
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");

  //Pick TTree
  std::vector<std::string> listOfHLTTrees = returnRootFileContentsList(inFile_p, "TTree", "hlt");
  if(listOfHLTTrees.size() == 0){
    std::cout << "Given inFileName \'" << inFileName << "\' contains no TTree w/ \'hlt\' in name. return 1" << std::endl;
    inFile_p->Close();
    delete inFile_p;
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

  inFile_p->Close();
  delete inFile_p;

  std::cout << "Listing " << finalListOfBranches.size() << " HLT triggers..." << std::endl;
  for(unsigned int bI = 0; bI < finalListOfBranches.size(); ++bI){
    std::cout << " " << bI << "/" << finalListOfBranches.size() << ": " << finalListOfBranches.at(bI) << std::endl;
  }

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/listHLT.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += listHLT(argv[1]);
  return retVal;
}

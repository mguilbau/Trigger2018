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
#include "include/stringUtil.h"

int checkHLTAgainstL1Xml(const std::string inFileName, const std::string l1XmlFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "Warning: Given inFileName \'" << inFileName << "\' is not a valid file. return 1" << std::endl;
    return 1;
  }
  else if(inFileName.find(".root") == std::string::npos){
    std::cout << "Warning: Given inFileName \'" << inFileName << "\' is not a valid .root file. return 1" << std::endl;
    return 1;
  }

  if(!checkFile(l1XmlFileName)){
    std::cout << "Warning: Given l1XmlFileName \'" << l1XmlFileName << "\' is not a valid file. return 1" << std::endl;
    return 1;
  }
  else if(l1XmlFileName.find(".xml") == std::string::npos){
    std::cout << "Warning: Given l1XmlFileName \'" << l1XmlFileName << "\' is not a valid .root file. return 1" << std::endl;
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

    branchName.replace(0, 4, "");
    branchName.replace(branchName.rfind("_v"), branchName.size(), "");

    while(branchName.find("_") != std::string::npos){branchName.replace(branchName.find("_"), 1, "");} 

    finalListOfBranches.push_back(branchName);
  }

  inFile_p->Close();
  delete inFile_p;

  std::vector<std::string> listOfL1Branches;
  std::vector<bool> listOfL1BranchesFound;
  std::ifstream inFile(l1XmlFileName.c_str());
  std::string tempStr;
  while(std::getline(inFile, tempStr)){
    if(tempStr.find("<name>L1_") == std::string::npos) continue;

    while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, "");}
    while(tempStr.find("_") != std::string::npos){tempStr.replace(tempStr.find("_"), 1, "");}
    tempStr.replace(tempStr.find("<name>"), std::string("<name>").size(), "");
    tempStr.replace(tempStr.find("</name>"), std::string("</name>").size(), "");
    
    listOfL1Branches.push_back(tempStr);
    listOfL1BranchesFound.push_back(false);
  }

  inFile.close();


  //Check both vectors for duplicates
  for(unsigned int bI = 0; bI < finalListOfBranches.size(); ++bI){
    std::string string1 = finalListOfBranches.at(bI);
    for(unsigned int bI2 = bI+1; bI2 < finalListOfBranches.size(); ++bI2){
      if(isStrSame(string1, finalListOfBranches.at(bI2))){
	std::cout << "Warning \'" << string1 << "\' is found twice in hlt vector. check and return 1" << std::endl;
	return 1;
      }
    }    
  } 

  for(unsigned int bI = 0; bI < listOfL1Branches.size(); ++bI){
    std::string string1 = listOfL1Branches.at(bI);
    for(unsigned int bI2 = bI+1; bI2 < listOfL1Branches.size(); ++bI2){
      if(isStrSame(string1, listOfL1Branches.at(bI2))){
	std::cout << "Warning \'" << string1 << "\' is found twice in hlt vector. check and return 1" << std::endl;
	return 1;
      }
    }    
  } 

  int missedHLT = 0;
  std::cout << "HLT Branches w/o corresponding L1" << std::endl;
  //Now do matching
  for(unsigned int bI = 0; bI < finalListOfBranches.size(); ++bI){
    bool isFound = false;
    std::string string1 = finalListOfBranches.at(bI);
    
    for(unsigned int bI2 = 0; bI2 < listOfL1Branches.size(); ++bI2){
      if(listOfL1BranchesFound.at(bI2)) continue;

      if(isStrSame(string1, listOfL1Branches.at(bI2))){
	isFound = true;
	listOfL1BranchesFound.at(bI2) = true;
	break;
      }
    }

    if(!isFound){
      std::cout << " " << missedHLT << ": " << finalListOfBranches.at(bI) << std::endl;
      ++missedHLT;
    }
  }

  int missedL1 = 0;
  std::cout << "L1 Branches w/o corresponding HLT" << std::endl;
  for(unsigned int bI = 0; bI < listOfL1Branches.size(); ++bI){
    if(listOfL1BranchesFound.at(bI)) continue;

    std::cout << " " << missedL1 << ": " << listOfL1Branches.at(bI) << std::endl;
    ++missedL1;
  }

  std::cout << "Total missedHLT: " << missedHLT << std::endl;
  std::cout << "Total missedL1: " << missedL1 << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/checkHLTAgainstL1Xml.exe <inFileName> <l1XmlFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += checkHLTAgainstL1Xml(argv[1], argv[2]);
  return retVal;
}

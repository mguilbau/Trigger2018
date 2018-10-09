//cpp dependencies
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <fstream>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TMath.h"

//Local dependencies
#include "include/checkMakeDir.h"
#include "include/listOfPrimes.h"
#include "include/returnRootFileContentsList.h"
#include "include/stringUtil.h"

int doPrescaling(const std::string inFileName, const std::string prescaleConfigName)
{
  if(!checkFile(inFileName)){
    std::cout << "Warning: Given inFileName \'" << inFileName << "\' is not a valid file. return 1" << std::endl;
    return 1;
  }
  else if(inFileName.find(".root") == std::string::npos){
    std::cout << "Warning: Given inFileName \'" << inFileName << "\' is not a valid .root file. return 1" << std::endl;
    return 1;
  }

  if(!checkFile(prescaleConfigName)){
    std::cout << "Warning: Given prescaleConfigName \'" << prescaleConfigName << "\' is not a valid file. return 1" << std::endl;
    return 1;
  }
  else if(prescaleConfigName.find(".txt") == std::string::npos){
    std::cout << "Warning: Given prescaleConfigName \'" << prescaleConfigName << "\' is not a valid .txt file. return 1" << std::endl;
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


  std::vector<std::string> trigNames;
  std::vector<int> trigPrescale;
  std::vector<std::string> trigPD;
  std::vector<std::string> trigSubPD;

  std::map<std::string, bool> pdMapIsFired;
  std::map<std::string, int> pdMapToFires;

  std::map<std::string, bool> subPDMapIsFired;
  std::map<std::string, int> subPDMapToFires;

  std::ifstream prescaleConfig(prescaleConfigName.c_str());
  std::string tempStr;
  while(std::getline(prescaleConfig, tempStr)){
    if(tempStr.size() == 0) continue;

    std::string trigName = tempStr.substr(0, tempStr.find(","));
    tempStr.replace(0, tempStr.find(",")+1, "");
    int prescale = std::stoi(tempStr.substr(0, tempStr.find(",")));
    tempStr.replace(0, tempStr.find(",")+1, "");
    std::string pd = tempStr.substr(0, tempStr.find(","));  
    tempStr.replace(0, tempStr.find(",")+1, "");
    std::string subPD = tempStr.substr(0, tempStr.find(","));  

    if(prescale < 0) continue;

    bool isInFile = false;
    for(unsigned int bI = 0; bI < finalListOfBranches.size(); ++bI){
      if(isStrSame(finalListOfBranches.at(bI), trigName)){
	isInFile = true;
	break;
      }
    }

    if(!isInFile){
      std::cout << "Warning: trigger \'" << trigName << "\' requested from file \'" << prescaleConfigName << "\' is not found in file \'" << inFileName << "\'. skipping" << std::endl;
      continue;
    }

    trigNames.push_back(trigName);
    trigPrescale.push_back(getNearestPrime(prescale));
    trigPD.push_back(pd);
    trigSubPD.push_back(subPD);

    if(pdMapToFires.count(pd) == 0) pdMapToFires[pd] = 0;
    if(subPDMapToFires.count(subPD) == 0) subPDMapToFires[subPD] = 0;
  }
  prescaleConfig.close();

  const Int_t nTrig = trigNames.size();
  Int_t trigVal[nTrig];
  Int_t trigFires[nTrig];
  Int_t trigPrescaledFires[nTrig];

  hltTree_p->SetBranchStatus("*", 0);
  for(Int_t bI = 0; bI < nTrig; ++bI){
    hltTree_p->SetBranchStatus(trigNames.at(bI).c_str(), 1);
    hltTree_p->SetBranchAddress(trigNames.at(bI).c_str(), &(trigVal[bI]));

    trigFires[bI] = 0;
    trigPrescaledFires[bI] = 0;
  }

  const Int_t nEntries = hltTree_p->GetEntries();
  const Int_t nDiv = TMath::Max(1, nEntries/20);

  Int_t totalFires = 0;

  std::cout << "Processing " << nEntries << "..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;

    hltTree_p->GetEntry(entry);

    for(auto const &iter : pdMapIsFired){
      pdMapIsFired[iter.first] = false;
    }

    for(auto const &iter : subPDMapIsFired){
      subPDMapIsFired[iter.first] = false;
    }
    
    bool doesGlobalFire = false;
    
    for(Int_t bI = 0; bI < nTrig; ++bI){
      if(trigVal[bI] == 1){
	if(trigFires[bI]%trigPrescale.at(bI) == 0){
	  pdMapIsFired[trigPD[bI]] = true;
	  subPDMapIsFired[trigSubPD[bI]] = true;
	  doesGlobalFire = true;
	  ++trigPrescaledFires[bI];
	}
	++trigFires[bI];
      }
    }
    
    for(auto const &iter : pdMapIsFired){
      if(iter.second) ++(pdMapToFires[iter.first]);
    }

    for(auto const &iter : subPDMapIsFired){
      if(iter.second) ++(subPDMapToFires[iter.first]);
    }

    if(doesGlobalFire) ++totalFires;
  }

  inFile_p->Close();
  delete inFile_p;

  
  std::cout << "#: Name, PD, SubPD, Final prescale, Fires, Prescaled Fires, Rate at 40kHz (Hz), Prescaled Rate at 40kHz (Hz)" << std::endl;
  for(Int_t bI = 0; bI < nTrig; ++bI){
    std::cout << " " << bI << "/" << nTrig << ": " << trigNames[bI] << ", " << trigPD[bI] << ", " << trigSubPD[bI] << ", " << trigPrescale[bI] << ", " << trigFires[bI] << ", " << trigPrescaledFires[bI] << ", " << ((Double_t)trigFires[bI])*40000./((Double_t)nEntries) << ", " << ((Double_t)trigPrescaledFires[bI])*40000./((Double_t)nEntries) << std::endl;
  }

  std::cout << std::endl;

  std::cout << "PD, Total Prescaled Fires, Rate at 40kHz (Hz)" << std::endl;
  for(auto const &iter : pdMapToFires){
    std::cout << iter.first << ", " << iter.second << ", " << ((Double_t)iter.second)*40000./((Double_t)nEntries) << std::endl;
  }

  std::cout << std::endl;

  std::cout << "SubPD, Total Prescaled Fires, Rate at 40kHz (Hz)" << std::endl;
  for(auto const &iter : subPDMapToFires){
    std::cout << iter.first << ", " << iter.second << ", " << ((Double_t)iter.second)*40000./((Double_t)nEntries) << std::endl;
  }

  std::cout << std::endl;

  std::cout << "Total fires, Rate at 40 kHz (Hz): " << totalFires << ", " << ((Double_t)totalFires)*40000./((Double_t)nEntries) << std::endl;

  std::cout << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/doPrescaling.exe <inFileName> <prescaleConfigName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += doPrescaling(argv[1], argv[2]);
  return retVal;
}

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
#include "include/listOfPrimes.h"
#include "include/returnRootFileContentsList.h"
#include "include/stringUtil.h"

int hltTurnOn(const std::string inFileName, const std::string prescaleConfigName)
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
  
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);

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

  std::vector<std::string> listOfAllTrees = returnRootFileContentsList(inFile_p, "TTree");

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
  
  std::vector<std::string> uniqueSubPD;

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

    bool isUniqueSubPD = true;
    for(unsigned int pI = 0; pI < uniqueSubPD.size(); ++pI){
      if(isStrSame(uniqueSubPD.at(pI), subPD)){
	isUniqueSubPD = false;
	break;
      }
    }
    if(isUniqueSubPD) uniqueSubPD.push_back(subPD);
  }
  prescaleConfig.close();

  const Int_t nTrig = trigNames.size();
  Int_t trigVal[nTrig];

  const Int_t nSubPD = uniqueSubPD.size();
  std::string subPDOfflieObj[nSubPD];
  
  //Editing here 
  PDDoubleEG;
    PDDoubleJet;
    PDDoubleJet30100;
    PDDoubleJet50100;
    PDDoubleMuon;
    PDMisc;
    PDSingleEG;
    PDSingleJet;
    PDSingleJet30100;
    PDSingleJet50100;
    PDSingleJetFWD;
    PDSingleJetFWD30100;
    PDSingleJetFWD50100;
    PDSingleMuon;
    PDXTrigger;

  for(Int_t sI = 0; sI < nSubPD; ++sI){
    if(isStrSame(uniqueSubPD.at(sI), )){

    }
  }
  

  hltTree_p->SetBranchStatus("*", 0);
  for(Int_t bI = 0; bI < nTrig; ++bI){
    hltTree_p->SetBranchStatus(trigNames.at(bI).c_str(), 1);
    hltTree_p->SetBranchAddress(trigNames.at(bI).c_str(), &(trigVal[bI]));
  }

  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".") != std::string::npos) outFileName.replace(outFileName.find("."), outFileName.size(), "");
  outFileName = "output/" + dateStr + "/" + outFileName + "_HLTTurnOn_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D* denom_p[nSubPD];
  TH1D* num_p[nTrig];

  //Init denom histograms
  for(Int_t dI = 0; dI < nSubPD; ++dI){
    denom_p[dI] = new TH1D(("denom_" + uniqueSubPD.at(dI) + "_h").c_str(), ";p_{T};Counts", 18, 20, 100);
  }

  for(Int_t tI = 0; tI < nSubPD; ++tI){
    std::string tempTrigName = trigNames.at(tI).substr(0, trigNames.rfind("_v"));
    tempTrigName.replace(0,4,"");
    num_p[tI] = new TH1D(("num_" + tempTrigName + "_h").c_str(), ";p_{T};Counts", 18, 20, 100);
  }

  const Int_t nEntries = hltTree_p->GetEntries();
  const Int_t nDiv = TMath::Max(1, nEntries/20);

  Int_t totalFires = 0;

  std::cout << "Processing " << nEntries << "..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;

    hltTree_p->GetEntry(entry);


  }

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();
  for(Int_t dI = 0; dI < nSubPD; ++dI){
    denom_p[dI]->Write("", TObject::kOverwrite);
    delete denom_p[dI];
  }

  for(Int_t tI = 0; tI < nSubPD; ++tI){
    num_p[tI]->Write("", TObject::kOverwrite);
    delete num_p[tI];
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/hltTurnOn.exe <inFileName> <prescaleConfigName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += hltTurnOn(argv[1], argv[2]);
  return retVal;
}

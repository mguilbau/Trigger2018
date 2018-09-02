//cpp dependencies
#include <fstream>
#include <iostream>
#include <map>
#include <string>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"

//Local dependcies
#include "include/checkMakeDir.h"
#include "include/L1AnalysisL1CaloTowerDataFormat.h"

int l1OfflineEvtDisp(const std::string inFileName, const int getEntry, const int etaLow = -41, const int etaHi = 41)
{
  if(!checkFile(inFileName) || inFileName.find(".root") == std::string::npos){
    std::cout << "Given inFileName \'" << inFileName << "\' is not valid. return 1" << std::endl;
    return 1;
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* l1CaloTree_p = (TTree*)inFile_p->Get("l1CaloTowerEmuTree/L1CaloTowerTree");
  L1Analysis::L1AnalysisL1CaloTowerDataFormat *towers_ = new L1Analysis::L1AnalysisL1CaloTowerDataFormat();

  l1CaloTree_p->SetBranchStatus("*", 0);
  l1CaloTree_p->SetBranchStatus("L1CaloTower", 1);
  l1CaloTree_p->SetBranchStatus("iet", 1);
  l1CaloTree_p->SetBranchStatus("ieta", 1);
  l1CaloTree_p->SetBranchStatus("iphi", 1);
  l1CaloTree_p->SetBranchStatus("iem", 1);
  l1CaloTree_p->SetBranchStatus("ihad", 1);

  l1CaloTree_p->SetBranchAddress("L1CaloTower", &towers_);

  l1CaloTree_p->GetEntry(getEntry);

  int minEtaVal = 10000;
  int maxEtaVal = -10000;
  int minPhiVal = 10000;
  int maxPhiVal = -10000;

  std::map<int, std::map<int, int> > phiToEtaToEtMap;

  for(unsigned int i = 0; i < towers_->iet.size(); ++i){
    if(towers_->ieta.at(i) < minEtaVal) minEtaVal = towers_->ieta.at(i);
    if(towers_->ieta.at(i) > maxEtaVal) maxEtaVal = towers_->ieta.at(i);

    if(towers_->iphi.at(i) < minPhiVal) minPhiVal = towers_->iphi.at(i);
    if(towers_->iphi.at(i) > maxPhiVal) maxPhiVal = towers_->iphi.at(i);

    std::map<int, int> etaToEtMap;
    if(phiToEtaToEtMap.count(towers_->iphi.at(i)) == 0) etaToEtMap[towers_->ieta.at(i)] = towers_->iet.at(i);
    else{
      etaToEtMap = phiToEtaToEtMap[towers_->iphi.at(i)];
      etaToEtMap[towers_->ieta.at(i)] = towers_->iet.at(i);
    }

    phiToEtaToEtMap[towers_->iphi.at(i)] = etaToEtMap;
  }

  int minPhi = 1;
  int maxPhi = 72;
  int minEta = -41;
  int maxEta = 41;

  std::cout << " Min/Max eta val: " << minEtaVal << ", " << maxEtaVal << std::endl;
  std::cout << " Min/Max phi val: " << minPhiVal << ", " << maxPhiVal << std::endl;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  std::ofstream outFile(("output/" + dateStr + "/unwrapCaloEvt" + std::to_string(getEntry) + ".txt").c_str());

  const Int_t minStringSize = 4;
  outFile << "Phi,";
  
  for(int eI = minEta; eI <= maxEta; ++eI){
    if(eI < etaLow) continue;
    if(eI > etaHi) continue;
    if(eI == 0) continue;

    std::string etaStr = std::to_string(eI) + ",";
    while(etaStr.size() < minStringSize){etaStr = " " + etaStr;}
    outFile << etaStr;
  }
  outFile << std::endl;

  for(int pI = minPhi; pI <= maxPhi; ++pI){
    std::string phiStr = std::to_string(pI) + ","; 
    while(phiStr.size() < minStringSize){phiStr = " " + phiStr;}
    outFile << phiStr;

    std::map<int, int> etaToEtMap;
    if(phiToEtaToEtMap.count(pI) != 0) etaToEtMap = phiToEtaToEtMap[pI];

    for(int eI = minEta; eI <= maxEta; ++eI){
      if(eI < etaLow) continue;
      if(eI > etaHi) continue;
      if(eI == 0) continue;

      int etVal = 0;
      if(etaToEtMap.count(eI) != 0) etVal = etaToEtMap[eI];
      std::string etStr = std::to_string(etVal) + ",";
      while(etStr.size() < minStringSize){etStr = " " + etStr;}
      outFile << etStr;
    }
    
    outFile << std::endl;
  }

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 3 || argc > 5){
    std::cout << "Usage: ./bin/l1OfflineEvtDisp.exe <inFileName> <entry> <etaLow-opt> <etaHi-opt>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 3) retVal += l1OfflineEvtDisp(argv[1], std::stoi(argv[2]));
  else if(argc == 4) retVal += l1OfflineEvtDisp(argv[1], std::stoi(argv[2]), std::stoi(argv[3]));
  else if(argc == 5) retVal += l1OfflineEvtDisp(argv[1], std::stoi(argv[2]), std::stoi(argv[3]), std::stoi(argv[4]));
  return retVal;
}

#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include "include/doGlobalDebug.h"
#include "include/runLumiEvtKey.h"

std::string getL1AlgoFromFileName(const std::string inFileName)
{
  const Int_t nValidAlgo = 4;
  const std::string validAlgo[nValidAlgo] = {"None", "Donut", "ChunkyDonut", "PhiRingPP"};

  Int_t algoPos = -1;
  for(Int_t i = 0; i < nValidAlgo; ++i){
    if(inFileName.find(validAlgo[i]) != std::string::npos){
      algoPos = i;
      break;
    }
  }

  std::string outStr = "NOVALIDALGO";
  if(algoPos == 1 && inFileName.find("ChunkyDonut") != std::string::npos) outStr = "ChunkyDonut";
  else if(algoPos != -1) outStr = validAlgo[algoPos];

  return outStr;
}


int l1Comp(const std::string inL1Algo1Name, const std::string inL1Algo2Name, const std::string inForestName, std::string outFileName = "")
{
  if(outFileName.size() == 0){
    const Int_t nInStr = 3;
    std::string inStr[nInStr] = {inL1Algo1Name, inL1Algo2Name, inForestName};
    for(Int_t i = 0; i < nInStr; ++i){
      while(inStr[i].find("/") != std::string::npos){inStr[i].replace(0, inStr[i].find("/")+1, "");}
      while(inStr[i].find(".root") != std::string::npos){inStr[i].replace(inStr[i].find(".root"), std::string(".root").size(), "");}
    }

    outFileName = "l1Comp_" + getL1AlgoFromFileName(inStr[0]) + "_" + getL1AlgoFromFileName(inStr[1]) + "_" + inStr[2];
  }
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");
  TDatime* date = new TDatime();
  outFileName = outFileName + "_" + std::to_string(date->GetDate()) + ".root";

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t nL1Algo = 2;
  const std::string l1AlgoStr[nL1Algo] = {getL1AlgoFromFileName(inL1Algo1Name), getL1AlgoFromFileName(inL1Algo2Name)};
  runLumiEvtKey l1AlgoMap[nL1Algo] = {runLumiEvtKey(), runLumiEvtKey()};
  

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* l1JetEt_h[nL1Algo];
  for(Int_t i = 0; i < nL1Algo; ++i){
    l1JetEt_h[i] = new TH1F(("l1JetEt_" + l1AlgoStr[i] + "_h").c_str(), (";L1 Jet E_{T} (" + l1AlgoStr[i] + ");Counts").c_str(), 50, 20, 120);
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TFile* inL1Algo1_p = new TFile(inL1Algo1Name.c_str(), "READ");
  TTree* inL1Algo1EvtTree_p = (TTree*)inL1Algo1_p->Get("l1UpgradeEmuTree/L1UpgradeTree");
  TTree* inL1Algo1UpgradeTree_p = (TTree*)inL1Algo1_p->Get("l1UpgradeEmuTree/L1UpgradeTree");

  UInt_t run1, lumi1;
  ULong64_t event1;

  std::vector<float>* jetEt1_p=0;
  std::vector<float>* jetEta1_p=0;
  std::vector<float>* jetPhi1_p=0;

  inL1Algo1EvtTree_p->SetBranchStatus("*", 0);
  inL1Algo1EvtTree_p->SetBranchStatus("run", 1);
  inL1Algo1EvtTree_p->SetBranchStatus("lumi", 1);
  inL1Algo1EvtTree_p->SetBranchStatus("event", 1);

  inL1Algo1EvtTree_p->SetBranchAddress("run", &run1);
  inL1Algo1EvtTree_p->SetBranchAddress("lumi", &lumi1);
  inL1Algo1EvtTree_p->SetBranchAddress("event", &event1);

  inL1Algo1UpgradeTree_p->SetBranchStatus("*", 0);
  inL1Algo1UpgradeTree_p->SetBranchStatus("jetEt", 1);
  inL1Algo1UpgradeTree_p->SetBranchStatus("jetEta", 1);
  inL1Algo1UpgradeTree_p->SetBranchStatus("jetPhi", 1);

  inL1Algo1UpgradeTree_p->SetBranchAddress("jetEt", &jetEt1_p);
  inL1Algo1UpgradeTree_p->SetBranchAddress("jetEta", &jetEta1_p);
  inL1Algo1UpgradeTree_p->SetBranchAddress("jetPhi", &jetPhi1_p);

  const Int_t nEntriesL1Algo1 = inL1Algo1UpgradeTree_p->GetEntries();

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t entry = 0; entry < nEntriesL1Algo1; ++entry){
    inL1Algo1EvtTree_p->GetEntry(entry);
    inL1Algo1UpgradeTree_p->GetEntry(entry);

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    if(doGlobalDebug) std::cout << jetEt1_p->size() << std::endl;
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    l1AlgoMap[0]->addKey(run1, lumi1, event1, entry);

    for(unsigned int jI = 0; jI < jetEt1_p->size(); ++jI){
      l1JetEt_h[0]->Fill(jetEt1_p->at(jI));
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TFile* inL1Algo2_p = new TFile(inL1Algo2Name.c_str(), "READ");
  TTree* inL1Algo2EvtTree_p = (TTree*)inL1Algo2_p->Get("l1UpgradeEmuTree/L1UpgradeTree");
  TTree* inL1Algo2UpgradeTree_p = (TTree*)inL1Algo2_p->Get("l1UpgradeEmuTree/L1UpgradeTree");

  UInt_t run2, lumi2;
  ULong64_t event2;

  std::vector<float>* jetEt2_p=0;
  std::vector<float>* jetEta2_p=0;
  std::vector<float>* jetPhi2_p=0;

  inL1Algo2EvtTree_p->SetBranchStatus("*", 0);
  inL1Algo2EvtTree_p->SetBranchStatus("run", 1);
  inL1Algo2EvtTree_p->SetBranchStatus("lumi", 1);
  inL1Algo2EvtTree_p->SetBranchStatus("event", 1);

  inL1Algo2EvtTree_p->SetBranchAddress("run", &run2);
  inL1Algo2EvtTree_p->SetBranchAddress("lumi", &lumi2);
  inL1Algo2EvtTree_p->SetBranchAddress("event", &event2);

  inL1Algo2UpgradeTree_p->SetBranchStatus("*", 0);
  inL1Algo2UpgradeTree_p->SetBranchStatus("jetEt", 1);
  inL1Algo2UpgradeTree_p->SetBranchStatus("jetEta", 1);
  inL1Algo2UpgradeTree_p->SetBranchStatus("jetPhi", 1);

  inL1Algo2UpgradeTree_p->SetBranchAddress("jetEt", &jetEt2_p);
  inL1Algo2UpgradeTree_p->SetBranchAddress("jetEta", &jetEta2_p);
  inL1Algo2UpgradeTree_p->SetBranchAddress("jetPhi", &jetPhi2_p);


  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t nEntriesL1Algo2 = inL1Algo2UpgradeTree_p->GetEntries();

  for(Int_t entry = 0; entry < nEntriesL1Algo2; ++entry){
    inL1Algo2EvtTree_p->GetEntry(entry);
    inL1Algo2UpgradeTree_p->GetEntry(entry);

    l1AlgoMap[1]->addKey(run2, lumi2, event2, entry);

    for(unsigned int jI = 0; jI < jetEt2_p->size(); ++jI){
      l1JetEt_h[1]->Fill(jetEt2_p->at(jI));
    }
  }

  TFile* forestFile_p = new TFile(inForestName.c_str(), "READ");
  

  forestFile_p->Close();
  delete forestFile_p;

  inL1Algo1_p->Close();
  delete inL1Algo1_p;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  inL1Algo2_p->Close();
  delete inL1Algo2_p;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  outFile_p->cd();
  for(Int_t i = 0; i < nL1Algo; ++i){
    l1JetEt_h[i]->Write("", TObject::kOverwrite);
    delete l1JetEt_h[i];
  }
  outFile_p->Close();
  delete outFile_p;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 4 && argc != 5){
    std::cout << "Usage: ./l1Comp.exe <inL1Algo1Name> <inL1Algo2Name> <inForestName> <outFileName-optional>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 4) retVal += l1Comp(argv[1], argv[2], argv[3]);
  else if(argc == 5) retVal += l1Comp(argv[1], argv[2], argv[3], argv[4]);
  return retVal;
}

#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"
#include "TNamed.h"

#include "include/doGlobalDebug.h"
#include "include/runLumiEvtKey.h"
#include "include/returnRootFileContentsList.h"
#include "include/getLinBins.h"
#include "include/plotUtilities.h"
#include "include/histDefUtility.h"
#include "include/etaPhiFunc.h"

std::vector<std::string> removeDuplicates(std::vector<std::string> inStr)
{
  unsigned int pos = 0;
  while(inStr.size() > pos){
    bool isGood = true;

    for(unsigned int i = pos+1; i < inStr.size(); ++i){
      if(inStr.at(pos).size() == inStr.at(i).size() && inStr.at(pos).find(inStr.at(i)) != std::string::npos){
	isGood = false;
	inStr.erase(inStr.begin()+pos);
	break;
      }
    }

    if(isGood) ++pos;
  }

  return inStr;
}

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
  //handling the output file in event none is specified
  if(outFileName.size() == 0){
    const Int_t nInStr = 3;
    std::string inStr[nInStr] = {inL1Algo1Name, inL1Algo2Name, inForestName};
    for(Int_t i = 0; i < nInStr; ++i){
      while(inStr[i].find("/") != std::string::npos){inStr[i].replace(0, inStr[i].find("/")+1, "");}
      while(inStr[i].find(".root") != std::string::npos){inStr[i].replace(inStr[i].find(".root"), std::string(".root").size(), "");}
    }

    outFileName = "l1Comp_" + getL1AlgoFromFileName(inStr[0]) + "_" + getL1AlgoFromFileName(inStr[1]) + "_" + inStr[2];
  }
  //append date to ouptut for simple versioning
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");
  TDatime* date = new TDatime();
  outFileName = outFileName + "_" + std::to_string(date->GetDate()) + ".root";
  delete date;

  //Quickly grab number of jet algorithms + tree names in forest
  TFile* forestFile_p = new TFile(inForestName.c_str(), "READ");
  std::vector<std::string> jetTreeStr = removeDuplicates(returnRootFileContentsList(forestFile_p, "TTree", "JetAnalyzer"));
  forestFile_p->Close();
  delete forestFile_p;

  const Int_t nJetAlgos = (Int_t)jetTreeStr.size();
  std::string jetAlgos[nJetAlgos];
  for(Int_t jI = 0; jI < nJetAlgos; ++jI){
    jetAlgos[jI] = jetTreeStr.at(jI).substr(0, jetTreeStr.at(jI).find("Jet"));
  }

  //number of algos to compare will always be two but just maintain a baseline
  const Int_t nL1Algo = 2;
  
  const Int_t nL1JetThresholds = 5;
  const Float_t l1JetThresholds[nL1JetThresholds] = {8., 16., 24., 32., 40.};

  const std::string l1AlgoStr[nL1Algo] = {getL1AlgoFromFileName(inL1Algo1Name), getL1AlgoFromFileName(inL1Algo2Name)};
  runLumiEvtKey l1AlgoMap[nL1Algo] = {runLumiEvtKey(), runLumiEvtKey()};
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* l1JetEt_h[nL1Algo];
  TH1F* l1JetEta_h[nL1Algo];
  TH1F* l1JetPhi_h[nL1Algo];
  TH1F* jetPt_h[nJetAlgos];
  TH1F* jetPt_Trig_h[nJetAlgos][nL1Algo][nL1JetThresholds];

  for(Int_t i = 0; i < nL1Algo; ++i){
    l1JetEt_h[i] = new TH1F(("l1JetEt_" + l1AlgoStr[i] + "_h").c_str(), (";L1 Jet E_{T} (" + l1AlgoStr[i] + ");Counts").c_str(), 50, 20, 120);
    l1JetEta_h[i] = new TH1F(("l1JetEta_" + l1AlgoStr[i] + "_h").c_str(), (";L1 Jet E_{T} (" + l1AlgoStr[i] + ");Counts").c_str(), 50, 20, 120);
    l1JetPhi_h[i] = new TH1F(("l1JetPhi_" + l1AlgoStr[i] + "_h").c_str(), (";L1 Jet E_{T} (" + l1AlgoStr[i] + ");Counts").c_str(), 50, 20, 120);
    centerTitles({l1JetEt_h[i], l1JetEta_h[i], l1JetPhi_h[i]});
  }

  const Int_t nJtPtBins = 14;
  const Float_t jtPtLow = 30;
  const Float_t jtPtHi = 100;
  Double_t jtPtBins[nJtPtBins+1];
  getLinBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);

  for(Int_t jI = 0; jI < nJetAlgos; ++jI){
    jetPt_h[jI] = new TH1F(("jetPt_" + jetAlgos[jI] + "_h").c_str(), (";Jet p_{T} (" + jetAlgos[jI] + ");Counts").c_str(), nJtPtBins, jtPtBins);
    centerTitles({jetPt_h[jI]});

    for(Int_t lI = 0; lI < nL1Algo; ++lI){
      for(Int_t tI = 0; tI < nL1JetThresholds; ++tI){
	jetPt_Trig_h[jI][lI][tI] = new TH1F(("jetPt_" + jetAlgos[jI] + "_" + l1AlgoStr[lI] + "_Pt" + prettyString(l1JetThresholds[tI], 1, true) + "_h").c_str(), (";Jet p_{T} w/ Trigger " + l1AlgoStr[lI] + ", L1 p_{T} > " + prettyString(l1JetThresholds[tI], 1, false) + " (" + jetAlgos[jI] + ");Counts").c_str(), nJtPtBins, jtPtBins);
	centerTitles({jetPt_Trig_h[jI][lI][tI]});
      }
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TFile* inL1Algo1_p = new TFile(inL1Algo1Name.c_str(), "READ");
  TTree* inL1Algo1EvtTree_p = (TTree*)inL1Algo1_p->Get("l1EventTree/L1EventTree");
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

    l1AlgoMap[0].addKey(run1, lumi1, event1, entry);

    for(unsigned int jI = 0; jI < jetEt1_p->size(); ++jI){
      l1JetEt_h[0]->Fill(jetEt1_p->at(jI));
      l1JetEta_h[0]->Fill(jetEta1_p->at(jI));
      l1JetPhi_h[0]->Fill(jetPhi1_p->at(jI));
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TFile* inL1Algo2_p = new TFile(inL1Algo2Name.c_str(), "READ");
  TTree* inL1Algo2EvtTree_p = (TTree*)inL1Algo2_p->Get("l1EventTree/L1EventTree");
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

    l1AlgoMap[1].addKey(run2, lumi2, event2, entry);

    for(unsigned int jI = 0; jI < jetEt2_p->size(); ++jI){
      l1JetEt_h[1]->Fill(jetEt2_p->at(jI));
      l1JetEta_h[1]->Fill(jetEta2_p->at(jI));
      l1JetPhi_h[1]->Fill(jetPhi2_p->at(jI));
    }
  }

  UInt_t runF, lumiF;
  ULong64_t eventF;

  const Int_t nMaxJets = 500;
  Int_t nref_[nJetAlgos];
  Float_t jtpt_[nJetAlgos][nMaxJets];
  Float_t jteta_[nJetAlgos][nMaxJets];
  Float_t jtphi_[nJetAlgos][nMaxJets];

  forestFile_p = new TFile(inForestName.c_str(), "READ");  
  TTree* hiTree_p = (TTree*)forestFile_p->Get("hiEvtAnalyzer/HiTree");
  TTree* jetTree_p[nJetAlgos];

  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("run", 1);
  hiTree_p->SetBranchStatus("lumi", 1);
  hiTree_p->SetBranchStatus("evt", 1);

  hiTree_p->SetBranchAddress("run", &runF);
  hiTree_p->SetBranchAddress("lumi", &lumiF);
  hiTree_p->SetBranchAddress("evt", &eventF);

  for(unsigned int i = 0; i < jetTreeStr.size(); ++i){
    jetTree_p[i] = (TTree*)forestFile_p->Get(jetTreeStr.at(i).c_str());

    jetTree_p[i]->SetBranchStatus("*", 0);
    jetTree_p[i]->SetBranchStatus("nref", 1);
    jetTree_p[i]->SetBranchStatus("jtpt", 1);
    jetTree_p[i]->SetBranchStatus("jteta", 1);
    jetTree_p[i]->SetBranchStatus("jtphi", 1);

    jetTree_p[i]->SetBranchAddress("nref", &nref_[i]);
    jetTree_p[i]->SetBranchAddress("jtpt", jtpt_[i]);
    jetTree_p[i]->SetBranchAddress("jteta", jteta_[i]);
    jetTree_p[i]->SetBranchAddress("jtphi", jtphi_[i]);
  }

  const Int_t nEntriesForest = hiTree_p->GetEntries();

  Int_t totalFound = 0;

  std::cout << "Processing forest..." << std::endl;
  for(Int_t entry = 0; entry < nEntriesForest; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntriesForest << std::endl;

    hiTree_p->GetEntry(entry);
    
    Int_t entryL1Algo1 = l1AlgoMap[0].getEntryFromKey(runF, lumiF, eventF);
    if(entryL1Algo1 < 0) continue;
    Int_t entryL1Algo2 = l1AlgoMap[1].getEntryFromKey(runF, lumiF, eventF);
    if(entryL1Algo2 < 0) continue;

    totalFound++;
    //    std::cout << "  Found entry!" << std::endl;

    for(unsigned int i = 0; i < jetTreeStr.size(); ++i){jetTree_p[i]->GetEntry(entry);}
    inL1Algo1UpgradeTree_p->GetEntry(entryL1Algo1);
    inL1Algo2UpgradeTree_p->GetEntry(entryL1Algo2);
    
    for(unsigned int i = 0; i < jetTreeStr.size(); ++i){
      Float_t leadingJtPt_ = -999;
      Float_t leadingJtPhi_ = -999;
      Float_t leadingJtEta_ = -999;

      for(Int_t jI = 0; jI < nref_[i]; ++jI){
	if(TMath::Abs(jteta_[i][jI]) > 2.) continue;

	if(jtpt_[i][jI] > leadingJtPt_){
	  leadingJtPt_ = jtpt_[i][jI];
	  leadingJtPhi_ = jtphi_[i][jI];
	  leadingJtEta_ = jteta_[i][jI];
	}
      }

      if(leadingJtPt_ > jtPtLow && leadingJtPt_ < jtPtHi){
	jetPt_h[i]->Fill(leadingJtPt_);

	Float_t matchedL1Pt_[nL1Algo] = {-999, -999};

	for(unsigned int lI = 0; lI < jetEt1_p->size(); ++lI){
	  if(getDR(jetEta1_p->at(lI), jetPhi1_p->at(lI), leadingJtEta_, leadingJtPhi_) > 2.) continue;
	  if(matchedL1Pt_[0] < jetEt1_p->at(lI)) matchedL1Pt_[0] = jetEt1_p->at(lI);
	}

	for(unsigned int lI = 0; lI < jetEt2_p->size(); ++lI){
	  if(getDR(jetEta2_p->at(lI), jetPhi2_p->at(lI), leadingJtEta_, leadingJtPhi_) > 2.) continue;
	  if(matchedL1Pt_[1] < jetEt2_p->at(lI)) matchedL1Pt_[1] = jetEt2_p->at(lI);
	}

	for(Int_t lI = 0; lI < nL1Algo; ++lI){
	  for(Int_t tI = 0; tI < nL1JetThresholds; ++tI){
	    if(matchedL1Pt_[lI] > l1JetThresholds[tI]) jetPt_Trig_h[i][lI][tI]->Fill(leadingJtPt_);
	    else break;
	  }
	}
      }
    }
  }

  std::cout << "Total found: " << totalFound << std::endl;

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

    l1JetEta_h[i]->Write("", TObject::kOverwrite);
    delete l1JetEta_h[i];

    l1JetPhi_h[i]->Write("", TObject::kOverwrite);
    delete l1JetPhi_h[i];
  }

  for(Int_t jI = 0; jI < nJetAlgos; ++jI){
    jetPt_h[jI]->Write("", TObject::kOverwrite);

    for(Int_t lI = 0; lI < nL1Algo; ++lI){
      for(Int_t tI = 0; tI < nL1JetThresholds; ++tI){
	jetPt_Trig_h[jI][lI][tI]->Write("", TObject::kOverwrite);
	delete jetPt_Trig_h[jI][lI][tI];
      }
    }

    delete jetPt_h[jI];
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

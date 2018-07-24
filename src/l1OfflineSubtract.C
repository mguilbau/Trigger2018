#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <fstream>
//#include <pair>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TMath.h"
#include "TDatime.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"

#include "include/checkMakeDir.h"
#include "include/doGlobalDebug.h"
#include "include/etaPhiFunc.h"
#include "include/runLumiEvtKey.h"
#include "include/getLinBins.h"
#include "include/L1Tools.h"
#include "include/mntToXRootdFileString.h"
#include "include/L1AnalysisEventDataFormat.h"
#include "include/L1AnalysisL1CaloTowerDataFormat.h"
#include "include/L1AnalysisL1UpgradeDataFormat.h"
#include "include/returnRootFileContentsList.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/vanGoghPalette.h"

const int mask_[9][9] = {
  { 1,2,2,2,2,2,2,2,2 },
  { 1,1,2,2,2,2,2,2,2 },
  { 1,1,1,2,2,2,2,2,2 },
  { 1,1,1,1,2,2,2,2,2 },
  { 1,1,1,1,0,2,2,2,2 },
  { 1,1,1,1,1,2,2,2,2 },
  { 1,1,1,1,1,1,2,2,2 },
  { 1,1,1,1,1,1,1,2,2 },
  { 1,1,1,1,1,1,1,1,2 },
};


void getXYBins(const int nBins, int* x, int* y)
{
  if(nBins == 1){
    (*x) = 1;
    (*y) = 1;
  }
  else if(nBins == 2){
    (*x) = 2;
    (*y) = 1;
  }
  else if(nBins == 3){
    (*x) = 3;
    (*y) = 1;
  }
  else if(nBins == 4){
    (*x) = 4;
    (*y) = 1;
  }
  else if(nBins == 5){
    (*x) = 3;
    (*y) = 2;
  }
  else if(nBins == 6){
    (*x) = 3;
    (*y) = 2;
  }
  else if(nBins == 7){
    (*x) = 4;
    (*y) = 2;
  }
  else if(nBins == 8){
    (*x) = 4;
    (*y) = 2;
  }
  
  return;
}


int l1OfflineSubtract(std::vector<std::string> inFileName, std::vector<std::string> forestFileName, const bool doWeights, const std::string extraTag, std::vector<std::string> subTypeStr)
{
  std::cout << "DO WEIGHTS: " << doWeights << std::endl;

  if(inFileName.size() != forestFileName.size()){
    std::cout << "Size of inFileName \'" << inFileName.size() << "\', but forestFileName is size \'" << forestFileName.size() << "\'... return 1." << std::endl;
    return 1;
  }

  bool isInFile = true;
  bool isForestFile = true;
  bool doForest = true;

  for(unsigned int fI = 0; fI < inFileName.size(); ++fI){
    isInFile = isInFile && checkFile(inFileName.at(fI));
    isForestFile = isForestFile && (checkFile(forestFileName.at(fI)) || forestFileName.at(fI).find("/mnt") != std::string::npos);
    doForest = doForest && checkFile(forestFileName.at(fI)) && forestFileName.at(fI).size() != 0;
  }

  if(!isInFile){
    std::cout << "Given input file \'" << inFileName.at(0) << "\', etc. is invalid. return 1." << std::endl;
    return 1;
  }

  std::vector<std::string> jetTrees;

  const Int_t nForestFile = forestFileName.size();

  UInt_t runF, lumiF;
  ULong64_t eventF;
  Float_t vz_;
  Float_t rho_;

  bool hasHI = false;
  bool hasRho = false;

  //std::cout << __LINE__ << std::endl;
  runLumiEvtKey* forestMap[nForestFile];
  TTree* hiTree_p = NULL;
  for(Int_t fI = 0; fI < nForestFile; ++fI){
    forestMap[fI] = NULL;
  }

  if(doForest){
    for(unsigned int fI = 0; fI < forestFileName.size(); ++fI){
      //std::cout << __LINE__ << std::endl;
      
      TFile* forestFile_p = TFile::Open(mntToXRootdFileString(forestFileName.at(fI)).c_str(), "READ");
      jetTrees = returnRootFileContentsList(forestFile_p, "TTree", "ak");
      unsigned int pos = 0;
      while(pos < jetTrees.size()){
	if(jetTrees.at(pos).substr(0,2).find("ak") != std::string::npos) pos++;
	else jetTrees.erase(jetTrees.begin()+pos);
      }
      
      //std::cout << __LINE__ << std::endl;
      
      forestMap[fI] = new runLumiEvtKey();
      
      std::vector<std::string> allTrees = returnRootFileContentsList(forestFile_p, "TTree", "");
      
      for(unsigned int aI = 0; aI < allTrees.size(); ++aI){
	if(allTrees.at(aI).find("hiEvtAnalyzer") != std::string::npos) hasHI = true;
	else if(allTrees.at(aI).find("rcRhoR4N11HiNtuplizer") != std::string::npos) hasRho = true;
      }
      
      //std::cout << __LINE__ << std::endl;
      
      Int_t HBHENoiseFilterResultRun2Loose_;
      Int_t pcollisionEventSelection_;
      
      TTree* skimTree_p = (TTree*)forestFile_p->Get("skimanalysis/HltTree");
      skimTree_p->SetBranchStatus("*", 0);
      skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
      skimTree_p->SetBranchStatus("pcollisionEventSelection", 1);
      
      skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);
      skimTree_p->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection_);
      
      if(hasHI) hiTree_p = (TTree*)forestFile_p->Get("hiEvtAnalyzer/HiTree");
      else if(hasRho) hiTree_p = (TTree*)forestFile_p->Get("rcRhoR4N11HiNtuplizer/EventTree");
      
      hiTree_p->SetBranchStatus("*", 0);
      hiTree_p->SetBranchStatus("run", 1);
      hiTree_p->SetBranchStatus("lumi", 1);
      hiTree_p->SetBranchStatus("evt", 1);
      hiTree_p->SetBranchStatus("vz", 1);
      hiTree_p->SetBranchStatus("rho", 1);
      
      hiTree_p->SetBranchAddress("run", &runF);
      hiTree_p->SetBranchAddress("lumi", &lumiF);
      hiTree_p->SetBranchAddress("evt", &eventF);
      hiTree_p->SetBranchAddress("vz", &vz_);
      hiTree_p->SetBranchAddress("rho", &rho_);
      
      const Int_t nEntries = hiTree_p->GetEntries();
      
      //std::cout << __LINE__ << std::endl;
      for(Int_t entry = 0; entry < nEntries; ++entry){
	hiTree_p->GetEntry(entry);
	skimTree_p->GetEntry(entry);
	
	if(!pcollisionEventSelection_) continue;
	if(TMath::Abs(vz_) > 15.) continue;
	
	forestMap[fI]->addKey(runF, lumiF, eventF, entry);
      }
      
      //std::cout << __LINE__ << std::endl;
      
      forestFile_p->Close();
      delete forestFile_p;
    }
  }

  std::cout << "Do Forest: " << doForest << std::endl;
  
  //std::cout << __LINE__ << std::endl;

  const Int_t nJetTrees = jetTrees.size();
  std::vector<std::string> jetNames;
  for(unsigned int jI = 0; jI < jetTrees.size(); ++jI){
    std::string tempStr = jetTrees.at(jI);
    while(tempStr.find("/") != std::string::npos){tempStr.replace(tempStr.find("/"), 1, "_");}
    jetNames.push_back(tempStr);
  }

  const Int_t nSubType = 14;
  const std::string subType[nSubType] = {"None", "PhiRingHITower", "PhiRingPPTower", "PhiRingHIRegion", "ChunkyDonut", "ChunkyDonutLUT", "ChunkyDonutZero", "PhiRingPPTower72", "PhiRingPPTower72NoDiv", "PhiRingPPTower72NoDivCorr", "PhiRingHITower144", "PhiRingHIRegion288", "ChunkyDonutHIMod", "PhiRingPPTower72NoDivNoNeg"};
  const Int_t styles[nSubType] = {20, 47, 34, 21, 43, 29, 24, 25, 46, 45, 20, 27, 28, 49};
  vanGoghPalette vg;
  bool doType[nSubType] = {false, false, false, false, false, false, false, false, false, false, false, false, false, false};
  bool isFound = false;

  int firstTypePos = -1;

  bool isFileType[nSubType] = {false, false, false, false, false, false, false, false, false, false, false, false, false, false};
  
  std::string globalAlgoString = "";

  const Int_t nL1Thresh = 13;
  const Float_t l1ThreshPt[nL1Thresh] = {8., 16., 24., 32., 40., 48., 56., 64., 72., 80., 96., 112., 128.};
  Double_t l1ThreshPtBins[nL1Thresh+1];
  l1ThreshPtBins[0] = l1ThreshPt[0] - (l1ThreshPt[1] - l1ThreshPt[0])/2.;
  l1ThreshPtBins[nL1Thresh] = l1ThreshPt[nL1Thresh-1] + (l1ThreshPt[nL1Thresh-1] - l1ThreshPt[nL1Thresh-2])/2.;
 
  TH1F* dummys_p[nSubType];
  TH1F* dummys2_p[nJetTrees];
  TH1F* dummys4_p[nSubType];
  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);

  TLegend* dummyLeg_p = new TLegend(0.4, 0.2, 0.9, 0.5);
  dummyLeg_p->SetBorderSize(0);
  dummyLeg_p->SetFillStyle(0);
  dummyLeg_p->SetFillColor(0);
  dummyLeg_p->SetTextFont(43);
  dummyLeg_p->SetTextSize(14);

  TLegend* dummyLeg2_p = new TLegend(0.15, 0.65, 0.4, 0.9);
  dummyLeg2_p->SetBorderSize(0);
  dummyLeg2_p->SetFillStyle(0);
  dummyLeg2_p->SetFillColor(0);
  dummyLeg2_p->SetTextFont(43);
  dummyLeg2_p->SetTextSize(14);


  TLegend* dummyLeg3_p = new TLegend(0.2, 0.2, 0.6, 0.5);
  dummyLeg3_p->SetBorderSize(0);
  dummyLeg3_p->SetFillStyle(0);
  dummyLeg3_p->SetFillColor(0);
  dummyLeg3_p->SetTextFont(43);
  dummyLeg3_p->SetTextSize(14);

  TLegend* dummyLeg4_p = new TLegend(0.4, 0.2, 0.9, 0.5);
  dummyLeg4_p->SetBorderSize(0);
  dummyLeg4_p->SetFillStyle(0);
  dummyLeg4_p->SetFillColor(0);
  dummyLeg4_p->SetTextFont(43);
  dummyLeg4_p->SetTextSize(14);

  for(Int_t i = 0; i < nSubType; ++i){
    std::string probeStr = "_" + subType[i] + "_";
    if(inFileName.at(0).find(probeStr) != std::string::npos){
      isFileType[i] = true;
      std::cout << "infile \'" << inFileName.at(0) << "\' is given type \'" << subType[i] << "\'" << std::endl;
      break;
    }
  }  


  for(Int_t i = 0; i < nSubType; ++i){
    dummys_p[i] = new TH1F(("dummy_" + std::to_string(i)).c_str(),"", 10, 0, 1);
    dummys_p[i]->SetMarkerColor(vg.getColor(i));
    dummys_p[i]->SetLineColor(vg.getColor(i));
    dummys_p[i]->SetMarkerStyle(styles[i]);
    dummys_p[i]->SetMarkerSize(1.2);
    setSumW2(dummys_p[i]);
    centerTitles(dummys_p[i]);

    for(unsigned int j = 0; j < subTypeStr.size(); ++j){
      if(subTypeStr.at(j).size() == subType[i].size() && subTypeStr.at(j).find(subType[i]) != std::string::npos){
	isFound = true;
	doType[i] = true;
	if(firstTypePos < 0) firstTypePos = i;

	globalAlgoString = globalAlgoString + "_" + subType[i];
	dummyLeg_p->AddEntry(dummys_p[i], subType[i].c_str(), "P L");
	dummyLeg3_p->AddEntry(dummys_p[i], subType[i].c_str(), "P L");
      }
    }
  }

  for(Int_t i = 0; i < nL1Thresh; ++i){
    dummys4_p[i] = new TH1F(("dummy4_" + std::to_string(i)).c_str(),"", 10, 0, 1);
    dummys4_p[i]->SetMarkerColor(vg.getColor(i));
    dummys4_p[i]->SetLineColor(vg.getColor(i));
    dummys4_p[i]->SetMarkerStyle(styles[i]);
    dummys4_p[i]->SetMarkerSize(1.2);
    setSumW2(dummys4_p[i]);
    centerTitles(dummys4_p[i]);

    std::string l1Str = "p_{T,L1} > " + prettyString(l1ThreshPt[i], 1, false);

    dummyLeg4_p->AddEntry(dummys4_p[i], l1Str.c_str(), "P L");
  }

  if(globalAlgoString.size() != 0) globalAlgoString.replace(0,1,"");


  for(Int_t i = 0; i < nJetTrees; ++i){
    dummys2_p[i] = new TH1F(("dummy2_" + std::to_string(i)).c_str(),"", 10, 0, 1);
    dummys2_p[i]->SetMarkerColor(vg.getColor(i));
    dummys2_p[i]->SetLineColor(vg.getColor(i));
    dummys2_p[i]->SetMarkerStyle(styles[i]);
    dummys2_p[i]->SetMarkerSize(1.2);
    setSumW2(dummys2_p[i]);
    centerTitles(dummys2_p[i]);

    dummyLeg2_p->AddEntry(dummys2_p[i], jetNames.at(i).c_str(), "P L");
  }

  if(!isFound){
    std::cout << "Given subtraction types are not found. Please pick one of: " << std::endl;
    std::cout << " ";
    for(Int_t i = 0; i < nSubType; ++i){
      std::cout << subType[i] << ",";
    }
    std::cout << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);

  //std::cout << __LINE__ << std::endl;

  const Double_t zeroedEtaLow = towerEta(24);
  const Double_t zeroedEtaHi = towerEta(29);

  const Int_t nRhoBins = 3;
  const Int_t rhoBinsLow[nRhoBins] = {150, 50, 0};
  const Int_t rhoBinsHi[nRhoBins] = {1000, 150, 50};

 
  //std::cout << __LINE__ << std::endl;

  for(Int_t tI = 1; tI < nL1Thresh; ++tI){
    l1ThreshPtBins[tI] = l1ThreshPt[tI] - (l1ThreshPt[tI] - l1ThreshPt[tI-1])/2;
  }

  Int_t triggerOpp[nRhoBins+1][nSubType][nL1Thresh];
  Int_t triggerFires[nRhoBins+1][nSubType][nL1Thresh];
  

  //std::cout << __LINE__ << std::endl;

  for(Int_t rI = 0; rI < nRhoBins+1; ++rI){
    for(Int_t sI = 0; sI < nSubType; ++sI){
      for(Int_t lI = 0; lI < nL1Thresh; ++lI){
	triggerOpp[rI][sI][lI] = 0;
	triggerFires[rI][sI][lI] = 0;
      }
    }
  }

  //std::cout << __LINE__ << std::endl;

  const int minSeedThresh = 8;
  const int nIEta = 83;
  const int nIPhi = 72;

  std::string subtractType = "";
  for(Int_t sI = 0; sI < nSubType; ++sI){
    if(doType[sI]) subtractType = subtractType + "_" + subType[sI] + "_";
  }
  subtractType.replace(0,1,"");
  subtractType.replace(subtractType.size()-1, 1, "");


  TFile* outFile_p = new TFile(("output/l1OfflineSubtract_" + subtractType + "_" + extraTag + "_" + dateStr + ".root").c_str(), "RECREATE");

  TH1F* pthat_p = new TH1F("pthat_h", ";ptHat;Counts", 100, 15, 215);
  TH1F* pthatWeighted_p = new TH1F("pthatWeighted_h", ";ptHat (Weighted);Counts", 100, 15, 215);

  TH1F* etaMiss_p = new TH1F("etaMiss_h", ";HW #eta;Counts (Misses)", 103, -51.5, 51.5);
  TH1F* phiMiss_p = new TH1F("phiMiss_h", ";HW #phi;Counts (Misses)", 150, -0.5, 149.5);

  TH1F* etaMiss_CMSSWJet_p = new TH1F("etaMiss_CMSSWJet_h", ";HW #eta;Counts (Misses)", 103, -51.5, 51.5);
  TH1F* phiMiss_CMSSWJet_p = new TH1F("phiMiss_CMSSWJet_h", ";HW #phi;Counts (Misses)", 150, -0.5, 149.5);

  setSumW2({pthat_p, pthatWeighted_p, etaMiss_p, phiMiss_p, etaMiss_CMSSWJet_p, phiMiss_CMSSWJet_p});
  centerTitles({pthat_p, pthatWeighted_p, etaMiss_p, phiMiss_p, etaMiss_CMSSWJet_p, phiMiss_CMSSWJet_p});

  //std::cout << __LINE__ << std::endl;

  Int_t tempNJtPtBins = 18;
  Int_t tempJtPtLow = 10;
  Int_t tempJtPtHi = 100;

  Int_t tempNJtPtBinsFor = 10;
  Int_t tempJtPtLowFor = 10;
  Int_t tempJtPtHiFor = 60;

  if(extraTag.find("80") != std::string::npos){
   tempNJtPtBins = 26;
   tempJtPtLow = 10;
   tempJtPtHi = 140;

   tempNJtPtBinsFor = 18;
   tempJtPtLowFor = 10;
   tempJtPtHiFor = 100;
  }
  
  const Int_t nJtPtBins = tempNJtPtBins;
  const Float_t jtPtLow = tempJtPtLow;
  const Float_t jtPtHi = tempJtPtHi;

  const Int_t nJtPtBinsFor = tempNJtPtBinsFor;
  const Float_t jtPtLowFor = tempJtPtLowFor;
  const Float_t jtPtHiFor = tempJtPtHiFor;

  Double_t jtPtBins[nJtPtBins+1];
  Double_t jtPtBinsFor[nJtPtBinsFor+1];
  getLinBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);
  getLinBins(jtPtLowFor, jtPtHiFor, nJtPtBinsFor, jtPtBinsFor);

  const Int_t nJtAbsEtaBins = 2;
  const Float_t jtAbsEtaBinsLow[nJtAbsEtaBins] = {0.0, 3.0};
  const Float_t jtAbsEtaBinsHi[nJtAbsEtaBins] = {3.0, 5.0};

  TH1F* leadingJetPt_h[nJetTrees][nRhoBins][nJtAbsEtaBins+1];
  TH1F* leadingJetPt_L1Pt_h[nJetTrees][nRhoBins][nJtAbsEtaBins+1][nSubType][nL1Thresh];
  TH1F* leadingJetPt_L1Pt_DR_h[nJetTrees][nRhoBins][nJtAbsEtaBins+1][nSubType][nL1Thresh];

  TH1F* triggerEta_L1Pt_h[nRhoBins+1][nSubType][nL1Thresh];
  TH1F* jetEta_L1Pt_h[nRhoBins+1][nJetTrees][nL1Thresh];
  
  const Int_t nL1ResPtBins = 5;
  const Double_t l1ResPtBins[nL1ResPtBins+1] = {20, 30, 40, 60, 100, 160};
  TH1F* jetL1PtOverOffline_Mean_h[nJetTrees][nRhoBins+1][nSubType];
  TH1F* jetL1PtOverOffline_Sigma_h[nJetTrees][nRhoBins+1][nSubType];
  TH1F* jetL1PtOverOffline_h[nJetTrees][nRhoBins+1][nSubType][nL1ResPtBins];

  const Double_t maxDR = 0.6;

  if(doForest){
    for(Int_t jI = 0; jI < nJetTrees; ++jI){
      std::string jetStr = "";
      if(jetTrees.at(jI).find("Calo") != std::string::npos) jetStr = "Calo. ";
      else if(jetTrees.at(jI).find("PF") != std::string::npos) jetStr = "PF ";
      else if(jetTrees.at(jI).find("Gen") != std::string::npos) jetStr = "Gen. ";

      for(Int_t cI = 0; cI < nRhoBins; ++cI){
	const std::string centStr = "Rho" + std::to_string(rhoBinsLow[cI]) + "to" + std::to_string(rhoBinsHi[cI]);
	const std::string centStr2 = std::to_string(rhoBinsLow[cI]) + " < #rho <" + std::to_string(rhoBinsHi[cI]);
	
	for(Int_t aI = 0; aI < nJtAbsEtaBins+1; ++aI){
	  std::string jtAbsEtaBinsStr = "JtAbsEtaAll";
	  std::string jtAbsEtaBinsStr2 = "0.0 < |#eta| < 5.0";
	  if(aI != nJtAbsEtaBins){
	    jtAbsEtaBinsStr = "JtAbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	    jtAbsEtaBinsStr2 = prettyString(jtAbsEtaBinsLow[aI], 1, false) + " < |#eta| < " + prettyString(jtAbsEtaBinsHi[aI], 1, false);
	  }

	  leadingJetPt_h[jI][cI][aI] = NULL;
	  if(jtAbsEtaBinsStr.find("AbsEta3p0") != std::string::npos) leadingJetPt_h[jI][cI][aI] = new TH1F(("leadingJetPt_" + jetNames.at(jI) + "_" + centStr + "_" + jtAbsEtaBinsStr + "_h").c_str(), (";" + jetStr + "Jet p_{T} (" + jtAbsEtaBinsStr2 + ");Counts (" + centStr2 + ")").c_str(), nJtPtBinsFor, jtPtBinsFor);
	  else leadingJetPt_h[jI][cI][aI] = new TH1F(("leadingJetPt_" + jetNames.at(jI) + "_" + centStr + "_" + jtAbsEtaBinsStr + "_h").c_str(), (";Jet p_{T};Counts (" + centStr2 + ")").c_str(), nJtPtBins, jtPtBins);
      
	  setSumW2({leadingJetPt_h[jI][cI][aI]});
	  centerTitles({leadingJetPt_h[jI][cI][aI]});

	  for(Int_t sI = 0; sI < nSubType; ++sI){
	    for(Int_t lI = 0; lI < nL1Thresh; ++lI){
	      
	      leadingJetPt_L1Pt_h[jI][cI][aI][sI][lI] = NULL;
	      leadingJetPt_L1Pt_DR_h[jI][cI][aI][sI][lI] = NULL;

	      if(jtAbsEtaBinsStr.find("AbsEta3p0") != std::string::npos){
		leadingJetPt_L1Pt_h[jI][cI][aI][sI][lI] = new TH1F(("leadingJetPt_" + jetNames.at(jI) + "_" + centStr + "_" + jtAbsEtaBinsStr + "_" + subType[sI] + "_L1Pt" + prettyString(l1ThreshPt[lI], 1, true) + "_h").c_str(), (";" + jetStr + "Jet p_{T} (" + jtAbsEtaBinsStr2 + ");Counts (" + centStr2 + ")").c_str(), nJtPtBinsFor, jtPtBinsFor);
		leadingJetPt_L1Pt_DR_h[jI][cI][aI][sI][lI] = new TH1F(("leadingJetPt_" + jetNames.at(jI) + "_" + centStr + "_" + jtAbsEtaBinsStr + "_" + subType[sI] + "_L1Pt" + prettyString(l1ThreshPt[lI], 1, true) + "_DR_h").c_str(), (";" +  jetStr + "Jet p_{T} (" + jtAbsEtaBinsStr2 + ");Counts (" + centStr2 + ")").c_str(), nJtPtBinsFor, jtPtBinsFor);
	      }
	      else{
		leadingJetPt_L1Pt_h[jI][cI][aI][sI][lI] = new TH1F(("leadingJetPt_" + jetNames.at(jI) + "_" + centStr + "_" + jtAbsEtaBinsStr + "_" + subType[sI] + "_L1Pt" + prettyString(l1ThreshPt[lI], 1, true) + "_h").c_str(), (";Jet p_{T};Counts (" + centStr2 + ")").c_str(), nJtPtBins, jtPtBins);
		leadingJetPt_L1Pt_DR_h[jI][cI][aI][sI][lI] = new TH1F(("leadingJetPt_" + jetNames.at(jI) + "_" + centStr + "_" + jtAbsEtaBinsStr + "_" + subType[sI] + "_L1Pt" + prettyString(l1ThreshPt[lI], 1, true) + "_DR_h").c_str(), (";Jet p_{T};Counts (" + centStr2 + ")").c_str(), nJtPtBins, jtPtBins);
	      }

	      setSumW2({leadingJetPt_L1Pt_h[jI][cI][aI][sI][lI]});
	      setSumW2({leadingJetPt_L1Pt_DR_h[jI][cI][aI][sI][lI]});
	      centerTitles({leadingJetPt_L1Pt_h[jI][cI][aI][sI][lI]});
	      centerTitles({leadingJetPt_L1Pt_DR_h[jI][cI][aI][sI][lI]});
	    }
	  }
	}
      }
    }
  }

  for(Int_t cI = 0; cI < nRhoBins+1; ++cI){
    std::string centStr = "AllRho";
    if(cI != nRhoBins) centStr = "Rho" + std::to_string(rhoBinsLow[cI]) + "to" + std::to_string(rhoBinsHi[cI]);
	
    for(Int_t sI = 0; sI < nSubType; ++sI){
      for(Int_t lI = 0; lI < nL1Thresh; ++lI){
	
	triggerEta_L1Pt_h[cI][sI][lI] = new TH1F(("triggerEta_" + centStr + "_" + subType[sI] + "_L1Pt" + prettyString(l1ThreshPt[lI], 1, true) + "_h").c_str(), (centStr + ";L1 Jet #eta (L1 p_{T} > " + prettyString(l1ThreshPt[lI], 1, false) + ");Counts").c_str(), 85, -42.5, 42.5);
	centerTitles(triggerEta_L1Pt_h[cI][sI][lI]);
	setSumW2(triggerEta_L1Pt_h[cI][sI][lI]);
      }
    }

    for(Int_t jI = 0; jI < nJetTrees; ++jI){
      for(Int_t lI = 0; lI < nL1Thresh; ++lI){       
	jetEta_L1Pt_h[cI][jI][lI] = new TH1F(("triggerEta_" + centStr + "_" + jetNames.at(jI) + "_JetPt" + prettyString(l1ThreshPt[lI], 1, true) + "_h").c_str(), (centStr + ";Jet #eta (p_{T} > " + prettyString(l1ThreshPt[lI], 1, false) + ");Counts").c_str(), 101, -5.1, 5.1);
	centerTitles(jetEta_L1Pt_h[cI][jI][lI]);
	setSumW2(jetEta_L1Pt_h[cI][jI][lI]);
      }

      for(Int_t sI = 0; sI < nSubType; ++sI){
	jetL1PtOverOffline_Mean_h[jI][cI][sI] = new TH1F(("jetL1PtOverOffline_" + centStr + "_" + subType[sI] + "_" + jetNames.at(jI) + "_Mean_h").c_str(), ";L1/Offline;Counts", nL1ResPtBins, l1ResPtBins);
	jetL1PtOverOffline_Sigma_h[jI][cI][sI] = new TH1F(("jetL1PtOverOffline_" + centStr + "_" + subType[sI] + "_" + jetNames.at(jI) + "_Sigma_h").c_str(), ";L1/Offline;Counts", nL1ResPtBins, l1ResPtBins);
	
	centerTitles({jetL1PtOverOffline_Mean_h[jI][cI][sI], jetL1PtOverOffline_Sigma_h[jI][cI][sI] });
	setSumW2({jetL1PtOverOffline_Mean_h[jI][cI][sI], jetL1PtOverOffline_Sigma_h[jI][cI][sI] });

	for(Int_t rI = 0; rI < nL1ResPtBins; ++rI){
	  jetL1PtOverOffline_h[jI][cI][sI][rI] = new TH1F(("jetL1PtOverOffline_" + centStr + "_" + subType[sI] + "_" + jetNames.at(jI) + "_OfflinePt" + prettyString(l1ResPtBins[rI], 1, true) + "to" + prettyString(l1ResPtBins[rI+1], 1, true) + "_h").c_str(), ";L1/Offline;Counts", 51, 0., 5.);

	  centerTitles(jetL1PtOverOffline_h[jI][cI][sI][rI]);
	  setSumW2(jetL1PtOverOffline_h[jI][cI][sI][rI]);
	}
      }
    } 
  }



  int totalMatch = 0;
  int foundCount = 0;
  int missCount = 0;

  const Int_t nMaxJet = 500;
  TFile* forestFile_p=NULL;
  TTree* jetTrees_p[nJetTrees];
  
  Float_t pthat_[nJetTrees];
  
  Int_t nref_[nJetTrees];
  Float_t jtpt_[nJetTrees][nMaxJet];
  Float_t jtphi_[nJetTrees][nMaxJet];
  Float_t jteta_[nJetTrees][nMaxJet];
  
  std::vector<std::vector<float>* > jtpt_p, jtphi_p, jteta_p;
  jtpt_p.reserve(nJetTrees);
  jtphi_p.reserve(nJetTrees);
  jteta_p.reserve(nJetTrees);

  //std::cout << __LINE__ << std::endl;
  
  for(Int_t fI = 0; fI < nJetTrees; ++fI){
    jtpt_p.push_back(NULL);
    jtphi_p.push_back(NULL);
    jteta_p.push_back(NULL);

    jtpt_p.at(fI) = new std::vector<float>;
    jtphi_p.at(fI) = new std::vector<float>;
    jteta_p.at(fI) = new std::vector<float>;
  }

  //std::cout << __LINE__ << std::endl;
  


  for(unsigned int fI = 0; fI < inFileName.size(); ++fI){
    TFile* inFile_p = TFile::Open(mntToXRootdFileString(inFileName.at(fI)).c_str(), "READ");
    TTree* l1CaloTree_p = (TTree*)inFile_p->Get("l1CaloTowerEmuTree/L1CaloTowerTree");
    L1Analysis::L1AnalysisL1CaloTowerDataFormat *towers_ = new L1Analysis::L1AnalysisL1CaloTowerDataFormat();
    TTree* inL1AlgoUpgradeTree_p = (TTree*)inFile_p->Get("l1UpgradeEmuTree/L1UpgradeTree");
    L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
    TTree* inL1AlgoEvtTree_p = (TTree*)inFile_p->Get("l1EventTree/L1EventTree");
    L1Analysis::L1AnalysisEventDataFormat* evt = new L1Analysis::L1AnalysisEventDataFormat();
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    l1CaloTree_p->SetBranchStatus("*", 0);
    l1CaloTree_p->SetBranchStatus("L1CaloTower", 1);
    l1CaloTree_p->SetBranchStatus("iet", 1);
    l1CaloTree_p->SetBranchStatus("ieta", 1);
    l1CaloTree_p->SetBranchStatus("iphi", 1);
    l1CaloTree_p->SetBranchStatus("iem", 1);
    l1CaloTree_p->SetBranchStatus("ihad", 1);
    
    l1CaloTree_p->SetBranchAddress("L1CaloTower", &towers_);
    
    //std::cout << __LINE__ << std::endl;
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    inL1AlgoUpgradeTree_p->SetBranchStatus("*", 0);
    inL1AlgoUpgradeTree_p->SetBranchStatus("L1Upgrade", 1);
    inL1AlgoUpgradeTree_p->SetBranchStatus("nJets", 1);
    inL1AlgoUpgradeTree_p->SetBranchStatus("jetEt", 1);
    inL1AlgoUpgradeTree_p->SetBranchStatus("jetIEt", 1);
    inL1AlgoUpgradeTree_p->SetBranchStatus("jetEta", 1);
    inL1AlgoUpgradeTree_p->SetBranchStatus("jetIEta", 1);
    inL1AlgoUpgradeTree_p->SetBranchStatus("jetPhi", 1);
    inL1AlgoUpgradeTree_p->SetBranchStatus("jetIPhi", 1);
    inL1AlgoUpgradeTree_p->SetBranchStatus("jetSeedEt", 1);
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    inL1AlgoUpgradeTree_p->SetBranchAddress("L1Upgrade", &upgrade);
    
    inL1AlgoEvtTree_p->SetBranchStatus("*", 0);
    inL1AlgoEvtTree_p->SetBranchStatus("Event", 1);
    inL1AlgoEvtTree_p->SetBranchStatus("run", 1);
    inL1AlgoEvtTree_p->SetBranchStatus("lumi", 1);
    inL1AlgoEvtTree_p->SetBranchStatus("event", 1);
    
    inL1AlgoEvtTree_p->SetBranchAddress("Event", &evt);

    //std::cout << __LINE__ << std::endl;
    
    //std::cout << __LINE__ << std::endl;
  
    if(doForest){          
      //std::cout << __LINE__ << std::endl;

      hiTree_p = NULL;      
      forestFile_p = TFile::Open(mntToXRootdFileString(forestFileName.at(fI)).c_str(), "READ");
      
      //std::cout << __LINE__ << std::endl;
      if(hasHI) hiTree_p = (TTree*)forestFile_p->Get("hiEvtAnalyzer/HiTree");
      else if(hasRho) hiTree_p = (TTree*)forestFile_p->Get("rcRhoR4N11HiNtuplizer/EventTree");
      hiTree_p->SetBranchStatus("*", 0);
      hiTree_p->SetBranchStatus("rho", 1);
      
      hiTree_p->SetBranchAddress("rho", &rho_);
      
      //std::cout << __LINE__ << std::endl;
      
      for(Int_t jI = 0; jI < nJetTrees; ++jI){
	jetTrees_p[jI] = (TTree*)forestFile_p->Get(jetTrees.at(jI).c_str());
	
	jetTrees_p[jI]->SetBranchStatus("*", 0);

	//std::cout << __LINE__ << std::endl;
	
	if(doWeights) jetTrees_p[jI]->SetBranchStatus("pthat", 1);
	if(jetTrees.at(jI).find("/t") != std::string::npos) jetTrees_p[jI]->SetBranchStatus("nref", 1);
	jetTrees_p[jI]->SetBranchStatus("jtpt", 1);
	jetTrees_p[jI]->SetBranchStatus("jtphi", 1);
	jetTrees_p[jI]->SetBranchStatus("jteta", 1);
	
	//std::cout << __LINE__ << std::endl;

	if(doWeights) jetTrees_p[jI]->SetBranchAddress("pthat", &(pthat_[jI]));
	if(jetTrees.at(jI).find("/t") != std::string::npos){
	  jetTrees_p[jI]->SetBranchAddress("nref", &(nref_[jI]));
	  jetTrees_p[jI]->SetBranchAddress("jtpt", jtpt_[jI]);
	  jetTrees_p[jI]->SetBranchAddress("jtphi", jtphi_[jI]);
	  jetTrees_p[jI]->SetBranchAddress("jteta", jteta_[jI]);
	}
	else{	  

	  //std::cout << __LINE__ << std::endl;
	  jetTrees_p[jI]->SetBranchAddress("jtpt", &(jtpt_p.at(jI)));
	  jetTrees_p[jI]->SetBranchAddress("jtphi", &(jtphi_p.at(jI)));
	  jetTrees_p[jI]->SetBranchAddress("jteta", &(jteta_p.at(jI)));

	  //std::cout << __LINE__ << std::endl;

	}
      }
      //std::cout << __LINE__ << std::endl;
    }
      
    
    //std::cout << __LINE__ << std::endl;
    
    //    const Int_t nEntries = TMath::Min(1000000, (Int_t)l1CaloTree_p->GetEntries());
    const Int_t nEntries = TMath::Min(10000000, (Int_t)l1CaloTree_p->GetEntries());
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
                
    std::cout << "Processing " << nEntries << "..." << std::endl;
    for(Int_t entry = 0; entry < nEntries; ++entry){
      //      if(entry != 94) continue;

      if(entry%1000 == 0){
	std::cout << " Entry " << entry << "/" << nEntries << std::endl;

	std::cout << triggerFires[0][firstTypePos][0] << "/" << triggerOpp[0][firstTypePos][0] << std::endl;
      }
      
      //    if(entry != 17113 && entry != 21212) continue;

      inFile_p->cd();
      l1CaloTree_p->GetEntry(entry);
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      inL1AlgoUpgradeTree_p->GetEntry(entry);
      inL1AlgoEvtTree_p->GetEntry(entry);
      
      bool firstForestCall = false;
      Int_t forestEntry = -1;
      bool firstLeadingFill[nJetTrees];
      for(Int_t tI = 0; tI < nJetTrees; ++tI){
	firstLeadingFill[tI] = false;
      }

      Double_t weight = 1.;
      
      if(doForest){
	forestFile_p->cd();
	if(!firstForestCall){	  
	  forestEntry = forestMap[fI]->getEntryFromKey(evt->run, evt->lumi, evt->event);
	  firstForestCall = true;
	}
	
	if(forestEntry == -1){
	  //	std::cout << "Continuing on..." << evt->run << ", " << evt->lumi << ", " <<  evt->event << std::endl;
	  continue;
	}
	totalMatch++;
	  
	hiTree_p->GetEntry(forestEntry);
	
	for(Int_t tI = 0; tI < nJetTrees; ++tI){
	  jetTrees_p[tI]->GetEntry(forestEntry);
	}

	if(doWeights){
	  if(pthat_[0] > 80.) weight = 0.012747;
	  else weight = 1.;
	  
	  pthat_p->Fill(pthat_[0], 1.);
	  pthatWeighted_p->Fill(pthat_[0], weight);
	}
      }

      int leadJtIEt = -999;
      int leadJtIEta = -999;
      int leadJtIPhi = -999;
      int leadJtSeedPt = -999;
      
      //    if(evt->run != 1) continue;
      //    if(evt->lumi != 405) continue;
      //    if(evt->event != 40469) continue;
      
      for(unsigned int jI = 0; jI < upgrade->jetIEt.size(); ++jI){
	//      if(TMath::Abs(upgrade->jetIEta.at(jI)) >= gtEta(25)) continue;
	if(upgrade->jetIEt.at(jI) > leadJtIEt){
	  leadJtIEt = upgrade->jetIEt.at(jI);
	  leadJtIEta = upgrade->jetIEta.at(jI);
	  leadJtIPhi = upgrade->jetIPhi.at(jI);
	  leadJtSeedPt = upgrade->jetSeedEt.at(jI);
	}
      }
       
      //std::cout << __LINE__ << std::endl;
      
      if(leadJtIEt == -999) leadJtIEt = 0;
      
      std::vector<std::map<int, std::map<int, int> > > newTowIEt;
      newTowIEt.reserve(nSubType);
      std::vector<std::map<int, std::map<int, bool> > > newTowIEtUsed;
      newTowIEtUsed.reserve(nSubType);
      
      for(Int_t sI = 0; sI < nSubType; ++sI){
	std::map<int, std::map<int, int> > temp;
	std::map<int, std::map<int, bool> > tempUsed;
	(temp[-41])[0] = 0;
	(tempUsed[-41])[0] = false;
	newTowIEt.push_back(temp);
	newTowIEtUsed.push_back(tempUsed);
	
	for(int etaI = -41; etaI <= 41; ++etaI){
	  for(int phiI = 0; phiI <= nIPhi; ++phiI){
	    ((newTowIEt.at(sI))[etaI])[phiI] = 0;
	    ((newTowIEtUsed.at(sI))[etaI])[phiI] = false;
	  }
	}
      }
      
      std::vector<std::map<Int_t, Int_t> > uePerEta;
      uePerEta.reserve(nSubType);
      for(Int_t sI = 0; sI < nSubType; ++sI){
	std::map<Int_t, Int_t> temp;
	temp[1] = 0;
	uePerEta.push_back(temp);
	
	for(Int_t i = 1; i <= 41; ++i){
	  (uePerEta.at(sI))[i] = 0;
	  (uePerEta.at(sI))[-i] = 0;
	}
      }
      
      for(unsigned int i = 0; i < towers_->iet.size(); ++i){
	if(towers_->ieta[i] == -25) continue;
	if(towers_->ieta[i] == -26) continue;
	if(towers_->ieta[i] == -27) continue;
	if(towers_->ieta[i] == -28) continue;
	
	if(towers_->ieta[i] == 25) continue;
	if(towers_->ieta[i] == 26) continue;
	if(towers_->ieta[i] == 27) continue;
	if(towers_->ieta[i] == 28) continue;
	
	for(Int_t sI = 0; sI < nSubType; ++sI){
	  (uePerEta.at(sI))[towers_->ieta[i]] += towers_->iet[i];
	}
      }
      
      
      for(Int_t sI = 0; sI < nSubType; ++sI){
	for(Int_t etaI = 29; etaI < 41; ++etaI){
	  (uePerEta.at(sI))[-etaI] = (uePerEta.at(sI))[-(etaI+1)];
	  (uePerEta.at(sI))[etaI] = (uePerEta.at(sI))[etaI+1];
	}
	(uePerEta.at(sI))[-41] = 0;
	(uePerEta.at(sI))[41] = 0;
	
	
	for(Int_t phiI = 0; phiI <= nIPhi; ++phiI){	  	
	  (uePerEta.at(sI))[-28] = 0;
	  (uePerEta.at(sI))[-27] = 0;
	  (uePerEta.at(sI))[-26] = 0;
	  (uePerEta.at(sI))[-25] = 0;
	  
	  (uePerEta.at(sI))[28] = 0;
	  (uePerEta.at(sI))[27] = 0;
	  (uePerEta.at(sI))[26] = 0;
	  (uePerEta.at(sI))[25] = 0;	
	}
            
      
	bool isNone = subType[sI].size() == std::string("None").size() && subType[sI].find("None") != std::string::npos;
	bool isChunkyDonut = subType[sI].size() == std::string("ChunkyDonut").size() && subType[sI].find("ChunkyDonut") != std::string::npos;
	bool isChunkyDonutHIMod = subType[sI].size() == std::string("ChunkyDonutHIMod").size() && subType[sI].find("ChunkyDonutHIMod") != std::string::npos;
	bool isChunkyDonutLUT = subType[sI].size() == std::string("ChunkyDonutLUT").size() && subType[sI].find("ChunkyDonutLUT") != std::string::npos;
	bool isChunkyDonutZero = subType[sI].size() == std::string("ChunkyDonutZero").size() && subType[sI].find("ChunkyDonutZero") != std::string::npos;
	bool isPhiRingPPTower = subType[sI].size() == std::string("PhiRingPPTower").size() && subType[sI].find("PhiRingPPTower") != std::string::npos;
	bool isPhiRingPPTower72 = subType[sI].size() == std::string("PhiRingPPTower72").size() && subType[sI].find("PhiRingPPTower72") != std::string::npos;
	bool isPhiRingPPTower72NoDiv = subType[sI].size() == std::string("PhiRingPPTower72NoDiv").size() && subType[sI].find("PhiRingPPTower72NoDiv") != std::string::npos;
	bool isPhiRingPPTower72NoDivNoNeg = subType[sI].size() == std::string("PhiRingPPTower72NoDivNoNeg").size() && subType[sI].find("PhiRingPPTower72NoDivNoNeg") != std::string::npos;
	bool isPhiRingPPTower72NoDivCorr = subType[sI].size() == std::string("PhiRingPPTower72NoDivCorr").size() && subType[sI].find("PhiRingPPTower72NoDivCorr") != std::string::npos;
	bool isPhiRingHITower = subType[sI].size() == std::string("PhiRingHITower").size() && subType[sI].find("PhiRingHITower") != std::string::npos;
	bool isPhiRingHITower144 = subType[sI].size() == std::string("PhiRingHITower144").size() && subType[sI].find("PhiRingHITower144") != std::string::npos;
	bool isPhiRingHIRegion = subType[sI].size() == std::string("PhiRingHIRegion").size() && subType[sI].find("PhiRingHIRegion") != std::string::npos;
	bool isPhiRingHIRegion288 = subType[sI].size() == std::string("PhiRingHIRegion288").size() && subType[sI].find("PhiRingHIRegion288") != std::string::npos;
      
	bool doEvt = doType[sI] && (isNone || isChunkyDonut || isChunkyDonutHIMod || isChunkyDonutZero || isChunkyDonutLUT || isPhiRingPPTower || isPhiRingPPTower72NoDiv || isPhiRingPPTower72 || isPhiRingPPTower72NoDivNoNeg || isPhiRingPPTower72NoDivCorr || isPhiRingHITower144 || isPhiRingHITower || isPhiRingHIRegion || isPhiRingHIRegion288);

	if(!doEvt) continue;
	
	int multFact = 1;
	int threshMultFact = 1;
	if(isPhiRingPPTower72 || isPhiRingPPTower72NoDiv || isPhiRingPPTower72NoDivCorr || isPhiRingPPTower72NoDivNoNeg) multFact = 72;
	else if(isPhiRingHITower144) multFact = 144;
	else if(isPhiRingHIRegion288) multFact = 288;
	
	if(isPhiRingPPTower72NoDiv || isPhiRingPPTower72NoDivNoNeg) threshMultFact = 72;
	else if(isPhiRingPPTower72NoDivCorr) threshMultFact = (72 - 9);
	
	if(isNone || isChunkyDonut || isChunkyDonutLUT || isChunkyDonutZero || isChunkyDonutHIMod){
	  for(Int_t i = 1; i <= 41; ++i){
	    (uePerEta.at(sI))[i] = 0;
	    (uePerEta.at(sI))[-i] = 0;
	  }
	}
	else if(isPhiRingPPTower){
	  for(unsigned int i = 1; i < 41; ++i){
	    (uePerEta.at(sI))[-i] /= 72;
	    (uePerEta.at(sI))[i] /= 72;
	  }
	}
	else if(isPhiRingHITower || isPhiRingHITower144){
	  for(unsigned int i = 1; i < 41; ++i){
	    if(i%2 != 0) continue;
	    Float_t tempVal = ((uePerEta.at(sI))[i-1] + (uePerEta.at(sI))[i]);
	    if(isPhiRingHITower) tempVal /= 144;
	    (uePerEta.at(sI))[i-1] = tempVal;
	    (uePerEta.at(sI))[i] = tempVal;
	    
	    tempVal = ((uePerEta.at(sI))[-(i-1)] + (uePerEta.at(sI))[-i]);
	    if(isPhiRingHITower) tempVal /= 144;
	    (uePerEta.at(sI))[-(i-1)] = tempVal;
	    (uePerEta.at(sI))[-i] = tempVal;
	  }
	}
	else if(isPhiRingHIRegion || isPhiRingHIRegion288){
	  for(unsigned int i = 1; i < 41; ++i){
	    if(i%4 == 3) continue;
	    if(i%4 == 0) continue;
	    if(i%4 == 2) continue;
	    
	    Float_t tempVal = ((uePerEta.at(sI))[i] + (uePerEta.at(sI))[i+1] + (uePerEta.at(sI))[i+2] + (uePerEta.at(sI))[i+3]);
	    if(isPhiRingHIRegion) tempVal /= 288;
	    (uePerEta.at(sI))[i] = tempVal;
	    (uePerEta.at(sI))[i+1] = tempVal;
	    (uePerEta.at(sI))[i+2] = tempVal;
	    (uePerEta.at(sI))[i+3] = tempVal;
	    
	    tempVal = ((uePerEta.at(sI))[-i] + (uePerEta.at(sI))[-i-1] + (uePerEta.at(sI))[-i-2] + (uePerEta.at(sI))[-i-3]);
	    if(isPhiRingHIRegion) tempVal /= 288;
	    (uePerEta.at(sI))[-i] = tempVal;
	    (uePerEta.at(sI))[-i-1] = tempVal;
	    (uePerEta.at(sI))[-i-2] = tempVal;
	    (uePerEta.at(sI))[-i-3] = tempVal;
	  }
	}

    
      
	for(unsigned int i = 0; i < towers_->iet.size(); ++i){
	  int puEtaPos = towers_->ieta[i];
	  if(puEtaPos >= 29) puEtaPos--;
	  else if(puEtaPos <= -29) puEtaPos++;
	  
	  //	int iEt = TMath::Max(0, towers_->iet.at(i)*multFact - ((int)((uePerEta.at(sI))[puEtaPos]/72)));
	  int iEt = towers_->iet.at(i)*multFact;
	  if(isPhiRingPPTower72NoDivNoNeg) iEt += multFact*200;
	  
	  iEt -= ((int)((uePerEta.at(sI))[puEtaPos]));	

	  //	  if(isPhiRingPPTower72NoDiv || isPhiRingPPTower72NoDivCorr || isPhiRingPPTower72NoDivNoNeg) iEt -= 36;

	  if(isPhiRingPPTower72NoDivNoNeg && iEt < 0) iEt = 0;

	  ((newTowIEt.at(sI))[towers_->ieta[i]])[towers_->iphi[i]] = iEt;
	  ((newTowIEtUsed.at(sI))[towers_->ieta[i]])[towers_->iphi[i]] = true;
	}

	
	
	for(Int_t phiI = 0; phiI <= nIPhi; ++phiI){
	  for(Int_t etaI = -41; etaI <= 41; ++etaI){	    
	    
	    if(!((newTowIEtUsed.at(sI))[etaI])[phiI]){
	      
	      //	      if(isPhiRingPPTower72NoDivCorr) std::cout << "Orig: " << ((newTowIEt.at(sI))[etaI])[phiI] << ", " << etaI << ", " << phiI << std::endl;
	      int etTempVal = -(uePerEta.at(sI))[etaI];
	      //	      if(isPhiRingPPTower72NoDiv || isPhiRingPPTower72NoDivCorr || isPhiRingPPTower72NoDivNoNeg) etTempVal -= 36;;
	      if(isPhiRingPPTower72NoDivNoNeg){
		etTempVal += multFact*200;
		if(etTempVal < 0) etTempVal = 0;
	      }

	      ((newTowIEt.at(sI))[etaI])[phiI] = etTempVal;

	      /*	      
	      if(isPhiRingPPTower72NoDivCorr){
		std::cout << " Sub: " << (uePerEta.at(sI))[etaI] << ", " << etaI << std::endl;
		std::cout << " Sub 36..." << std::endl;
		std::cout << " Final: " << ((newTowIEt.at(sI))[etaI])[phiI] << std::endl;
	      }
	      */
	    }
	  }
	}     
      

	for(Int_t phiI = 0; phiI <= nIPhi; ++phiI){
	  for(Int_t etaI = 29; etaI < 41; ++etaI){
	    ((newTowIEt.at(sI))[-etaI])[phiI] = ((newTowIEt.at(sI))[-(etaI+1)])[phiI];
	    ((newTowIEt.at(sI))[etaI])[phiI] = ((newTowIEt.at(sI))[etaI+1])[phiI];
	  }
	  ((newTowIEt.at(sI))[-41])[phiI] = 0;
	  ((newTowIEt.at(sI))[41])[phiI] = 0;
	}
	for(Int_t phiI = 0; phiI <= nIPhi; ++phiI){
	  
	  if(!isChunkyDonut && !isChunkyDonutLUT && !isPhiRingPPTower72NoDivNoNeg){
	    ((newTowIEt.at(sI))[-28])[phiI] = 0;
	    ((newTowIEt.at(sI))[-27])[phiI] = 0;
	    ((newTowIEt.at(sI))[-26])[phiI] = 0;
	    ((newTowIEt.at(sI))[-25])[phiI] = 0;
	    
	    ((newTowIEt.at(sI))[28])[phiI] = 0;
	    ((newTowIEt.at(sI))[27])[phiI] = 0;
	    ((newTowIEt.at(sI))[26])[phiI] = 0;
	    ((newTowIEt.at(sI))[25])[phiI] = 0;		
	  }
	  else if(isPhiRingPPTower72NoDivNoNeg){
	    ((newTowIEt.at(sI))[-28])[phiI] = multFact*200;
            ((newTowIEt.at(sI))[-27])[phiI] = multFact*200;
            ((newTowIEt.at(sI))[-26])[phiI] = multFact*200;
            ((newTowIEt.at(sI))[-25])[phiI] = multFact*200;

            ((newTowIEt.at(sI))[28])[phiI] = multFact*200;
            ((newTowIEt.at(sI))[27])[phiI] = multFact*200;
            ((newTowIEt.at(sI))[26])[phiI] = multFact*200;
            ((newTowIEt.at(sI))[25])[phiI] = multFact*200;
	  }
	}
	
	

	/*      
	const std::string csvStr = "output/l1Towers_" + subType[sI] + "_Run" + std::to_string(evt->run) + "_Lumi" + std::to_string(evt->lumi) + "_Evt" + std::to_string(evt->event) + "_" + dateStr + ".csv";
	
	std::cout << csvStr << std::endl;
	std::ofstream file(csvStr.c_str());
	
	
	std::string disp = "  ,";
	file << disp;
	//std::cout << disp;
	for(int etaI = -41; etaI <= 41; ++etaI){
	  disp = std::to_string(TMath::Abs(etaI)) + ",";
	  if(disp.size() == 2) disp = " " + disp;
	  file << disp;
	  //std::cout << disp;
	}
	file << std::endl;
	//std::cout << std::endl;
	
	for(int phiI = 0; phiI <= nIPhi; ++phiI){	
	  disp = std::to_string(phiI) + ",";
	  if(disp.size() == 2) disp = " " + disp;
	  file << disp;
	  //std::cout << disp;
	  
	  for(int etaI = -41; etaI <= 41; ++etaI){
	    disp = std::to_string(((newTowIEt.at(sI))[etaI])[phiI]);
	    if(disp.size() == 1) disp = " " + disp + ",";
	    else if(disp.size() == 2) disp = disp + ",";
	    
	    file << disp;
	    //std::cout << disp;
	  }
	  file << std::endl;
	  //std::cout << std::endl;
	}    
	file.close();
	//std::cout << __LINE__ << std::endl;
	*/
	
	std::vector<int> jetPt;
	std::vector<int> jetPhi;
	std::vector<int> jetEta;
	std::vector<int> jetPhiPrev;
	std::vector<int> jetEtaPrev;
	std::vector<int> seedPt;
	
	for(int etaI = -41; etaI <= 41; ++etaI){
	  
	  if(etaI == 0) continue;
	  
	  if(!isChunkyDonut && !isChunkyDonutLUT){
	    if(etaI == -25) continue;
	    if(etaI == -26) continue;
	    if(etaI == -27) continue;
	    if(etaI == -28) continue;
	    
	    if(etaI == 25) continue;
	    if(etaI == 26) continue;
	    if(etaI == 27) continue;
	    if(etaI == 28) continue;
	  }
	  
	  for(int phiI = 0; phiI <= nIPhi; ++phiI){
	    if(((newTowIEt.at(sI))[etaI])[phiI] >= minSeedThresh*multFact){
	      int tempJetPt = 0;
	      bool goodSeed = true;
	      
	      for(Int_t etaI2 = TMath::Max(-41, etaI - 4); etaI2 <= TMath::Min(41, etaI + 4); ++etaI2){	
		int etaPos2 = etaI2;
		if(etaI > 0 && etaPos2 <= 0) etaPos2--;
		if(etaI < 0 && etaPos2 >= 0) etaPos2++;
		
		if(!isChunkyDonut && !isChunkyDonutLUT){
		  if(etaPos2 == -25) continue;
		  if(etaPos2 == -26) continue;
		  if(etaPos2 == -27) continue;
		  if(etaPos2 == -28) continue;
		  
		  if(etaPos2 == 25) continue;
		  if(etaPos2 == 26) continue;
		  if(etaPos2 == 27) continue;
		  if(etaPos2 == 28) continue;
		}
		
		for(Int_t phiI2 = phiI - 4; phiI2 <= phiI + 4; ++phiI2){ 	      
		  int phiPos = phiI2;
		  if(phiPos <= 0) phiPos += nIPhi;
		  if(phiPos > nIPhi) phiPos -= nIPhi;
		  
		  tempJetPt += ((newTowIEt.at(sI))[etaPos2])[phiPos];
		  
		  if(etaPos2 == etaI && phiPos == phiI) continue;
		  else{
		    int dPhi = phiI2 - phiI + 4;
		    int dEta = etaI2 - etaI + 4;
		    
		    if(mask_[8-dPhi][dEta] == 1 && ((newTowIEt.at(sI))[etaI])[phiI] < (newTowIEt.at(sI))[etaPos2][phiPos]) goodSeed = false;
		    else if(mask_[8-dPhi][dEta] == 2 && ((newTowIEt.at(sI))[etaI])[phiI] <= ((newTowIEt.at(sI))[etaPos2])[phiPos]) goodSeed = false;
		    
		    if(!goodSeed){
		      break;
		    }
		  }
		}
		if(!goodSeed) break;
	      }
	      
	      if(!goodSeed) continue;
	      
	      if(isChunkyDonut || isChunkyDonutLUT || isChunkyDonutZero || isChunkyDonutHIMod){
		int ptFlapPhiUp = 0;
		int ptFlapPhiDown = 0;
		int ptFlapEtaUp = 0;
		int ptFlapEtaDown = 0;
		
		for(Int_t etaI2 = TMath::Max(-41, etaI - 4); etaI2 <= TMath::Min(41, etaI + 4); ++etaI2){	
		  Int_t etaPos = etaI2;
		  if(etaI < 0 && etaPos >= 0) etaPos++;
		  if(etaI > 0 && etaPos <= 0) etaPos--;
		  
		  for(Int_t phiI2 = phiI + 5; phiI2 <= phiI + 7; ++phiI2){ 	      
		    int phiPos = phiI2;
		    if(phiPos <= 0) phiPos += nIPhi;
		    if(phiPos > nIPhi) phiPos -= nIPhi;
		    
		    ptFlapPhiUp += ((newTowIEt.at(sI))[etaPos])[phiPos];
		  }
		}
		
		for(Int_t etaI2 = TMath::Max(etaI - 4, -41); etaI2 <= TMath::Min(41, etaI + 4); ++etaI2){	
		  Int_t etaPos = etaI2;
		  if(etaI < 0 && etaPos >= 0) etaPos++;
		  if(etaI > 0 && etaPos <= 0) etaPos--;
		  
		  for(Int_t phiI2 = phiI - 7; phiI2 <= phiI - 5; ++phiI2){ 	      
		    int phiPos = phiI2;
		    if(phiPos <= 0) phiPos += nIPhi;
		    if(phiPos > nIPhi) phiPos -= nIPhi;
		    
		    ptFlapPhiDown += ((newTowIEt.at(sI))[etaPos])[phiPos];
		  }
		}
		
		for(Int_t phiI2 = phiI - 4; phiI2 <= phiI + 4; ++phiI2){ 	      
		  int phiPos = phiI2;
		  if(phiPos <= 0) phiPos += nIPhi;
		  if(phiPos > nIPhi) phiPos -= nIPhi;
		  
		  for(Int_t etaI2 = TMath::Max(-41, etaI - 7); etaI2 <= TMath::Max(-41, etaI - 5); ++etaI2){
		    Int_t etaPos = etaI2;
		    if(etaI < 0 && etaPos >= 0) etaPos++;
		    if(etaI > 0 && etaPos <= 0) etaPos--;
		    
		    ptFlapEtaDown += ((newTowIEt.at(sI))[etaPos])[phiPos];
		  }
		}
		
		for(Int_t phiI2 = phiI - 4; phiI2 <= phiI + 4; ++phiI2){ 	      
		  int phiPos = phiI2;
		  if(phiPos <= 0) phiPos += nIPhi;
		  if(phiPos > nIPhi) phiPos -= nIPhi;
		  
		  for(Int_t etaI2 = TMath::Min(41, etaI + 5); etaI2 <= TMath::Min(41, etaI + 7); ++etaI2){
		    Int_t etaPos = etaI2;
		    if(etaI < 0 && etaPos >= 0) etaPos++;
		    if(etaI > 0 && etaPos <= 0) etaPos--;
		    
		    ptFlapEtaUp += ((newTowIEt.at(sI))[etaPos])[phiPos];
		  }
		}
		
		std::vector<int> flaps = {ptFlapEtaUp, ptFlapEtaDown, ptFlapPhiUp, ptFlapPhiDown};
		std::sort(std::begin(flaps), std::end(flaps));
		
		if(isChunkyDonut || isChunkyDonutLUT || isChunkyDonutZero || TMath::Abs(etaI) < 30) tempJetPt -= (flaps.at(0) + flaps.at(1) + flaps.at(2));
		else tempJetPt -= (flaps.at(2) + flaps.at(2) + flaps.at(3));
		
		if(tempJetPt < 0) tempJetPt = 0;
		if(tempJetPt > 0){
		  //	      std::cout << " Flaps: " << flaps.at(0) << ", " << flaps.at(1) << ", " << flaps.at(2) << ", " << flaps.at(3) << std::endl;
		}
	      }
	      
	      if((newTowIEt.at(sI))[etaI][phiI] == 510 || (newTowIEt.at(sI))[etaI][phiI] == 509 || (newTowIEt.at(sI))[etaI][phiI] == 511) tempJetPt = 2047;
	      
	      if(isPhiRingPPTower72 || isPhiRingHITower144 || isPhiRingHIRegion288) tempJetPt /= multFact;
	      
	      if(tempJetPt <= 0) continue;
	      
	      if(TMath::Abs(etaI) < 25 || TMath::Abs(etaI) > 29 || !isChunkyDonutLUT){
		if(isPhiRingPPTower72NoDivCorr){
		  Double_t corrPt = 0;
	      
		  for(Int_t etaI2 = TMath::Max(-41, etaI - 4); etaI2 <= TMath::Min(41, etaI + 4); ++etaI2){	
		    int etaPos2 = etaI2;
		    if(etaI > 0 && etaPos2 <= 0) etaPos2--;
		    if(etaI < 0 && etaPos2 >= 0) etaPos2++;
		    
		    if(etaPos2 == -25) continue;
		    if(etaPos2 == -26) continue;
		    if(etaPos2 == -27) continue;
		    if(etaPos2 == -28) continue;
		    
		    if(etaPos2 == 25) continue;
		    if(etaPos2 == 26) continue;
		    if(etaPos2 == 27) continue;
		    if(etaPos2 == 28) continue;

		    Int_t origNineTow = 0;

		    int puEtaPos = etaPos2;
		    //		    if(puEtaPos >= 29) puEtaPos--;
		    //		    else if(puEtaPos <= -29) puEtaPos++;
		
		    for(Int_t phiI2 = phiI - 4; phiI2 <= phiI + 4; ++phiI2){ 	      
		      int phiPos = phiI2;
		      if(phiPos <= 0) phiPos += nIPhi;
		      if(phiPos > nIPhi) phiPos -= nIPhi;

		      if(isPhiRingPPTower72NoDivCorr){
			//			std::cout << "Compare " << ((newTowIEt.at(sI))[etaPos2])[phiPos] << ", " << etaPos2 << ", " << phiPos << std::endl;
			//			std::cout << " Add" << ((int)((uePerEta.at(sI))[puEtaPos])) << std::endl;
			//			std::cout << " Add 36..." << std::endl;
			//			std::cout << " Final " << ((newTowIEt.at(sI))[etaPos2])[phiPos] + 36 + ((int)((uePerEta.at(sI))[puEtaPos])) << std::endl;
		    
			if((((newTowIEt.at(sI))[etaPos2])[phiPos] + 36 + ((int)((uePerEta.at(sI))[puEtaPos])))%72 != 0){
			  //			  std::cout << "UH OH ERROR " << etaPos2 << ", " << phiPos << std::endl;
			}
		      }

		      //		      origNineTow += ((newTowIEt.at(sI))[etaPos2])[phiPos] + 36((int)((uePerEta.at(sI))[puEtaPos]));		      
		    }
		    origNineTow /= 72;
		    
		    Double_t tempUE = (uePerEta.at(sI))[etaPos2] - origNineTow;
		    
		    origNineTow *= (72 - 9);
		    origNineTow -= 9*tempUE - 9*31;
		    corrPt += origNineTow;
		  }
		  tempJetPt = corrPt;
		}
	      
		if(isPhiRingPPTower72NoDivNoNeg){
		  //		  std::cout << "Check: " << tempJetPt;

		  int rowsExclude = 0;
		  if(TMath::Abs(etaI) == 24 || TMath::Abs(etaI) == 29) rowsExclude = 4;
		  if(TMath::Abs(etaI) == 23 || TMath::Abs(etaI) == 30) rowsExclude = 3;
		  if(TMath::Abs(etaI) == 22 || TMath::Abs(etaI) == 31) rowsExclude = 2;
		  if(TMath::Abs(etaI) == 21 || TMath::Abs(etaI) == 32) rowsExclude = 1;

		  tempJetPt -= (81 - 9*rowsExclude)*multFact*200;
		  //		  std::cout << ", " << tempJetPt << std::endl;
		  if(tempJetPt < 0) continue;
		}

		jetPt.push_back(tempJetPt);
		jetPhi.push_back(phiI);
		jetEta.push_back(etaI);
		jetPhiPrev.push_back(phiI);
		jetEtaPrev.push_back(etaI);
		seedPt.push_back(((newTowIEt.at(sI))[etaI])[phiI]/multFact);
	      }
	    }
	  }
	}
	
	//  if(jetPt.size() == 0) continue;
	
	//std::cout << __LINE__ << std::endl;
	
	if(jetPt.size() > 0){
	  for(unsigned int jI = 0; jI < jetPt.size() - 1; ++jI){
	    for(unsigned int jI2 = jI+1; jI2 < jetPt.size(); ++jI2){
	      if(jetPt.at(jI) < jetPt.at(jI2)){
		int tempPt = jetPt.at(jI);
		int tempPhi = jetPhi.at(jI);
		int tempEta = jetEta.at(jI);
		int tempPhiPrev = jetPhiPrev.at(jI);
		int tempEtaPrev = jetEtaPrev.at(jI);
		int tempSeed = seedPt.at(jI);
		
		jetPt.at(jI) = jetPt.at(jI2);
		jetPhi.at(jI) = jetPhi.at(jI2);
		jetEta.at(jI) = jetEta.at(jI2);
		jetPhiPrev.at(jI) = jetPhiPrev.at(jI2);
		jetEtaPrev.at(jI) = jetEtaPrev.at(jI2);
		seedPt.at(jI) = seedPt.at(jI2);
	      
		jetPt.at(jI2) = tempPt;
		jetPhi.at(jI2) = tempPhi;
		jetEta.at(jI2) = tempEta;
		jetPhiPrev.at(jI2) = tempPhiPrev;
		jetEtaPrev.at(jI2) = tempEtaPrev;
		seedPt.at(jI2) = tempSeed;
	      } 
	    }
	  }
	}
	
	if(isFileType[sI]){
	  for(unsigned int jI = 0; jI < upgrade->jetIEt.size(); ++jI){
	    bool isJtFound = false;
	    
	    for(unsigned int jI2 = 0; jI2 < jetPt.size(); ++jI2){
	      if(jetPt.at(jI2) != upgrade->jetIEt.at(jI)) continue;
	      if(gtPhi(jetPhi.at(jI2)) != upgrade->jetIPhi.at(jI)) continue;
	      if(gtEta(jetEta.at(jI2)) != upgrade->jetIEta.at(jI)) continue;
	      
	      isJtFound = true;
	      break;
	    }
	    
	    if(!isJtFound){
	      std::cout << "Warning in event: " << evt->event << ", entry " << entry << std::endl;
	      std::cout << " Cannot find jet (pt, phi, eta): " << upgrade->jetIEt.at(jI) << ", " << upgrade->jetIPhi.at(jI) << ", " << upgrade->jetIEta.at(jI) << std::endl;
	      std::cout << " Candidates were (using sub type " << subType[sI] << "): " << std::endl;
	      for(unsigned int jI2 = 0; jI2 < jetPt.size(); ++jI2){
		std::cout << "  " << jetPt.at(jI2) << ", " << jetPhi.at(jI2) << ", " << jetEta.at(jI2) << std::endl;
	      }
	    }
	  }
	}
     

	for(Int_t lI = 0; lI < nL1Thresh; ++lI){
	  triggerOpp[nRhoBins][sI][lI]++;
	  
	  if(jetPt.size() > 0){
	    if(jetPt.at(0) > l1ThreshPt[lI]*threshMultFact){
	      triggerFires[nRhoBins][sI][lI]++;	  
	    }
	    
	    /*       
		     for(unsigned int jI = 0; jI < TMath::Min(1, jetPt.size()); ++jI){
		     if(jetPt.at(jI) > l1ThreshPt[lI]) triggerEta_L1Pt_h[nRhoBins][sI][lI]->Fill(jetEta.at(jI));
		     }
	    */
	    
	  }
	}
	
	
	
	if(doForest){
	  Int_t centPos = -1;
	  for(Int_t cI = 0; cI < nRhoBins; ++cI){
	    if(rhoBinsLow[cI] <= rho_ && rho_ < rhoBinsHi[cI]){
	      centPos = cI;
	      break;
	    }
	  }
	  
	  if(centPos == -1){
	    std::cout << "Warning: " << centPos << " for rho " << rho_ << std::endl;
	    continue;
	  }    
	  
	  for(Int_t lI = 0; lI < nL1Thresh; ++lI){
	    triggerOpp[centPos][sI][lI]++;
	    
	    if(jetPt.size() > 0){
	      if(jetPt.at(0) > l1ThreshPt[lI]*threshMultFact){
		triggerFires[centPos][sI][lI]++;	  
	      }
	    }
	  }
	  
	  for(Int_t lI = 0; lI < nL1Thresh; ++lI){
	    for(int jI = 0; jI < TMath::Min(10000, (int)jetPt.size()); ++jI){
	      if(jetPt.at(jI) > l1ThreshPt[lI]*threshMultFact){
		triggerEta_L1Pt_h[nRhoBins][sI][lI]->Fill(jetEta.at(jI), weight);
		if(centPos >= 0) triggerEta_L1Pt_h[centPos][sI][lI]->Fill(jetEta.at(jI), weight);
	      }
	    }
	  }
	  
	  /*
	    for(Int_t lI = 0; lI < nL1Thresh; ++lI){
	    triggerOpp[nRhoBins][sI][lI]++;
	    triggerOpp[centPos][sI][lI]++;
	    
	    if(jetPt.at(0) > l1ThreshPt[lI]){
	    triggerFires[nRhoBins][sI][lI]++;	  
	    triggerFires[centPos][sI][lI]++;
	    }
	    }
	  */
	  
	  
	  
	  for(Int_t tI = 0; tI < nJetTrees; ++tI){	    
	    Float_t tempLeadingJtPt_ = -999;
	    Float_t tempLeadingJtPhi_ = -999;
	    Float_t tempLeadingJtEta_ = -999;
	    
	    if(jetTrees.at(tI).find("/t") != std::string::npos){
	      for(Int_t jI = 0; jI < nref_[tI]; ++jI){
		if(TMath::Abs(jteta_[tI][jI]) > zeroedEtaLow && TMath::Abs(jteta_[tI][jI]) < zeroedEtaHi) continue;
		
		for(Int_t lI = 0; lI < nL1Thresh; ++lI){
		  if(jtpt_[tI][jI] > l1ThreshPt[lI]*threshMultFact){
		    jetEta_L1Pt_h[nRhoBins][tI][lI]->Fill(jteta_[tI][jI], weight);
		    if(centPos >= 0) jetEta_L1Pt_h[centPos][tI][lI]->Fill(jteta_[tI][jI], weight);
		  }
		}
		
		if(jtpt_[tI][jI] > tempLeadingJtPt_){
		  tempLeadingJtPt_ = jtpt_[tI][jI];
		  tempLeadingJtPhi_ = jtphi_[tI][jI];
		  tempLeadingJtEta_ = jteta_[tI][jI];
		}
	      }
	    }
	    else{
	      for(unsigned int jI = 0; jI < jtpt_p.at(tI)->size(); ++jI){
		if(TMath::Abs(jteta_p.at(tI)->at(jI)) > zeroedEtaLow && TMath::Abs(jteta_p.at(tI)->at(jI)) < zeroedEtaHi) continue;
		
		for(Int_t lI = 0; lI < nL1Thresh; ++lI){
		  if(jtpt_p.at(tI)->at(jI) > l1ThreshPt[lI]*threshMultFact){
		    jetEta_L1Pt_h[nRhoBins][tI][lI]->Fill(jteta_p.at(tI)->at(jI), weight);
		    if(centPos >= 0) jetEta_L1Pt_h[centPos][tI][lI]->Fill(jteta_p.at(tI)->at(jI), weight);
		  }
		}
		
		if(jtpt_p.at(tI)->at(jI) > tempLeadingJtPt_){
		  tempLeadingJtPt_ = jtpt_p.at(tI)->at(jI);
		  tempLeadingJtPhi_ = jtphi_p.at(tI)->at(jI);
		  tempLeadingJtEta_ = jteta_p.at(tI)->at(jI);
		}
	      }
	    }
	    
	    Int_t jtAbsEtaPos = -1;
	    for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	      if(TMath::Abs(tempLeadingJtEta_) >= jtAbsEtaBinsLow[aI] && TMath::Abs(tempLeadingJtEta_) < jtAbsEtaBinsHi[aI]){
		jtAbsEtaPos = aI;
		break;
	      }
	    }
	    
	    if(tempLeadingJtPt_ >= jtPtHi) tempLeadingJtPt_ = jtPtHi - 1;
	    
	    if(tempLeadingJtPt_ >= jtPtLow && tempLeadingJtPt_ < jtPtHi){
	      if(!firstLeadingFill[tI]){
		if(jtAbsEtaPos >= 0) leadingJetPt_h[tI][centPos][jtAbsEtaPos]->Fill(tempLeadingJtPt_, weight);
		leadingJetPt_h[tI][centPos][nJtAbsEtaBins]->Fill(tempLeadingJtPt_, weight);
		firstLeadingFill[tI] = true;
	      }
	      
	      Float_t leadL1Pt = -999;
	      Float_t leadL1PtDR = -999;
	      
	      for(unsigned int lI = 0; lI < jetPt.size(); ++lI){
		if(jetPt.at(lI) > leadL1Pt) leadL1Pt = jetPt.at(lI);
		
		if(getDR(tempLeadingJtEta_, tempLeadingJtPhi_, towerEta(jetEta.at(lI)), towerPhi(jetPhi.at(lI))) < maxDR && jetPt.at(lI) > leadL1PtDR) leadL1PtDR = jetPt.at(lI);
	      }
	      
	      for(Int_t lI = 0; lI < nL1Thresh; ++lI){
		if(leadL1Pt > l1ThreshPt[lI]*threshMultFact){
		  if(jtAbsEtaPos >= 0) leadingJetPt_L1Pt_h[tI][centPos][jtAbsEtaPos][sI][lI]->Fill(tempLeadingJtPt_, weight);
		  leadingJetPt_L1Pt_h[tI][centPos][nJtAbsEtaBins][sI][lI]->Fill(tempLeadingJtPt_, weight);
		}
		
		
		Int_t jtPos = -1;
		for(Int_t rI = 0; rI < nL1ResPtBins; ++rI){
		  if(l1ResPtBins[rI] <= tempLeadingJtPt_ && tempLeadingJtPt_ < l1ResPtBins[rI+1]){
		    jtPos = rI;
		    break;
		  }
		}
		
		if(jtPos >= 0){
		  if(leadL1PtDR <= 0) leadL1PtDR = 0.;
		  jetL1PtOverOffline_h[tI][centPos][sI][jtPos]->Fill(leadL1PtDR/tempLeadingJtPt_, weight);
		  jetL1PtOverOffline_h[tI][nRhoBins][sI][jtPos]->Fill(leadL1PtDR/tempLeadingJtPt_, weight);
		}
		
		if(leadL1PtDR > l1ThreshPt[lI]*threshMultFact){
		  if(jtAbsEtaPos >= 0) leadingJetPt_L1Pt_DR_h[tI][centPos][jtAbsEtaPos][sI][lI]->Fill(tempLeadingJtPt_, weight);
		  leadingJetPt_L1Pt_DR_h[tI][centPos][nJtAbsEtaBins][sI][lI]->Fill(tempLeadingJtPt_, weight);
		  
		  if(false/*tempLeadingJtPt_ < 30 && l1ThreshPt[lI] == 64 && jetTrees.at(tI).find("Calo") != std::string::npos*/){
		    std::cout << "Found < 30 jet with l1 64 in entry " << entry << std::endl;
		    std::cout << " run, lumi, evt: " << evt->run << ", " << evt->lumi << ", " << evt->event << std::endl;
		    std::cout << " jetpt, phi, eta: " << tempLeadingJtPt_ << ", " << tempLeadingJtPhi_ << ", " << tempLeadingJtEta_ << std::endl;
		    std::cout << "  options (pt, phi, eta): " << std::endl;
		    for(unsigned int lI2 = 0; lI2 < jetPt.size(); ++lI2){
		      std::cout << "   " << jetPt.at(lI2) << ", " << towerPhi(jetPhi.at(lI2)) << ", " << towerEta(jetEta.at(lI2)) << std::endl;
		      std::cout << "    " << jetPhi.at(lI2) << ", " << jetEta.at(lI2) << std::endl;
		    }
		    
		  }
		}
		else if(((Int_t)l1ThreshPt[lI]) == 24 && tempLeadingJtPt_ > 100 && jetTrees.at(tI).find("Calo") != std::string::npos){
		  std::cout << "Missed jet in entry, rhoPos " << entry << ", " << rhoBinsLow[centPos] << "-" << rhoBinsHi[centPos] << std::endl;
		  std::cout << " run, lumi, evt: " << evt->run << ", " << evt->lumi << ", " << evt->event << std::endl;
		  std::cout << " jetpt, phi, eta: " << tempLeadingJtPt_ << ", " << tempLeadingJtPhi_ << ", " << tempLeadingJtEta_ << std::endl;
		  //		std::cout << "  options (pt, phi, eta): " << std::endl;
		  for(unsigned int lI2 = 0; lI2 < jetPt.size(); ++lI2){
		    //		  std::cout << "   " << jetPt.at(lI2) << ", " << towerPhi(jetPhi.at(lI2)) << ", " << towerEta(jetEta.at(lI2)) << std::endl;
		    //		  std::cout << "    " << jetPhi.at(lI2) << ", " << jetEta.at(lI2) << std::endl;
		  }
		}
	      }
	    }
	  }
	} 
	
	
	//std::cout << __LINE__ << std::endl;
	
	if(jetPt.size() > 0){
	  if(leadJtIEt != jetPt.at(0)){
	    missCount++;
	    //      std::cout << "Lead iet, iphi, ieta, seedEt: " << leadJtIEt << ", " << leadJtIPhi << ", " << leadJtIEta << ", " << leadJtSeedPt << std::endl;
	    //      std::cout << "NJets: " << jetPt.size() << std::endl;
	    for(int i = 0; i < TMath::Min(8, (int)jetPt.size()); ++i){
	      //	std::cout << " " << i << "/" << jetPt.size() << ", pt, phi, eta, seed): " << jetPt.at(i) << ", " << jetPhi.at(i) << ", " << jetEta.at(i) << ", " << seedPt.at(i) << " (" << jetPhiPrev.at(i) << ", " << jetEtaPrev.at(i) << ")" << std::endl;
	    }
	  }
	  else foundCount++;
	  
	  if(leadJtIEt - jetPt.at(0) != 0 || jetPhi.at(0) != leadJtIPhi || jetEta.at(0) != leadJtIEta){
	    etaMiss_p->Fill(jetEta.at(0), weight);
	    phiMiss_p->Fill(jetPhi.at(0), weight);
	    
	    etaMiss_CMSSWJet_p->Fill(leadJtIEta, weight);
	    phiMiss_CMSSWJet_p->Fill(leadJtIPhi, weight);
	  }
	}
	//std::cout << __LINE__ << std::endl;

      }
    }

	//std::cout << __LINE__ << std::endl;
 
    inFile_p->Close();
    delete inFile_p;

    if(doForest){
      forestFile_p->Close();
      delete forestFile_p;
    }
  }


  
  //std::cout << __LINE__ << std::endl;

  outFile_p->cd();
  TH1F* triggerFrac_h[nRhoBins+1][nSubType];

  std::cout << "Total match: " << totalMatch << std::endl;

  std::cout << " With nRhoBins = " << nRhoBins << std::endl;
  for(Int_t rI = 0; rI < nRhoBins+1; ++rI){

    std::string centStr = "Rho" + std::to_string(rhoBinsLow[nRhoBins-1]) + "to" + std::to_string(rhoBinsHi[0]);
    std::string centStr2 =  std::to_string(rhoBinsLow[nRhoBins-1]) + " < #rho < " + std::to_string(rhoBinsHi[0]);
    if(rI < nRhoBins ){
      centStr = "Rho" + std::to_string(rhoBinsLow[rI]) + "to" + std::to_string(rhoBinsHi[rI]);
      centStr2 = std::to_string(rhoBinsLow[rI]) + " < #rho < " + std::to_string(rhoBinsHi[rI]);
    }

    if(rI < nRhoBins) std::cout << "Rho: " << prettyString(rhoBinsLow[rI], 3, false) << "-" << prettyString(rhoBinsHi[rI], 3, false) << std::endl;
    else std::cout << "All rho" << std::endl;

    for(Int_t sI = 0; sI < nSubType; ++sI){
      if(!doType[sI]) continue;

      triggerFrac_h[rI][sI] = new TH1F(("triggerFrac_" + centStr + "_" + subType[sI] + "_h").c_str(), ";L1 p_{T} Threshold;Fraction", nL1Thresh, l1ThreshPtBins);

      std::cout << " Sub type: " << subType[sI] << std::endl;
	
      for(Int_t lI = 0; lI < nL1Thresh; ++lI){
	if(triggerOpp[rI][sI][lI] == 0) continue;

	Double_t val1 = triggerFires[rI][sI][lI];
	Double_t val2 = triggerOpp[rI][sI][lI];
	triggerFrac_h[rI][sI]->SetBinContent(lI+1, val1/val2);
	triggerFrac_h[rI][sI]->SetBinError(lI+1, TMath::Sqrt(val1)/val2);

	std::cout << "  " << l1ThreshPt[lI] << ": " << triggerFires[rI][sI][lI] << "/" << triggerOpp[rI][sI][lI] << " ( " << prettyString(((Double_t)triggerFires[rI][sI][lI])/((Double_t)triggerOpp[rI][sI][lI]), 3, false) << ")" << std::endl;
      }
      
      triggerFrac_h[rI][sI]->Write("", TObject::kOverwrite);
    }

    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetLeftMargin(0.14);
    canv_p->SetBottomMargin(0.14);

    bool isDrawn = false;
    for(Int_t sI = 0; sI < nSubType; ++sI){
      if(!doType[sI]) continue;

      triggerFrac_h[rI][sI]->SetMarkerStyle(styles[sI]);
      triggerFrac_h[rI][sI]->SetMarkerColor(vg.getColor(sI));
      triggerFrac_h[rI][sI]->SetLineColor(vg.getColor(sI));
      triggerFrac_h[rI][sI]->SetMarkerSize(1.);
      
      centerTitles(triggerFrac_h[rI][sI]);

      triggerFrac_h[rI][sI]->SetMaximum(4.);
      triggerFrac_h[rI][sI]->SetMinimum(0.0001);

      if(!isDrawn){
	triggerFrac_h[rI][sI]->DrawCopy("E1 P HIST");
	isDrawn = true;
      }
      else triggerFrac_h[rI][sI]->DrawCopy("E1 P HIST SAME");
    }

    label_p->DrawLatex(0.2, 0.94, centStr2.c_str());
    dummyLeg3_p->Draw("SAME");
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    canv_p->SaveAs(("pdfDir/" + dateStr  + "/triggerFrac_" + centStr + "_" + extraTag + "_" + dateStr + ".pdf").c_str());

    delete canv_p;

    for(Int_t sI = 0; sI < nSubType; ++sI){
      if(!doType[sI]) continue;
      delete triggerFrac_h[rI][sI];
    }
  }

  std::cout << "End tower display" << std::endl;
  
  std::cout << "Found: " << foundCount << std::endl;
  std::cout << "Miss: " << missCount << std::endl;

  std::cout << "Processing complete!" << std::endl;
  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  outFile_p->cd();

  TCanvas* pthatCanv_p = new TCanvas("pthatCanv_p", "", 450, 450);
  pthatCanv_p->SetTopMargin(0.01);
  pthatCanv_p->SetRightMargin(0.01);
  pthatCanv_p->SetLeftMargin(0.14);
  pthatCanv_p->SetBottomMargin(0.14);

  pthat_p->SetMarkerStyle(styles[0]);
  pthat_p->SetMarkerColor(vg.getColor(0));
  pthat_p->SetLineColor(vg.getColor(0));
  pthat_p->SetMarkerSize(1);

  pthatWeighted_p->SetMarkerStyle(styles[1]);
  pthatWeighted_p->SetMarkerColor(vg.getColor(1));
  pthatWeighted_p->SetLineColor(vg.getColor(1));
  pthatWeighted_p->SetMarkerSize(1);

  pthatCanv_p->cd();
  pthat_p->SetMaximum(pthat_p->GetMaximum()*10.);

  Double_t min = 100000;
  for(Int_t pI = 0; pI < pthatWeighted_p->GetNbinsX(); ++pI){
    if(pthatWeighted_p->GetBinContent(pI+1) > 0 && pthatWeighted_p->GetBinContent(pI+1) < min) min = pthatWeighted_p->GetBinContent(pI+1);
  }

  pthat_p->SetMinimum(min/10.);

  pthat_p->DrawCopy("HIST E1 P");
  pthatWeighted_p->DrawCopy("HIST E1 P SAME");

  gPad->SetLogy();
  gStyle->SetOptStat(0);
  pthatCanv_p->SaveAs(("pdfDir/" + dateStr + "/pthatComp_" + extraTag + "_" + dateStr + ".pdf").c_str());

  delete pthatCanv_p;

  pthat_p->Write("", TObject::kOverwrite);
  pthatWeighted_p->Write("", TObject::kOverwrite);

  etaMiss_p->Write("", TObject::kOverwrite);
  phiMiss_p->Write("", TObject::kOverwrite);

  etaMiss_CMSSWJet_p->Write("", TObject::kOverwrite);
  phiMiss_CMSSWJet_p->Write("", TObject::kOverwrite);
  

  std::vector<std::vector<std::string> > algoCompInRhoBins;
  std::vector<std::vector<std::string> > algoCompInRhoBinsInv;
  std::vector<std::vector<std::string> > threshCompInRhoBins;
  std::vector<std::vector<std::string> > threshCompInRhoBinsInv;

  for(Int_t rI = 0; rI < nRhoBins; ++rI){
    algoCompInRhoBins.push_back({});
    algoCompInRhoBinsInv.push_back({});
    threshCompInRhoBins.push_back({});
    threshCompInRhoBinsInv.push_back({});
  }

  if(doForest){
    std::vector< std::vector<TH1*> > num, num2;
    std::vector<std::string> algo;
    std::vector<TH1*> denom, denom2;
    std::vector<std::string> rhoStr, rhoStr2;
    std::vector<int> rhoPos, rhoPos2;
    std::vector<std::vector<int> > subTypePos, subTypePos2;

    for(Int_t jI = 0; jI < nJetTrees; ++jI){
      for(Int_t cI = 0; cI < nRhoBins; ++cI){
	std::string centStr2 =  std::to_string(rhoBinsLow[nRhoBins-1]) + " < #rho < " + std::to_string(rhoBinsHi[0]);
	if(cI < nRhoBins ) centStr2 = std::to_string(rhoBinsLow[cI]) + " < #rho < " + std::to_string(rhoBinsHi[cI]);

	for(Int_t aI = 0; aI < nJtAbsEtaBins+1; ++aI){
	  leadingJetPt_h[jI][cI][aI]->Write("", TObject::kOverwrite);
	  	  
	  for(Int_t lI = 0; lI < nL1Thresh; ++lI){
	    std::vector<TH1*> temp1, temp2;
	    std::vector<int> type1, type2;
	    for(Int_t sI = 0; sI < nSubType; ++sI){
	      if(!doType[sI]) continue;
	      leadingJetPt_L1Pt_h[jI][cI][aI][sI][lI]->Write("", TObject::kOverwrite);
	      leadingJetPt_L1Pt_DR_h[jI][cI][aI][sI][lI]->Write("", TObject::kOverwrite);
	      
	      if(jI == 0 && aI == 0){
		triggerEta_L1Pt_h[cI][sI][lI]->Write("", TObject::kOverwrite);
	      }
	      
	      temp1.push_back(leadingJetPt_L1Pt_h[jI][cI][aI][sI][lI]);
	      temp2.push_back(leadingJetPt_L1Pt_DR_h[jI][cI][aI][sI][lI]);
	      type1.push_back(sI);
	      type2.push_back(sI);
	    }

	    denom.push_back(leadingJetPt_h[jI][cI][aI]);
            denom.push_back(leadingJetPt_h[jI][cI][aI]);
	    num.push_back(temp1);
	    num.push_back(temp2);
	    rhoStr.push_back(centStr2);
	    rhoStr.push_back(centStr2);
	    rhoPos.push_back(cI);
	    rhoPos.push_back(cI);
	    subTypePos.push_back(type1);
	    subTypePos.push_back(type2);

	    if(aI == 0) jetEta_L1Pt_h[cI][jI][lI]->Write("", TObject::kOverwrite);
	  }


	  for(Int_t sI = 0; sI < nSubType; ++sI){
	    if(!doType[sI]) continue;

	    std::vector<TH1*> temp1, temp2;

	    for(Int_t lI = 0; lI < nL1Thresh; ++lI){
	      temp1.push_back(leadingJetPt_L1Pt_h[jI][cI][aI][sI][lI]);
              temp2.push_back(leadingJetPt_L1Pt_DR_h[jI][cI][aI][sI][lI]);
	    }

	    denom2.push_back(leadingJetPt_h[jI][cI][aI]);
	    denom2.push_back(leadingJetPt_h[jI][cI][aI]);
	    num2.push_back(temp1);
	    num2.push_back(temp2);
	    algo.push_back(subType[sI]);
	    algo.push_back(subType[sI]);
	    rhoStr2.push_back(centStr2);
	    rhoStr2.push_back(centStr2);
	    rhoPos2.push_back(cI);
	    rhoPos2.push_back(cI);
	  }
	}
      }
    }

    for(Int_t jI = 0; jI < nJetTrees; ++jI){
      for(Int_t cI = 0; cI < nRhoBins+1; ++cI){
	std::string centStr = "AllRho";
	if(cI != nRhoBins) centStr = "Rho" + std::to_string(rhoBinsLow[cI]) + "to" + std::to_string(rhoBinsHi[cI]);

	TCanvas* meanOverlay_p = new TCanvas("meanOverlay_p", "", 450., 450.);
	meanOverlay_p->SetTopMargin(0.01);
	meanOverlay_p->SetRightMargin(0.01);
	meanOverlay_p->SetBottomMargin(0.08);
	meanOverlay_p->SetLeftMargin(0.08);

	TCanvas* sigmaOverlay_p = new TCanvas("sigmaOverlay_p", "", 450., 450.);
	sigmaOverlay_p->SetTopMargin(0.01);
	sigmaOverlay_p->SetRightMargin(0.01);
	sigmaOverlay_p->SetBottomMargin(0.08);
	sigmaOverlay_p->SetLeftMargin(0.08);

	Double_t maxMean = -1;
	Double_t minMean = 100000000;

	Double_t maxSigma = -1;
	Double_t minSigma = 100000000;

	for(Int_t sI = 0; sI < nSubType; ++sI){	  
	  if(!doType[sI]) continue;

	  Int_t x,y;
	  getXYBins(nL1ResPtBins, &x, &y);

	  TCanvas* canv_Res_p = new TCanvas("canv_Res_p", "", x*300., y*300.);
	  canv_Res_p->Divide(x, y);

	  for(Int_t rI = 0; rI < nL1ResPtBins; ++rI){
	    jetL1PtOverOffline_Mean_h[jI][cI][sI]->SetBinContent(rI+1, jetL1PtOverOffline_h[jI][cI][sI][rI]->GetMean());
	    jetL1PtOverOffline_Sigma_h[jI][cI][sI]->SetBinContent(rI+1, jetL1PtOverOffline_h[jI][cI][sI][rI]->GetStdDev());

	    jetL1PtOverOffline_Mean_h[jI][cI][sI]->SetBinError(rI+1, jetL1PtOverOffline_h[jI][cI][sI][rI]->GetMeanError());
	    jetL1PtOverOffline_Sigma_h[jI][cI][sI]->SetBinError(rI+1, jetL1PtOverOffline_h[jI][cI][sI][rI]->GetStdDevError());
	    
	    canv_Res_p->cd();
	    canv_Res_p->cd(rI+1);

	    jetL1PtOverOffline_h[jI][cI][sI][rI]->SetTitle((prettyString(l1ResPtBins[rI], false, 2) + " < p_{T,Offline} < " + prettyString(l1ResPtBins[rI+1], false, 2)).c_str());
	    jetL1PtOverOffline_h[jI][cI][sI][rI]->DrawCopy("HIST E1");
	    gStyle->SetOptStat(0);
	    
	    jetL1PtOverOffline_h[jI][cI][sI][rI]->Write("", TObject::kOverwrite);
	    delete jetL1PtOverOffline_h[jI][cI][sI][rI];
	  }
      
	  canv_Res_p->SaveAs(("pdfDir/" + dateStr + "/jetL1PtOverOffline_" + centStr + "_" + subType[sI] + "_" + jetNames.at(jI) + "_Points_" + extraTag + "_" + dateStr + ".pdf").c_str());
	  delete canv_Res_p;
	  
	  jetL1PtOverOffline_Mean_h[jI][cI][sI]->Write("", TObject::kOverwrite);
	  if(jetL1PtOverOffline_Mean_h[jI][cI][sI]->GetMaximum() > maxMean) maxMean = jetL1PtOverOffline_Mean_h[jI][cI][sI]->GetMaximum();
	  if(jetL1PtOverOffline_Mean_h[jI][cI][sI]->GetMinimum() < minMean) minMean = jetL1PtOverOffline_Mean_h[jI][cI][sI]->GetMinimum();


	  jetL1PtOverOffline_Sigma_h[jI][cI][sI]->Write("", TObject::kOverwrite);
	  if(jetL1PtOverOffline_Sigma_h[jI][cI][sI]->GetMaximum() > maxSigma) maxSigma = jetL1PtOverOffline_Sigma_h[jI][cI][sI]->GetMaximum();
	  if(jetL1PtOverOffline_Sigma_h[jI][cI][sI]->GetMinimum() < minSigma) minSigma = jetL1PtOverOffline_Sigma_h[jI][cI][sI]->GetMinimum();
	}

	double interval = maxMean - minMean;
	maxMean += interval/10.;
	minMean -= interval/10.;
	if(minMean < 0) minMean = 0;

	interval = maxSigma - minSigma;
	maxSigma += interval/10.;
	minSigma -= interval/10.;
	if(minSigma < 0) minSigma = 0;

	bool isDrawn = false;
	
	for(Int_t sI = 0; sI < nSubType; ++sI){	  
	  if(!doType[sI]) continue;

	  meanOverlay_p->cd();
	  jetL1PtOverOffline_Mean_h[jI][cI][sI]->SetMarkerStyle(styles[sI]);
	  jetL1PtOverOffline_Mean_h[jI][cI][sI]->SetMarkerSize(1.2);
	  jetL1PtOverOffline_Mean_h[jI][cI][sI]->SetMarkerColor(vg.getColor(sI));
	  jetL1PtOverOffline_Mean_h[jI][cI][sI]->SetLineColor(vg.getColor(sI));
	  jetL1PtOverOffline_Mean_h[jI][cI][sI]->SetMaximum(maxMean);
	  jetL1PtOverOffline_Mean_h[jI][cI][sI]->SetMinimum(minMean);

	  if(!isDrawn) jetL1PtOverOffline_Mean_h[jI][cI][sI]->DrawCopy("HIST E1 P");
	  else jetL1PtOverOffline_Mean_h[jI][cI][sI]->DrawCopy("HIST E1 P SAME");
	  
	  sigmaOverlay_p->cd();
	  jetL1PtOverOffline_Sigma_h[jI][cI][sI]->SetMarkerStyle(styles[sI]);
	  jetL1PtOverOffline_Sigma_h[jI][cI][sI]->SetMarkerSize(1.2);
	  jetL1PtOverOffline_Sigma_h[jI][cI][sI]->SetMarkerColor(vg.getColor(sI));
	  jetL1PtOverOffline_Sigma_h[jI][cI][sI]->SetLineColor(vg.getColor(sI));
	  jetL1PtOverOffline_Sigma_h[jI][cI][sI]->SetMaximum(maxSigma);
	  jetL1PtOverOffline_Sigma_h[jI][cI][sI]->SetMinimum(minSigma);

	  if(!isDrawn) jetL1PtOverOffline_Sigma_h[jI][cI][sI]->DrawCopy("HIST E1 P");
	  else jetL1PtOverOffline_Sigma_h[jI][cI][sI]->DrawCopy("HIST E1 P SAME");

	  isDrawn = true;

	  delete jetL1PtOverOffline_Mean_h[jI][cI][sI];
	  delete jetL1PtOverOffline_Sigma_h[jI][cI][sI];
	}

	meanOverlay_p->SaveAs(("pdfDir/" + dateStr + "/jetL1PtOverOffline_" + centStr + "_" + jetNames.at(jI) + "_Mean_" + extraTag + "_" + dateStr + ".pdf").c_str());
	delete meanOverlay_p;

	sigmaOverlay_p->SaveAs(("pdfDir/" + dateStr + "/jetL1PtOverOffline_" + centStr + "_" + jetNames.at(jI) + "_Sigma_" + extraTag + "_" + dateStr + ".pdf").c_str());
	delete sigmaOverlay_p;


      }
    }

  

  
    //    std::cout << "How many hist: " << num.size() << std::endl;
    for(unsigned int i = 0; i < num.size(); ++i){
      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);
      canv_p->SetLeftMargin(0.12);
      canv_p->SetBottomMargin(0.12);
     
      std::vector<TH1*> temp = num.at(i);
      std::vector<int> type = subTypePos.at(i);
      //      std::cout << " Sub hist: " << temp.size() << std::endl;


      const Int_t nBins = temp.at(0)->GetNbinsX();
      Double_t bins[nBins+1];
      
      for(Int_t bIX = 0; bIX < nBins+1; ++bIX){
	bins[bIX] = temp.at(0)->GetBinLowEdge(bIX+1);
      }

      std::string tempX = temp.at(0)->GetXaxis()->GetTitle();
      TH1F* tempHist_p = new TH1F("tempHist_h", (";" + tempX + ";Efficiency").c_str(), nBins, bins);
      tempHist_p->SetMinimum(0.0);
      tempHist_p->SetMaximum(1.2);     
      
      canv_p->cd();
      
      centerTitles(tempHist_p);
      tempHist_p->DrawCopy();

      const Int_t nAPt = temp.size();
      TGraphAsymmErrors* aPt_p[nAPt];
    
      for(unsigned int sI = 0; sI < temp.size(); ++sI){
	//	temp.at(sI)->Print("ALL");
	//	denom.at(i)->Print("ALL");

	aPt_p[sI] = new TGraphAsymmErrors();
	//	if(!doType[sI]) continue;
	aPt_p[sI]->BayesDivide(temp.at(sI), denom.at(i));
	aPt_p[sI]->SetMarkerColor(vg.getColor(type.at(sI)));
	aPt_p[sI]->SetLineColor(vg.getColor(type.at(sI)));
	aPt_p[sI]->SetMarkerStyle(styles[type.at(sI)]);
	aPt_p[sI]->SetMarkerSize(1.2);

	aPt_p[sI]->Draw("P");	
      }

      gStyle->SetOptStat(0);

      dummyLeg_p->Draw("SAME");
      line_p->DrawLine(bins[0], 1., bins[nBins], 1.);
    
      std::string histName = temp.at(0)->GetName();
      histName.replace(histName.find("leadingJetPt_"), std::string("leadingJetPt_").size(), "");
      histName.replace(histName.find("EventTree_"), std::string("EventTree_").size(), "");
      globalAlgoString = "";

      std::string thresh = histName.substr(histName.find("L1Pt")+4, histName.size());
      thresh.replace(thresh.find("_"), thresh.size(), "");
      thresh.replace(thresh.find("p"), 1, ".");
      thresh = "L1 p_{T} > " + thresh;
      label_p->DrawLatex(0.15, 0.94, thresh.c_str());
      label_p->DrawLatex(0.35, 0.94, rhoStr.at(i).c_str());
      if(histName.find("DR_") != std::string::npos) label_p->DrawLatex(0.65, 0.94, "#DeltaR Matched");
      else label_p->DrawLatex(0.65, 0.94, "No #DeltaR Match");

      std::string absEtaStr = histName.substr(histName.find("AbsEta"), histName.size());
      absEtaStr.replace(absEtaStr.find("_"), absEtaStr.size(), "");
      label_p->DrawLatex(0.15, 0.88, absEtaStr.c_str());


      std::string saveName = "turnOn_" + histName + globalAlgoString + "_" + extraTag + "_" + dateStr + ".pdf";
      algoCompInRhoBins.at(rhoPos.at(i)).push_back(saveName);
      
      canv_p->SaveAs(("pdfDir/" + dateStr + "/" + saveName).c_str());
      delete canv_p;

      for(unsigned int sI = 0; sI < temp.size(); ++sI){
	delete aPt_p[sI];
      }

      delete tempHist_p;
    }

    for(unsigned int i = 0; i < num.size(); ++i){
      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);
      canv_p->SetLeftMargin(0.12);
      canv_p->SetBottomMargin(0.12);
     
      std::vector<TH1*> temp = num.at(i);
      std::vector<int> type = subTypePos.at(i);

      const Int_t nBins = temp.at(0)->GetNbinsX();
      Double_t bins[nBins+1];
      
      for(Int_t bIX = 0; bIX < nBins+1; ++bIX){
	bins[bIX] = temp.at(0)->GetBinLowEdge(bIX+1);
      }

      std::string tempX = temp.at(0)->GetXaxis()->GetTitle();
      TH1F* tempHist_p = new TH1F("tempHist_h", (";" + tempX + ";Inefficiency").c_str(), nBins, bins);
      tempHist_p->SetMinimum(0.0);
      tempHist_p->SetMaximum(1.2);
            
      canv_p->cd();

      const Int_t nAPt = temp.size();
      TGraphAsymmErrors* aPt_p[nAPt];

      Double_t minVal = 100;
      Double_t maxVal = -1;
      
      for(unsigned int sI = 0; sI < temp.size(); ++sI){
	aPt_p[sI] = new TGraphAsymmErrors();
	//	if(!doType[sI]) continue;

	TH1F* tempClone = (TH1F*)temp.at(sI)->Clone("tempClone");
	for(Int_t bI = 0; bI < tempClone->GetNbinsX(); ++bI){
	  Double_t val = denom.at(i)->GetBinContent(bI+1) - tempClone->GetBinContent(bI+1);
	  Double_t err = TMath::Sqrt(val);
	  tempClone->SetBinContent(bI+1, val);
	  tempClone->SetBinError(bI+1, err);
	}

	aPt_p[sI]->BayesDivide(tempClone, denom.at(i));
	aPt_p[sI]->SetMarkerColor(vg.getColor(type.at(sI)));
	aPt_p[sI]->SetLineColor(vg.getColor(type.at(sI)));
	aPt_p[sI]->SetMarkerStyle(styles[type.at(sI)]);
	aPt_p[sI]->SetMarkerSize(1.2);

	Double_t xVal, yVal;
	for(Int_t pI = 0; pI < aPt_p[sI]->GetN(); ++pI){
	  aPt_p[sI]->GetPoint(pI, xVal, yVal);
 	  if(yVal > maxVal) maxVal = yVal;
 	  if(yVal < minVal && yVal > 0) minVal = yVal;
	  if(yVal <= 0) aPt_p[sI]->SetPointError(pI, 0, 0, 0, 0);
	}

	delete tempClone;
      }      

      tempHist_p->SetMaximum(4.);
      tempHist_p->SetMinimum(.000001);

      centerTitles(tempHist_p);
      tempHist_p->DrawCopy();

      for(unsigned int sI = 0; sI < temp.size(); ++sI){
	aPt_p[sI]->Draw("P");
      }


      gStyle->SetOptStat(0);

      dummyLeg_p->Draw("SAME");
      line_p->DrawLine(bins[0], 1., bins[nBins], 1.);
      gPad->SetLogy();
      std::string histName = temp.at(0)->GetName();
      histName.replace(histName.find("leadingJetPt_"), std::string("leadingJetPt_").size(), "");
      histName.replace(histName.find("EventTree_"), std::string("EventTree_").size(), "");
      globalAlgoString = "";
      
      std::string thresh = histName.substr(histName.find("L1Pt")+4, histName.size());
      thresh.replace(thresh.find("_"), thresh.size(), "");
      thresh.replace(thresh.find("p"), 1, ".");
      thresh = "L1 p_{T} > " + thresh;
      label_p->DrawLatex(0.15, 0.94, thresh.c_str());
      label_p->DrawLatex(0.35, 0.94, rhoStr.at(i).c_str());
      if(histName.find("DR_") != std::string::npos) label_p->DrawLatex(0.65, 0.94, "#DeltaR Matched");
      else label_p->DrawLatex(0.65, 0.94, "No #DeltaR Match");

      std::string absEtaStr = histName.substr(histName.find("AbsEta"), histName.size());
      absEtaStr.replace(absEtaStr.find("_"), absEtaStr.size(), "");
      label_p->DrawLatex(0.15, 0.88, absEtaStr.c_str());

      std::string saveName = "invTurnOn_" + histName +  globalAlgoString + "_" + extraTag + "_" + dateStr + ".pdf";
      algoCompInRhoBinsInv.at(rhoPos.at(i)).push_back(saveName);

      canv_p->SaveAs(("pdfDir/" + dateStr + "/" + saveName).c_str());
      delete canv_p;

      for(unsigned int sI = 0; sI < temp.size(); ++sI){
	delete aPt_p[sI];
      }

      delete tempHist_p;
    }
 

    // TRIG TURN ON
  
    for(unsigned int i = 0; i < num2.size(); ++i){
      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);
      canv_p->SetLeftMargin(0.12);
      canv_p->SetBottomMargin(0.12);
     
      std::vector<TH1*> temp = num2.at(i);

      const Int_t nBins = temp.at(0)->GetNbinsX();
      Double_t bins[nBins+1];
      
      for(Int_t bIX = 0; bIX < nBins+1; ++bIX){
	bins[bIX] = temp.at(0)->GetBinLowEdge(bIX+1);
      }

      std::string tempX = temp.at(0)->GetXaxis()->GetTitle();
      TH1F* tempHist_p = new TH1F("tempHist_h", (";" + tempX + ";Efficiency").c_str(), nBins, bins);
      tempHist_p->SetMinimum(0.0);
      tempHist_p->SetMaximum(1.2);
            
      canv_p->cd();
      
      centerTitles(tempHist_p);
      tempHist_p->DrawCopy();

      const Int_t nAPt = temp.size();
      TGraphAsymmErrors* aPt_p[nAPt];

      std::cout << "CHECKING" << std::endl;
      for(unsigned int sI = 0; sI < temp.size(); ++sI){
	std::cout << " " << temp.at(sI)->GetName() << std::endl;
	aPt_p[sI] = new TGraphAsymmErrors();
	aPt_p[sI]->BayesDivide(temp.at(sI), denom2.at(i));
	aPt_p[sI]->SetMarkerColor(vg.getColor(sI));
	aPt_p[sI]->SetLineColor(vg.getColor(sI));
	aPt_p[sI]->SetMarkerStyle(styles[sI%nL1Thresh]);
	aPt_p[sI]->SetMarkerSize(1.2);

	aPt_p[sI]->Draw("P");	
      }

      gStyle->SetOptStat(0);

      //      dummyLeg_p->Draw("SAME");
      line_p->DrawLine(bins[0], 1., bins[nBins], 1.);
    
      std::string histName = temp.at(0)->GetName();
      histName.replace(histName.find("leadingJetPt_"), std::string("leadingJetPt_").size(), "");
      histName.replace(histName.find("EventTree_"), std::string("EventTree_").size(), "");
      globalAlgoString = "_Algo" + algo.at(i);
     
     label_p->DrawLatex(0.65, 0.88, algo.at(i).c_str());
      label_p->DrawLatex(0.15, 0.94, rhoStr2.at(i).c_str());
      if(histName.find("DR_") != std::string::npos) label_p->DrawLatex(0.65, 0.94, "#DeltaR Matched");
      else label_p->DrawLatex(0.65, 0.94, "No #DeltaR Match");

      std::string absEtaStr = histName.substr(histName.find("AbsEta"), histName.size());
      absEtaStr.replace(absEtaStr.find("_"), absEtaStr.size(), "");
      label_p->DrawLatex(0.15, 0.88, absEtaStr.c_str());

      dummyLeg4_p->Draw("SAME");
      std::string saveName = "turnOn_" + histName + globalAlgoString + "_" + extraTag + "_" + dateStr + ".pdf";
      threshCompInRhoBins.at(rhoPos2.at(i)).push_back(saveName);
      canv_p->SaveAs(("pdfDir/" + dateStr + "/" + saveName).c_str());
      delete canv_p;

      for(unsigned int sI = 0; sI < temp.size(); ++sI){
	delete aPt_p[sI];
      }

      delete tempHist_p;
    }

    for(unsigned int i = 0; i < num2.size(); ++i){
      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);
      canv_p->SetLeftMargin(0.12);
      canv_p->SetBottomMargin(0.12);
     
      std::vector<TH1*> temp = num2.at(i);

      const Int_t nBins = temp.at(0)->GetNbinsX();
      Double_t bins[nBins+1];
      
      for(Int_t bIX = 0; bIX < nBins+1; ++bIX){
	bins[bIX] = temp.at(0)->GetBinLowEdge(bIX+1);
      }

      std::string tempX = temp.at(0)->GetXaxis()->GetTitle();
      TH1F* tempHist_p = new TH1F("tempHist_h", (";" + tempX + ";Inefficiency").c_str(), nBins, bins);
      tempHist_p->SetMinimum(0.0);
      tempHist_p->SetMaximum(1.2);
            
      canv_p->cd();

      const Int_t nAPt = temp.size();
      TGraphAsymmErrors* aPt_p[nAPt];

      Double_t minVal = 100;
      Double_t maxVal = -1;
      
      for(unsigned int sI = 0; sI < temp.size(); ++sI){
	aPt_p[sI] = new TGraphAsymmErrors();

	TH1F* tempClone = (TH1F*)temp.at(sI)->Clone("tempClone");
	for(Int_t bI = 0; bI < tempClone->GetNbinsX(); ++bI){
	  Double_t val = denom2.at(i)->GetBinContent(bI+1) - tempClone->GetBinContent(bI+1);
	  Double_t err = TMath::Sqrt(val);
	  tempClone->SetBinContent(bI+1, val);
	  tempClone->SetBinError(bI+1, err);
	}

	aPt_p[sI]->BayesDivide(tempClone, denom2.at(i));
	aPt_p[sI]->SetMarkerColor(vg.getColor(sI));
	aPt_p[sI]->SetLineColor(vg.getColor(sI));
	aPt_p[sI]->SetMarkerStyle(styles[sI%nL1Thresh]);
	aPt_p[sI]->SetMarkerSize(1.2);

	Double_t xVal, yVal;
	for(Int_t pI = 0; pI < aPt_p[sI]->GetN(); ++pI){
	  aPt_p[sI]->GetPoint(pI, xVal, yVal);
 	  if(yVal > maxVal) maxVal = yVal;
 	  if(yVal < minVal && yVal > 0) minVal = yVal;
	  if(yVal <= 0) aPt_p[sI]->SetPointError(pI, 0, 0, 0, 0);
	}

	delete tempClone;
      }      

      tempHist_p->SetMaximum(4.);
      tempHist_p->SetMinimum(.000001);

      centerTitles(tempHist_p);
      tempHist_p->DrawCopy();

      for(unsigned int sI = 0; sI < temp.size(); ++sI){
	aPt_p[sI]->Draw("P");
      }


      gStyle->SetOptStat(0);

      //      dummyLeg_p->Draw("SAME");
      line_p->DrawLine(bins[0], 1., bins[nBins], 1.);
      gPad->SetLogy();
      std::string histName = temp.at(0)->GetName();
      histName.replace(histName.find("leadingJetPt_"), std::string("leadingJetPt_").size(), "");
      histName.replace(histName.find("EventTree_"), std::string("EventTree_").size(), "");
      globalAlgoString = "_Algo" + algo.at(i);
    
      label_p->DrawLatex(0.65, 0.88, algo.at(i).c_str());
      label_p->DrawLatex(0.15, 0.94, rhoStr2.at(i).c_str());
      if(histName.find("DR_") != std::string::npos) label_p->DrawLatex(0.65, 0.94, "#DeltaR Matched");
      else label_p->DrawLatex(0.65, 0.94, "No #DeltaR Match");

      std::string absEtaStr = histName.substr(histName.find("AbsEta"), histName.size());
      absEtaStr.replace(absEtaStr.find("_"), absEtaStr.size(), "");
      label_p->DrawLatex(0.15, 0.88, absEtaStr.c_str());

      /*
      std::string thresh = histName.substr(histName.find("L1Pt")+4, histName.size());
      thresh.replace(thresh.find("_"), thresh.size(), "");
      thresh.replace(thresh.find("p"), 1, ".");
      thresh = "L1 p_{T} > " + thresh;
      label_p->DrawLatex(0.5, 0.94, thresh.c_str());
      */

      std::string saveName = "invTurnOn_" + histName + globalAlgoString + "_" + extraTag + "_" + dateStr + ".pdf";


      threshCompInRhoBinsInv.at(rhoPos2.at(i)).push_back(saveName);
      canv_p->SaveAs(("pdfDir/" + dateStr + "/" + saveName).c_str());
      delete canv_p;

      for(unsigned int sI = 0; sI < temp.size(); ++sI){
	delete aPt_p[sI];
      }

      delete tempHist_p;
    }
  }

  // END TRIG TURN ON

  dummyLeg_p->SetX1NDC(.15);
  dummyLeg_p->SetX2NDC(.4);
  dummyLeg_p->SetY1NDC(.65);
  dummyLeg_p->SetY2NDC(.9);

  std::vector<std::vector<std::string> > l1TrigCanvName_Rho;
  
  for(Int_t cI = 0; cI < nRhoBins; ++cI){
    l1TrigCanvName_Rho.push_back({});
  }
  
  for(Int_t cI = 0; cI < nRhoBins+1; ++cI){
    std::string centStr = "AllRho";
    if(cI != nRhoBins) centStr = "Rho" + std::to_string(rhoBinsLow[cI]) + "to" + std::to_string(rhoBinsHi[cI]);

    for(Int_t lI = 0; lI < nL1Thresh; ++lI){
      int firstPos = -1;
      for(Int_t sI = 0; sI < nSubType; ++sI){
	if(doType[sI]){
	  firstPos = sI;
	  break;
	}
      }


      if(triggerEta_L1Pt_h[cI][firstPos][lI]->GetEntries() == 0) continue;

      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);
      canv_p->SetLeftMargin(0.08);
      canv_p->SetBottomMargin(0.08);

      Double_t min = 0;
      Double_t max = -1;
      for(Int_t sI = 0; sI < nSubType; ++sI){
	centerTitles(triggerEta_L1Pt_h[cI][sI][lI]);
	triggerEta_L1Pt_h[cI][sI][lI]->SetMarkerColor(vg.getColor(sI));
	triggerEta_L1Pt_h[cI][sI][lI]->SetLineColor(vg.getColor(sI));
	triggerEta_L1Pt_h[cI][sI][lI]->SetMarkerStyle(styles[sI]);
	triggerEta_L1Pt_h[cI][sI][lI]->SetMarkerSize(1.2);
	
	if(triggerEta_L1Pt_h[cI][sI][lI]->GetMaximum() > max) max = triggerEta_L1Pt_h[cI][sI][lI]->GetMaximum();
      }
      
      max *= 1.3;
      
      for(Int_t sI = 0; sI < nSubType; ++sI){
	if(!doType[sI]) continue;

	triggerEta_L1Pt_h[cI][sI][lI]->SetMaximum(max);
	triggerEta_L1Pt_h[cI][sI][lI]->SetMinimum(min);
	
	if(sI == firstPos) triggerEta_L1Pt_h[cI][sI][lI]->DrawCopy("HIST E1 P");
	else triggerEta_L1Pt_h[cI][sI][lI]->DrawCopy("HIST E1 P SAME");
      }
      
      dummyLeg_p->Draw("SAME");

      
      const std::string partialSaveName = "triggerEta_L1Pt" + prettyString(l1ThreshPt[lI], 1, true) + "_" + centStr + "_" + globalAlgoString + "_" + extraTag + "_" + dateStr + ".pdf";
      const std::string saveName = "pdfDir/" + dateStr + "/" + partialSaveName;
      canv_p->SaveAs(saveName.c_str());
      delete canv_p;
    
      if(cI != nRhoBins) l1TrigCanvName_Rho.at(cI).push_back(partialSaveName);
    }
  }

  for(Int_t cI = 0; cI < nRhoBins+1; ++cI){
    std::string centStr = "AllRho";
    if(cI != nRhoBins) centStr = "Rho" + std::to_string(rhoBinsLow[cI]) + "to" + std::to_string(rhoBinsHi[cI]);

    for(Int_t lI = 0; lI < nL1Thresh; ++lI){
      if(jetEta_L1Pt_h[cI][0][lI]->GetEntries() == 0) continue;

      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);
      canv_p->SetLeftMargin(0.08);
      canv_p->SetBottomMargin(0.08);
      
      Double_t min = 0;
      Double_t max = -1;

      for(Int_t jI = 0; jI < nJetTrees; ++jI){
	centerTitles(jetEta_L1Pt_h[cI][jI][lI]);
	jetEta_L1Pt_h[cI][jI][lI]->SetMarkerColor(vg.getColor(jI));
	jetEta_L1Pt_h[cI][jI][lI]->SetLineColor(vg.getColor(jI));
	jetEta_L1Pt_h[cI][jI][lI]->SetMarkerStyle(styles[jI]);
	jetEta_L1Pt_h[cI][jI][lI]->SetMarkerSize(1.2);
	
	if(jetEta_L1Pt_h[cI][jI][lI]->GetMaximum() > max) max = jetEta_L1Pt_h[cI][jI][lI]->GetMaximum();
      }
      
      max *= 1.3;
      
      for(Int_t jI = 0; jI < nJetTrees; ++jI){
	jetEta_L1Pt_h[cI][jI][lI]->SetMaximum(max);
	jetEta_L1Pt_h[cI][jI][lI]->SetMinimum(min);
	
	if(jI == 0) jetEta_L1Pt_h[cI][jI][lI]->DrawCopy("HIST E1 P");
	else jetEta_L1Pt_h[cI][jI][lI]->DrawCopy("HIST E1 P SAME");
      }
      
      dummyLeg2_p->Draw("SAME");
      
      const std::string saveName = "pdfDir/" + dateStr + "/jetEta_L1Pt" + prettyString(l1ThreshPt[lI], 1, true) + "_" + centStr + "_" + globalAlgoString + "_" + extraTag + "_" + dateStr + ".pdf";
      canv_p->SaveAs(saveName.c_str());
      delete canv_p;
    }
  }

  for(unsigned int sI = 0; sI < nSubType; ++sI){
    delete dummys_p[sI];
  }

  for(Int_t jI = 0; jI < nJetTrees; ++jI){
    delete dummys2_p[jI];
  }
  
  delete dummyLeg_p;
  delete dummyLeg2_p;
  delete line_p;
  delete label_p;
  
  delete pthat_p;
  delete pthatWeighted_p;

  delete etaMiss_p;
  delete phiMiss_p;

  delete etaMiss_CMSSWJet_p;
  delete phiMiss_CMSSWJet_p;
  
  if(doForest){
    for(Int_t jI = 0; jI < nJetTrees; ++jI){
      for(Int_t cI = 0; cI < nRhoBins; ++cI){
	for(Int_t aI = 0; aI < nJtAbsEtaBins+1; ++aI){
	  delete leadingJetPt_h[jI][cI][aI];
	  
	  for(Int_t sI = 0; sI < nSubType; ++sI){
	    for(Int_t lI = 0; lI < nL1Thresh; ++lI){
	      delete leadingJetPt_L1Pt_h[jI][cI][aI][sI][lI];
	      delete leadingJetPt_L1Pt_DR_h[jI][cI][aI][sI][lI];
	      
	      if(jI == 0 && aI == 0) delete triggerEta_L1Pt_h[cI][sI][lI];
	    }
	  }
	}
      }
    }
  }
  
  outFile_p->Close();
  delete outFile_p;;

  std::string fileNameTex = "pdfDir/" + dateStr + "/l1OfflineSubtract_" + extraTag + "_" + dateStr + ".tex";
  std::ofstream texFile(fileNameTex.c_str());

  texFile << "\\RequirePackage{xspace}" << std::endl;
  texFile << "\\RequirePackage{amsmath}" << std::endl;
  texFile << std::endl;

  texFile << "\\documentclass[xcolor=dvipsnames]{beamer}" << std::endl;
  texFile << "\\usetheme{Warsaw}" << std::endl;
  texFile << "\\setbeamercolor{structure}{fg=NavyBlue!90!NavyBlue}" << std::endl;
  texFile << "\\setbeamercolor{footlinecolor}{fg=white,bg=lightgray}" << std::endl;
  texFile << std::endl;

  texFile << "\\newcommand{\\pt}{\\ensuremath{p_{\\mathrm{T}}}\\xspace}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamersize{text margin left=5pt,text margin right=5pt}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamertemplate{frametitle}" << std::endl;
  texFile << "{" << std::endl;
  texFile << "    \\nointerlineskip" << std::endl;
  texFile << "    \\begin{beamercolorbox}[sep=0.3cm, ht=1.8em, wd=\\paperwidth]{frametitle}" << std::endl;
  texFile << "        \\vbox{}\\vskip-2ex%" << std::endl;
  texFile << "        \\strut\\insertframetitle\\strut" << std::endl;
  texFile << "        \\vskip-0.8ex%" << std::endl;
  texFile << "    \\end{beamercolorbox}" << std::endl;
  texFile << "}" << std::endl;
  texFile << std::endl;
  texFile << "\\setbeamertemplate{footline}{%" << std::endl;
  texFile << "  \\begin{beamercolorbox}[sep=.8em,wd=\\paperwidth,leftskip=0.5cm,rightskip=0.5cm]{footlinecolor}" << std::\
    endl;
  texFile << "    \\hspace{0.3cm}%" << std::endl;
  texFile << "    \\hfill\\insertauthor \\hfill\\insertpagenumber" << std::endl;
  texFile << "  \\end{beamercolorbox}%" << std::endl;
  texFile << "}" << std::endl;
  texFile << "\\setbeamertemplate{navigation symbols}{}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamertemplate{itemize item}[circle]" << std::endl;
  texFile << "\\setbeamertemplate{itemize subitem}[circle]" << std::endl;
  texFile << "\\setbeamertemplate{itemize subsubitem}[circle]" << std::endl;
  texFile << "\\setbeamercolor{itemize item}{fg=black}" << std::endl;
  texFile << "\\setbeamercolor{itemize subitem}{fg=black}" << std::endl;
  texFile << "\\setbeamercolor{itemize subsubitem}{fg=black}" << std::endl;
  texFile << std::endl;

  texFile << "\\definecolor{links}{HTML}{00BFFF}" << std::endl;
  texFile << "\\hypersetup{colorlinks,linkcolor=,urlcolor=links}" << std::endl;
  texFile << std::endl;

  texFile << "\\author[CM]{Placeholder}" << std::endl;
  texFile << std::endl;

  texFile << "\\begin{document}" << std::endl;
  texFile << "\\begin{frame}" << std::endl;
  texFile << "\\frametitle{\\centerline{L1 Trigger (" << dateStr << ")}}" << std::endl;
  texFile << " \\begin{itemize}" << std::endl;
  texFile << "  \\fontsize{10}{10}\\selectfont" << std::endl;
  texFile << "  \\item{Placeholder}" << std::endl;
  texFile << "  \\begin{itemize}" << std::endl;
  texFile << "   \\fontsize{10}{10}\\selectfont" << std::endl;
  texFile << "   \\item{Placeholder}" << std::endl;
  texFile << "  \\end{itemize}" << std::endl;
  texFile << " \\end{itemize}" << std::endl;
  texFile << "\\end{frame}" << std::endl;

  texFile << std::endl;

  for(unsigned int aI = 0; aI < algoCompInRhoBins.at(0).size(); ++aI){
    std::string algo = algoCompInRhoBins.at(0).at(aI);
    std::string ptThresh = algo.substr(algo.find("L1Pt")+4, algo.size());
    ptThresh.replace(ptThresh.find("_"), ptThresh.size(), "");
    ptThresh.replace(ptThresh.find("p"), 1, ".");
    ptThresh = "$p_{T,L1} > $" + ptThresh;

    std::string drStr = "Not $\\Delta$R Matched";
    if(algo.find("DR_") != std::string::npos) drStr = "$\\Delta$R Matched";

    std::string offlineJet = "";
    if(algo.find("Calo") != std::string::npos) offlineJet = "Offline Calo. Jets";
    else if(algo.find("PF") != std::string::npos) offlineJet = "Offline PF Jets";
    else if(algo.find("Gen") != std::string::npos) offlineJet = "Offline Gen. Jets";

    std::string absEtaStr = "";
    if(algo.find("AbsEta0p0to3p0") != std::string::npos) absEtaStr = "0.0 $< |\\eta| <$ 3.0";
    else if(algo.find("AbsEta3p0to5p0") != std::string::npos) absEtaStr = "3.0 $< |\\eta| <$ 5.0";
    else absEtaStr = "0.0 $< |\\eta| <$ 5.0";

    texFile << "\\begin{frame}" << std::endl;
    texFile << "\\frametitle{\\centerline{Algo. comp., " << ptThresh << ", " << offlineJet << "}}" << std::endl;
    texFile << "\\begin{center}" << std::endl;
    for(unsigned int lI = 0; lI < algoCompInRhoBins.size(); ++lI){
      texFile << "\\includegraphics[width=" << 0.33 << "\\textwidth]{" << algoCompInRhoBins.at(algoCompInRhoBins.size() - 1 - lI).at(aI) << "}" << std::endl;
    }
    texFile << "\\end{center}" << std::endl;
    texFile << "\\begin{itemize}" << std::endl;
    texFile << "\\fontsize{8}{8}\\selectfont" << std::endl;
    texFile << "\\item{Left to Right: Increasing centrality.}" << std::endl;
    texFile << "\\item{" << ptThresh << "}" << std::endl;
    texFile << "\\item{" << drStr << "}" << std::endl;
    texFile << "\\item{" << absEtaStr << "}" << std::endl;
    texFile << "\\item{" << offlineJet << "}" << std::endl;
    texFile << "\\end{itemize}" << std::endl;
    texFile << "\\end{frame}" << std::endl;
  }

  std::cout << "Plotting " << l1TrigCanvName_Rho.at(0).size() << std::endl;
  for(unsigned int aI = 0; aI < l1TrigCanvName_Rho.at(0).size(); ++aI){
    std::string algo = l1TrigCanvName_Rho.at(0).at(aI);
    std::string ptThresh = algo.substr(algo.find("L1Pt")+4, algo.size());
    ptThresh.replace(ptThresh.find("_"), ptThresh.size(), "");
    ptThresh.replace(ptThresh.find("p"), 1, ".");
    ptThresh = "$p_{T,L1} > $" + ptThresh;

    texFile << "\\begin{frame}" << std::endl;
    texFile << "\\frametitle{\\centerline{Algo. $\eta$ comp., " << ptThresh << "}}" << std::endl;
    texFile << "\\begin{center}" << std::endl;
    for(unsigned int lI = 0; lI < l1TrigCanvName_Rho.size(); ++lI){
      texFile << "\\includegraphics[width=" << 0.33 << "\\textwidth]{" << l1TrigCanvName_Rho.at(l1TrigCanvName_Rho.size() - 1 - lI).at(aI) << "}" << std::endl;
    }
    texFile << "\\end{center}" << std::endl;
    texFile << "\\begin{itemize}" << std::endl;
    texFile << "\\fontsize{8}{8}\\selectfont" << std::endl;
    texFile << "\\item{Left to Right: Increasing centrality.}" << std::endl;
    texFile << "\\item{" << ptThresh << "}" << std::endl;
    texFile << "\\end{itemize}" << std::endl;
    texFile << "\\end{frame}" << std::endl;
  }


  for(unsigned int aI = 0; aI < threshCompInRhoBins.at(0).size(); ++aI){
    std::string algo = threshCompInRhoBins.at(0).at(aI);

    std::string subAlgo = algo;
    for(Int_t i = 0; i < nSubType; ++i){
      if(subAlgo.find(std::string("_" + subType[i] + "_")) != std::string::npos){
	subAlgo = subType[i];
	break;
      }
    }

    std::string drStr = "Not $\\Delta$R Matched";
    if(algo.find("DR_") != std::string::npos) drStr = "$\\Delta$R Matched";

    std::string offlineJet = "";
    if(algo.find("Calo") != std::string::npos) offlineJet = "Offline Calo. Jets";
    else if(algo.find("PF") != std::string::npos) offlineJet = "Offline PF Jets";
    else if(algo.find("Gen") != std::string::npos) offlineJet = "Offline Gen. Jets";

    std::string absEtaStr = "";
    if(algo.find("AbsEta0p0to3p0") != std::string::npos) absEtaStr = "0.0 $< |\\eta| <$ 3.0";
    else if(algo.find("AbsEta3p0to5p0") != std::string::npos) absEtaStr = "3.0 $< |\\eta| <$ 5.0";
    else absEtaStr = "0.0 $< |\\eta| <$ 5.0";

    texFile << "\\begin{frame}" << std::endl;
    texFile << "\\frametitle{\\centerline{Threshold comp., " << subAlgo << ", " << offlineJet << "}}" << std::endl;
    texFile << "\\begin{center}" << std::endl;
    for(unsigned int lI = 0; lI < threshCompInRhoBins.size(); ++lI){
      texFile << "\\includegraphics[width=" << 0.33 << "\\textwidth]{" << threshCompInRhoBins.at(threshCompInRhoBins.size() - 1 - lI).at(aI) << "}" << std::endl;
    }
    texFile << "\\end{center}" << std::endl;
    texFile << "\\begin{itemize}" << std::endl;
    texFile << "\\fontsize{8}{8}\\selectfont" << std::endl;
    texFile << "\\item{Left to Right: Increasing centrality.}" << std::endl;
    texFile << "\\item{" << subAlgo << "}" << std::endl;
    texFile << "\\item{" << drStr << "}" << std::endl;
    texFile << "\\item{" << absEtaStr << "}" << std::endl;
    texFile << "\\item{" << offlineJet << "}" << std::endl;
    texFile << "\\end{itemize}" << std::endl;
    texFile << "\\end{frame}" << std::endl;
  }

  texFile << "\\end{document}" << std::endl;

  texFile.close();

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc < 6 || argc > 13){
    std::cout << "Usage ./bin/l1OfflineSubtract.exe <inFileName> <forestName> <doWeights> <extraTag> <subTypeStr-upTo5>" << std::endl;
    return 1;
  }

  std::vector<std::string> subTypeStr;
  if(argc >= 6) subTypeStr.push_back(argv[5]);
  if(argc >= 7) subTypeStr.push_back(argv[6]);
  if(argc >= 8) subTypeStr.push_back(argv[7]);
  if(argc >= 9) subTypeStr.push_back(argv[8]);
  if(argc >= 10) subTypeStr.push_back(argv[9]);
  if(argc >= 11) subTypeStr.push_back(argv[10]);
  if(argc >= 12) subTypeStr.push_back(argv[11]);
  if(argc >= 13) subTypeStr.push_back(argv[12]);

  //std::cout << __LINE__ << std::endl;

  std::cout << "Arguments " << std::endl;
  for(Int_t aI = 0; aI < argc; ++aI){
    std::cout << " " << aI << "/" << argc << ": " << argv[aI] << std::endl;
  }

  std::vector<std::string> vect1, vect2;
  std::string argv1 = argv[1];
  std::string argv2 = argv[2];
  while(argv1.find(",") != std::string::npos){
    vect1.push_back(argv1.substr(0, argv1.find(",")));
    argv1.replace(0, argv1.find(",")+1,"");
  }
  if(argv1.size() != 0) vect1.push_back(argv1);

  //std::cout << __LINE__ << std::endl;

  while(argv2.find(",") != std::string::npos){
    vect2.push_back(argv2.substr(0, argv2.find(",")));
    argv2.replace(0, argv2.find(",")+1,"");
  }
  if(argv2.size() != 0) vect2.push_back(argv2);

  //std::cout << __LINE__ << std::endl;
  
  if(vect1.size() == 1 && vect2.size() == 0) vect2.push_back("");

  //std::cout << __LINE__ << std::endl;

  int retVal = 0;
  retVal += l1OfflineSubtract(vect1, vect2, std::stoi(argv[3]), argv[4], subTypeStr);
  return retVal;
}

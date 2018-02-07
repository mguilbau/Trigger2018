#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1F.h"
#include "TMath.h"
#include "TNamed.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"

#include "include/L1AnalysisEventDataFormat.h"
#include "include/L1AnalysisL1UpgradeDataFormat.h"
#include "include/mntToXRootdFileString.h"
#include "include/doGlobalDebug.h"
#include "include/runLumiEvtKey.h"
#include "include/returnRootFileContentsList.h"
#include "include/getLinBins.h"
#include "include/plotUtilities.h"
#include "include/histDefUtility.h"
#include "include/etaPhiFunc.h"
#include "include/kirchnerPalette.h"
#include "include/inToOutFileString.h"
#include "include/textFileToVector.h"
#include "include/cppWatch.h"

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
  const Int_t nValidAlgo = 9;
  const std::string validAlgo[nValidAlgo] = {"None", "Donut", "ChunkyDonut", "PhiRingPPExclude", "PhiRingPPTower", "PhiRingPPTowerMask", "PhiRingPP", "PhiRingHITower", "PhiRingHIRegion"};

  const Int_t nParam = 1;
  const std::string param[nParam] = {"SeedThresh2"};

  Int_t algoPos = -1;
  for(Int_t i = 0; i < nValidAlgo; ++i){
    if(inFileName.find(validAlgo[i]) != std::string::npos){
      algoPos = i;
      break;
    }
  }

  Int_t paramPos = -1;
  for(Int_t i = 0; i < nParam; ++i){
    if(inFileName.find(param[i]) != std::string::npos){
      paramPos = i;
      break;
    }
  }

  std::string outStr = "NOVALIDALGO";
  if(algoPos == 1 && inFileName.find("ChunkyDonut") != std::string::npos) outStr = "ChunkyDonut";
  else if(algoPos != -1) outStr = validAlgo[algoPos];

  if(paramPos != -1) outStr = outStr + "_" + param[paramPos];

  return outStr;
}


void turnOnToCanv(TFile* inFile_p, TH1* inHistDenom_p, std::vector<TH1*> inHistNum_p, std::string ptThresh, std::string jetLabel, std::string jetLabel2)
{
  if(inHistNum_p.size() == 0){
    std::cout << "WARNING: input vector of numerators is empty. return" << std::endl;
    return;
  }

  kirchnerPalette col;
  TLegend* leg_p = new TLegend(.6, .2, .9, .5);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(14);
  leg_p->SetBorderSize(0);
  leg_p->SetFillStyle(0);

  inFile_p->cd();
  TDatime* date = new TDatime();
  TH1F* hist_p = new TH1F("hist_p", ";Raw jet p_{T} (GeV/c);Efficiency", inHistDenom_p->GetNbinsX(), inHistDenom_p->GetBinLowEdge(1), inHistDenom_p->GetBinLowEdge(inHistDenom_p->GetNbinsX()+1));
  centerTitles({hist_p});
  hist_p->SetMinimum(0.0);
  hist_p->SetMaximum(1.1);

  std::string histName = inHistNum_p.at(0)->GetName();
  std::string canvStr = "canv";
  for(unsigned int i = 0; i < inHistNum_p.size(); ++i){
    canvStr = canvStr + "_" + getL1AlgoFromFileName(histName);
  }

  if(histName.find("Match") != std::string::npos) canvStr = canvStr + "_DRMatch";
  if(histName.find("Cent") != std::string::npos){
    std::string tempCentStr = histName.substr(histName.find("Cent"), histName.size() - histName.find("Cent"));
    tempCentStr = tempCentStr.substr(0, tempCentStr.find("_"));
    canvStr = canvStr + "_" + tempCentStr;
  }

  canvStr = canvStr + "_Pt" + ptThresh + "_TurnOn_c";
  TCanvas* canv_p = new TCanvas(canvStr.c_str(), canvStr.c_str(), 500, 500);
  gStyle->SetOptStat(0);
  canv_p->SetLeftMargin(canv_p->GetLeftMargin()*1.5);
  canv_p->SetBottomMargin(canv_p->GetLeftMargin());
  canv_p->SetTopMargin(canv_p->GetLeftMargin()/2.);
  canv_p->SetRightMargin(0.01);

  hist_p->DrawCopy("HIST");
  const Int_t nA = inHistNum_p.size();
  TGraphAsymmErrors* aPt_p[nA];
  TH1F* aPt_Leg_p[nA];

  for(unsigned int i = 0; i < inHistNum_p.size(); ++i){
    aPt_Leg_p[i] = new TH1F();
    aPt_p[i] = new TGraphAsymmErrors();
    aPt_p[i]->BayesDivide(inHistNum_p[i], inHistDenom_p);
    aPt_p[i]->SetMarkerSize(.6);
    aPt_p[i]->SetMarkerColor(col.getColor(i));
    aPt_p[i]->SetLineColor(col.getColor(i));
    aPt_p[i]->SetMarkerStyle(20);

    aPt_Leg_p[i]->SetMarkerSize(.6);
    aPt_Leg_p[i]->SetMarkerColor(col.getColor(i));
    aPt_Leg_p[i]->SetLineColor(col.getColor(i));
    aPt_Leg_p[i]->SetMarkerStyle(20);

    leg_p->AddEntry(aPt_Leg_p[i], getL1AlgoFromFileName(inHistNum_p[i]->GetName()).c_str(), "P L");

    aPt_p[i]->Draw("P");
  }

  leg_p->Draw("SAME");

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(14);
  label_p->SetNDC();
  label_p->DrawLatex(.15, .95, jetLabel.c_str());
  delete label_p;

  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);
  line_p->DrawLine(inHistDenom_p->GetBinLowEdge(1), 1., inHistDenom_p->GetBinLowEdge(inHistDenom_p->GetNbinsX()+1), 1.);
  delete line_p;

  canv_p->SaveAs(("pdfDir/" + canvStr + "_" + jetLabel2 + "_" + std::to_string(date->GetDate()) + ".pdf").c_str());

  delete canv_p;

  delete leg_p;

  for(unsigned int i = 0; i < inHistNum_p.size(); ++i){
    delete aPt_p[i];
  }

  delete hist_p;

  delete date;


  return;
}


int l1Comp(const std::string inForestName, std::vector<std::string> inL1AlgoNames)
{
  //handling the output file in event none is specified
  const Int_t nInStr = 1+inL1AlgoNames.size();
  std::string inStr[nInStr];
  for(Int_t i = 0; i < nInStr-1; ++i){inStr[i] = inL1AlgoNames.at(i);}
  inStr[nInStr-1] = inForestName;
  
  for(Int_t i = 0; i < nInStr; ++i){
    while(inStr[i].find("/") != std::string::npos){inStr[i].replace(0, inStr[i].find("/")+1, "");}
    while(inStr[i].find(".root") != std::string::npos){inStr[i].replace(inStr[i].find(".root"), std::string(".root").size(), "");}
  }
  std::string outFileName = "L1ALGOCOMP_";
  for(Int_t i = 0; i < nInStr-1; ++i){
    outFileName = outFileName + getL1AlgoFromFileName(inStr[i]) + "_";
  }
  outFileName.replace(outFileName.size()-1, 1, "");
  outFileName = "output/" + inToOutFileString(inForestName, outFileName);
  
  //Quickly grab number of jet algorithms + tree names in forest
  TFile* forestFile_p = TFile::Open(mntToXRootdFileString(inForestName).c_str(), "READ");
  std::vector<std::string> jetTreeStr = removeDuplicates(returnRootFileContentsList(forestFile_p, "TTree", "JetAnalyzer"));
  forestFile_p->Close();
  delete forestFile_p;

  const Int_t nJetAlgos = (Int_t)jetTreeStr.size();
  std::string jetAlgos[nJetAlgos];
  for(Int_t jI = 0; jI < nJetAlgos; ++jI){
    jetAlgos[jI] = jetTreeStr.at(jI).substr(0, jetTreeStr.at(jI).find("Jet"));
  }

  //number of algos to compare will always be two but just maintain a baseline
  const Int_t nL1Algo = inL1AlgoNames.size();
  
  const Int_t nL1JetThresholds = 5;
  const Float_t l1JetThresholds[nL1JetThresholds] = {8., 16., 24., 32., 40.};

  std::string l1AlgoStr[nL1Algo];
  runLumiEvtKey* l1AlgoMap[nL1Algo];
  for(Int_t i = 0; i < nL1Algo; ++i){
    l1AlgoStr[i] = getL1AlgoFromFileName(inL1AlgoNames.at(i));

    l1AlgoMap[i] = NULL;
    l1AlgoMap[i] = new runLumiEvtKey();
  }
  
  TFile* outFile_p = TFile::Open(mntToXRootdFileString(outFileName).c_str(), "RECREATE");

  const Int_t nCentBins = 3;
  const Int_t centBinsLow[nCentBins] = {0, 10, 30};
  const Int_t centBinsHi[nCentBins] = {10, 30, 100};

  const Int_t nL1PtBins2D = 3;
  const Double_t l1PtBins2DLow[nL1PtBins2D+1] = { 8.,  8., 24.,  40.};
  const Double_t l1PtBins2DHi[nL1PtBins2D+1] = {200., 24., 40., 200.};

  TH1F* l1JetEt_h[nL1Algo][nCentBins];
  TH1F* l1JetEta_h[nL1Algo][nCentBins][nL1PtBins2D+1];
  TH1F* l1JetPhi_h[nL1Algo][nCentBins][nL1PtBins2D+1];
  TH1F* jetPt_h[nJetAlgos][nCentBins];
  TH1F* jetPt_Trig_h[nJetAlgos][nL1Algo][nCentBins][nL1JetThresholds];
  TH1F* jetPt_MatchTrig_h[nJetAlgos][nL1Algo][nCentBins][nL1JetThresholds];
  TH1F* nJet_MissedTrig_h[nJetAlgos][nL1Algo][nCentBins][nL1JetThresholds];
  TH1F* hiBin_MissedTrig_h[nJetAlgos][nL1Algo][nL1JetThresholds];
  TH1F* allJetEta_MissedTrig_h[nJetAlgos][nL1Algo][nCentBins][nL1JetThresholds];

  for(Int_t i = 0; i < nL1Algo; ++i){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      const std::string centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";
   
      l1JetEt_h[i][cI] = new TH1F(("l1JetEt_" + l1AlgoStr[i] + "_" + centStr + "_h").c_str(), (";L1 Jet E_{T} (" + l1AlgoStr[i] + ");Counts (" + centStr2 + ")").c_str(), 50, 20, 120);
      centerTitles({l1JetEt_h[i][cI]});
      
      for(Int_t j = 0; j < nL1PtBins2D+1; ++j){
	std::string ptStr = "Pt" + prettyString(l1PtBins2DLow[j], 1, true) + "to" + prettyString(l1PtBins2DHi[j], 1, true);
	std::string ptStr2 = prettyString(l1PtBins2DLow[j], 1,false) + " < p_{T,L1} < " + prettyString(l1PtBins2DHi[j], 1, false);
	
	l1JetEta_h[i][cI][j] = new TH1F(("l1JetEta_" + l1AlgoStr[i] + "_" + centStr + "_" + ptStr + "_h").c_str(), (";L1 Jet #eta (" + l1AlgoStr[i] + ", " + ptStr2 + ");Counts (" + centStr2 + ")").c_str(), 50, -5.1, 5.1);
	l1JetPhi_h[i][cI][j] = new TH1F(("l1JetPhi_" + l1AlgoStr[i] + "_" + centStr + "_" + ptStr + "_h").c_str(), (";L1 Jet #phi (" + l1AlgoStr[i] + ", " + ptStr2 + ");Counts (" + centStr2 + ")").c_str(), 50, -TMath::Pi(), TMath::Pi());
	
	centerTitles({l1JetEta_h[i][cI][j], l1JetPhi_h[i][cI][j]});
      }
    }
  }

  const Int_t nJtPtBins = 19;
  const Float_t jtPtLow = 25;
  const Float_t jtPtHi = 120;
  Double_t jtPtBins[nJtPtBins+1];
  getLinBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);

  for(Int_t jI = 0; jI < nJetAlgos; ++jI){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      const std::string centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";

      jetPt_h[jI][cI] = new TH1F(("jetPt_" + jetAlgos[jI] + "_" + centStr + "_h").c_str(), (";Jet p_{T} (" + jetAlgos[jI] + ");Counts (" + centStr2 + ")").c_str(), nJtPtBins, jtPtBins);
      centerTitles({jetPt_h[jI][cI]});

      for(Int_t lI = 0; lI < nL1Algo; ++lI){
	for(Int_t tI = 0; tI < nL1JetThresholds; ++tI){
	  jetPt_Trig_h[jI][lI][cI][tI] = new TH1F(("jetPt_" + jetAlgos[jI] + "_" + l1AlgoStr[lI] + "_" + centStr + "_Pt" + prettyString(l1JetThresholds[tI], 1, true) + "_h").c_str(), (";Jet p_{T} w/ Trigger " + l1AlgoStr[lI] + ", L1 p_{T} > " + prettyString(l1JetThresholds[tI], 1, false) + " (" + jetAlgos[jI] + ");Counts (" + centStr2 + ")").c_str(), nJtPtBins, jtPtBins);
	  
	  jetPt_MatchTrig_h[jI][lI][cI][tI] = new TH1F(("jetPt_Match_" + jetAlgos[jI] + "_" + l1AlgoStr[lI] + "_" + centStr + "_Pt" + prettyString(l1JetThresholds[tI], 1, true) + "_h").c_str(), (";Jet p_{T} w/ Trigger " + l1AlgoStr[lI] + ", L1 p_{T} > " + prettyString(l1JetThresholds[tI], 1, false) + " (" + jetAlgos[jI] + ");Counts (" + centStr2 + ")").c_str(), nJtPtBins, jtPtBins);
	  
	  nJet_MissedTrig_h[jI][lI][cI][tI] = new TH1F(("nJet_Missed_" + jetAlgos[jI] + "_" + l1AlgoStr[lI] + "_" + centStr + "_Pt" + prettyString(l1JetThresholds[tI], 1, true) + "_h").c_str(), (";nJets for missed w/" + l1AlgoStr[lI] + ", L1 p_{T} > " + prettyString(l1JetThresholds[tI], 1, false) + " (" + jetAlgos[jI] + ");Counts (" + centStr2 + ")").c_str(), 15, -1.5, 13.5);

	  if(cI == 0) hiBin_MissedTrig_h[jI][lI][tI] = new TH1F(("hiBin_Missed_" + jetAlgos[jI] + "_" + l1AlgoStr[lI] + "_Pt" + prettyString(l1JetThresholds[tI], 1, true) + "_h").c_str(), (";hiBin for missed w/" + l1AlgoStr[lI] + ", L1 p_{T} > " + prettyString(l1JetThresholds[tI], 1, false) + " (" + jetAlgos[jI] + ");Counts").c_str(), 200, -0.5, 199.5);
	  
	  allJetEta_MissedTrig_h[jI][lI][cI][tI] = new TH1F(("allJetEta_Missed_" + jetAlgos[jI] + "_" + l1AlgoStr[lI] + "_" + centStr + "_Pt" + prettyString(l1JetThresholds[tI], 1, true) + "_h").c_str(), (";All L1 Jet #eta for missed w/" + l1AlgoStr[lI] + ", L1 p_{T} > " + prettyString(l1JetThresholds[tI], 1, false) + " (" + jetAlgos[jI] + ");Counts (" + centStr2 + ")").c_str(), 201, -100.5, 100.5);
	  
	  centerTitles({jetPt_Trig_h[jI][lI][cI][tI], jetPt_MatchTrig_h[jI][lI][cI][tI], nJet_MissedTrig_h[jI][lI][cI][tI], hiBin_MissedTrig_h[jI][lI][tI], allJetEta_MissedTrig_h[jI][lI][cI][tI]});
	}
      }
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << "Processing l1 algos..." << std::endl;

  TFile* inL1Algo_p[nL1Algo];
  TTree* inL1AlgoEvtTree_p[nL1Algo];
  TTree* inL1AlgoUpgradeTree_p[nL1Algo];

  L1Analysis::L1AnalysisEventDataFormat* evt[nL1Algo];
  L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade[nL1Algo];
  TBranch* b_upgrade[nL1Algo];

  for(Int_t i = 0; i < nL1Algo; ++i){
    std::cout << " Input " << i << ": " << inL1AlgoNames.at(i) << std::endl;

    inL1Algo_p[i] = NULL;
    inL1Algo_p[i] = TFile::Open(mntToXRootdFileString(inL1AlgoNames.at(i)).c_str(), "READ");
    inL1AlgoEvtTree_p[i] = (TTree*)inL1Algo_p[i]->Get("l1EventTree/L1EventTree");
    inL1AlgoUpgradeTree_p[i] = (TTree*)inL1Algo_p[i]->Get("l1UpgradeEmuTree/L1UpgradeTree");
    inL1AlgoEvtTree_p[i]->SetMaxVirtualSize(4000000000);
    inL1AlgoUpgradeTree_p[i]->SetMaxVirtualSize(4000000000);

    evt[i] = new L1Analysis::L1AnalysisEventDataFormat();
    upgrade[i] = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
    b_upgrade[i]=NULL;

    inL1AlgoEvtTree_p[i]->SetBranchStatus("*", 0);
    inL1AlgoEvtTree_p[i]->SetBranchStatus("Event", 1);
    inL1AlgoEvtTree_p[i]->SetBranchStatus("run", 1);
    inL1AlgoEvtTree_p[i]->SetBranchStatus("lumi", 1);
    inL1AlgoEvtTree_p[i]->SetBranchStatus("event", 1);

    inL1AlgoEvtTree_p[i]->SetBranchAddress("Event", &(evt[i]));


    inL1AlgoUpgradeTree_p[i]->SetBranchStatus("*", 0);
    inL1AlgoUpgradeTree_p[i]->SetBranchStatus("L1Upgrade", 1);
    inL1AlgoUpgradeTree_p[i]->SetBranchStatus("nJets", 1);
    inL1AlgoUpgradeTree_p[i]->SetBranchStatus("jetEt", 1);
    inL1AlgoUpgradeTree_p[i]->SetBranchStatus("jetEta", 1);
    inL1AlgoUpgradeTree_p[i]->SetBranchStatus("jetIEta", 1);
    inL1AlgoUpgradeTree_p[i]->SetBranchStatus("jetPhi", 1);

    inL1AlgoUpgradeTree_p[i]->SetBranchAddress("L1Upgrade", &(upgrade[i]), &(b_upgrade[i]));
    //    b_upgrade[i]->LoadBaskets();
    //    inL1AlgoUpgradeTree_p[i]->GetBranch("L1Upgrade")->LoadBaskets();
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


  for(Int_t i = 0; i < nL1Algo; ++i){
    std::cout << "Processing " << i << ": " << inL1AlgoNames.at(i) << std::endl;

    for(Int_t entry = 0; entry < inL1AlgoEvtTree_p[i]->GetEntries(); ++entry){
      if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << inL1AlgoEvtTree_p[i]->GetEntries() << std::endl;

      inL1AlgoEvtTree_p[i]->GetEntry(entry);

      l1AlgoMap[i]->addKey(evt[i]->run, evt[i]->lumi, evt[i]->event, entry);
    }
  }

  std::cout << "Processing forests..." << std::endl;

  UInt_t runF, lumiF;
  ULong64_t eventF;
  Float_t vz_;
  Int_t hiBin_;

  Int_t pBeamScrapingFilter_;
  Int_t pPAprimaryVertexFilter_;
  Int_t HBHENoiseFilterResultRun2Loose_;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t nMaxJets = 500;
  Int_t nref_[nJetAlgos];
  Float_t jtpt_[nJetAlgos][nMaxJets];
  Float_t rawpt_[nJetAlgos][nMaxJets];
  Float_t jteta_[nJetAlgos][nMaxJets];
  Float_t jtphi_[nJetAlgos][nMaxJets];

  Float_t jtPfCHF_[nJetAlgos][nMaxJets];
  Float_t jtPfCEF_[nJetAlgos][nMaxJets];
  Float_t jtPfNHF_[nJetAlgos][nMaxJets];
  Float_t jtPfNEF_[nJetAlgos][nMaxJets];
  Float_t jtPfMUF_[nJetAlgos][nMaxJets];

  forestFile_p = TFile::Open(mntToXRootdFileString(inForestName).c_str(), "READ");  
  TTree* hiTree_p = (TTree*)forestFile_p->Get("hiEvtAnalyzer/HiTree");
  hiTree_p->SetMaxVirtualSize(4000000000);
  TTree* skimTree_p = (TTree*)forestFile_p->Get("skimanalysis/HltTree");
  skimTree_p->SetMaxVirtualSize(4000000000);
  TTree* jetTree_p[nJetAlgos];

  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("run", 1);
  hiTree_p->SetBranchStatus("lumi", 1);
  hiTree_p->SetBranchStatus("evt", 1);
  hiTree_p->SetBranchStatus("vz", 1);
  hiTree_p->SetBranchStatus("hiBin", 1);

  hiTree_p->SetBranchAddress("run", &runF);
  hiTree_p->SetBranchAddress("lumi", &lumiF);
  hiTree_p->SetBranchAddress("evt", &eventF);
  hiTree_p->SetBranchAddress("vz", &vz_);
  hiTree_p->SetBranchAddress("hiBin", &hiBin_);

  skimTree_p->SetBranchStatus("*", 0);
  skimTree_p->SetBranchStatus("pBeamScrapingFilter", 1);
  skimTree_p->SetBranchStatus("pPAprimaryVertexFilter", 1);
  skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);

  skimTree_p->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter_);
  skimTree_p->SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter_);
  skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int i = 0; i < jetTreeStr.size(); ++i){
    jetTree_p[i] = (TTree*)forestFile_p->Get(jetTreeStr.at(i).c_str());
    jetTree_p[i]->SetMaxVirtualSize(4000000000);

    jetTree_p[i]->SetBranchStatus("*", 0);
    jetTree_p[i]->SetBranchStatus("nref", 1);
    jetTree_p[i]->SetBranchStatus("jtpt", 1);
    jetTree_p[i]->SetBranchStatus("rawpt", 1);
    jetTree_p[i]->SetBranchStatus("jteta", 1);
    jetTree_p[i]->SetBranchStatus("jtphi", 1);

    jetTree_p[i]->SetBranchAddress("nref", &nref_[i]);
    jetTree_p[i]->SetBranchAddress("jtpt", jtpt_[i]);
    jetTree_p[i]->SetBranchAddress("rawpt", rawpt_[i]);
    jetTree_p[i]->SetBranchAddress("jteta", jteta_[i]);
    jetTree_p[i]->SetBranchAddress("jtphi", jtphi_[i]);

    if(jetTreeStr.at(i).find("Cs") != std::string::npos || jetTreeStr.at(i).find("CS") != std::string::npos){
      jetTree_p[i]->SetBranchStatus("jtPfCHF", 1);
      jetTree_p[i]->SetBranchStatus("jtPfCEF", 1);
      jetTree_p[i]->SetBranchStatus("jtPfNHF", 1);
      jetTree_p[i]->SetBranchStatus("jtPfNEF", 1);
      jetTree_p[i]->SetBranchStatus("jtPfMUF", 1);

      jetTree_p[i]->SetBranchAddress("jtPfCHF", jtPfCHF_[i]);
      jetTree_p[i]->SetBranchAddress("jtPfCEF", jtPfCEF_[i]);
      jetTree_p[i]->SetBranchAddress("jtPfNHF", jtPfNHF_[i]);
      jetTree_p[i]->SetBranchAddress("jtPfNEF", jtPfNEF_[i]);
      jetTree_p[i]->SetBranchAddress("jtPfMUF", jtPfMUF_[i]);
    }
  }

  const Int_t nEntriesForest = hiTree_p->GetEntries();

  Int_t totalFound = 0;

  std::cout << "Processing forest..." << std::endl;
  for(Int_t entry = 0; entry < nEntriesForest; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntriesForest << std::endl;

    skimTree_p->GetEntry(entry);
    hiTree_p->GetEntry(entry);

    if(TMath::Abs(vz_) > 15.) continue;
    if(!pPAprimaryVertexFilter_) continue;
    if(!pBeamScrapingFilter_) continue;
    if(!HBHENoiseFilterResultRun2Loose_) continue;

    bool isGood = true;
    Int_t entryL1Algos[nL1Algo];
    for(Int_t lI = 0; lI < nL1Algo; ++lI){entryL1Algos[lI] = -1;}
    for(Int_t lI = 0; lI < nL1Algo; ++lI){
      entryL1Algos[lI] = l1AlgoMap[lI]->getEntryFromKey(runF, lumiF, eventF);
      if(entryL1Algos[lI] < 0){isGood = false; break;}
    }
    if(!isGood) continue;
    
    Int_t centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(hiBin_/2 >= centBinsLow[cI] && hiBin_/2 < centBinsHi[cI]){
	centPos = cI;
	break;
      }
    }


    totalFound++;
    //    std::cout << "  Found entry!" << std::endl;

    for(unsigned int i = 0; i < jetTreeStr.size(); ++i){jetTree_p[i]->GetEntry(entry);}

    //    std::cout << runF << ", " << lumiF << ", " << eventF << std::endl;

    for(Int_t lI = 0; lI < nL1Algo; ++lI){
      //      std::cout << " " << entryL1Algos[lI] << ",";
      inL1AlgoUpgradeTree_p[lI]->GetEntry(entryL1Algos[lI]);
    
      for(unsigned int jI = 0; jI < upgrade[lI]->jetEt.size(); ++jI){
	l1JetEt_h[lI][centPos]->Fill(upgrade[lI]->jetEt.at(jI));
	
	for(Int_t ptlI = 0; ptlI < nL1PtBins2D+1; ++ptlI){
	  if(l1PtBins2DLow[ptlI] <= upgrade[lI]->jetEt.at(jI) && upgrade[lI]->jetEt.at(jI) < l1PtBins2DHi[ptlI]){
	    l1JetEta_h[lI][centPos][ptlI]->Fill(upgrade[lI]->jetEta.at(jI));
	    l1JetPhi_h[lI][centPos][ptlI]->Fill(upgrade[lI]->jetPhi.at(jI));
	  }
	}
      }
    }
    //    std::cout << std::endl;

    for(unsigned int i = 0; i < jetTreeStr.size(); ++i){
      //      Float_t leadingJtPt_ = -999;
      Float_t leadingRawPt_ = -999;
      Float_t leadingJtPhi_ = -999;
      Float_t leadingJtEta_ = -999;

      for(Int_t jI = 0; jI < nref_[i]; ++jI){
	if(TMath::Abs(jteta_[i][jI]) > 2.) continue;

	if(jetTreeStr.at(i).find("Cs") != std::string::npos || jetTreeStr.at(i).find("CS") != std::string::npos){
	  if(jtPfMUF_[i][jI] > .8) continue;
	  if(jtPfCEF_[i][jI] > .8) continue;
	}
	
	if(rawpt_[i][jI] > leadingRawPt_){
	  //	  leadingJtPt_ = jtpt_[i][jI];
	  leadingRawPt_ = rawpt_[i][jI];
	  leadingJtPhi_ = jtphi_[i][jI];
	  leadingJtEta_ = jteta_[i][jI];
	}
      }

      if(leadingRawPt_ > jtPtLow && leadingRawPt_ < jtPtHi){
	jetPt_h[i][centPos]->Fill(leadingRawPt_);

	Float_t unmatchedL1Pt_[nL1Algo];
	Float_t matchedL1Pt_[nL1Algo];
	for(Int_t lI = 0; lI < nL1Algo; ++lI){
	  unmatchedL1Pt_[lI] = -999.;
	  matchedL1Pt_[lI] = -999.;

       	  for(unsigned int jI = 0; jI < upgrade[lI]->jetEt.size(); ++jI){
	    if(unmatchedL1Pt_[lI] < upgrade[lI]->jetEt.at(jI)) unmatchedL1Pt_[lI] = upgrade[lI]->jetEt.at(jI);
	    if(getDR(upgrade[lI]->jetEta.at(jI), upgrade[lI]->jetPhi.at(jI), leadingJtEta_, leadingJtPhi_) > 1.5) continue;
	    if(matchedL1Pt_[lI] < upgrade[lI]->jetEt.at(jI)) matchedL1Pt_[lI] = upgrade[lI]->jetEt.at(jI);
	  }
	}

	for(Int_t lI = 0; lI < nL1Algo; ++lI){
	  for(Int_t tI = 0; tI < nL1JetThresholds; ++tI){
	    if(unmatchedL1Pt_[lI] > l1JetThresholds[tI]) jetPt_Trig_h[i][lI][centPos][tI]->Fill(leadingRawPt_);

	    if(matchedL1Pt_[lI] > l1JetThresholds[tI]) jetPt_MatchTrig_h[i][lI][centPos][tI]->Fill(leadingRawPt_);
	    else if(leadingRawPt_ > 80 && matchedL1Pt_[lI] < l1JetThresholds[tI] && lI == nL1Algo-1){
	      std::cout << "Missed by " << l1AlgoStr[lI] << ", " << l1JetThresholds[tI] << ": " << leadingRawPt_ << ", " << leadingJtPhi_ << ", " << leadingJtEta_ << std::endl;
	      std::cout << " Entry Forest: " << entry << std::endl;
	      std::cout << " Entry L1: " << entryL1Algos[lI] << std::endl;

	      nJet_MissedTrig_h[i][lI][centPos][tI]->Fill(upgrade[lI]->nJets);
	      hiBin_MissedTrig_h[i][lI][tI]->Fill(hiBin_);

	      for(unsigned int jI = 0; jI < upgrade[lI]->jetIEta.size(); ++jI){
		allJetEta_MissedTrig_h[i][lI][centPos][tI]->Fill(upgrade[lI]->jetIEta.at(jI));
	      }

	    }
	  }
	}
      }
    }
  }

  std::cout << "Total found: " << totalFound << std::endl;

  forestFile_p->Close();
  delete forestFile_p;

  for(Int_t lI = 0; lI < nL1Algo; ++lI){
    inL1Algo_p[lI]->Close();
    delete inL1Algo_p[lI];
  }
  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  outFile_p->cd();
  for(Int_t i = 0; i < nL1Algo; ++i){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      l1JetEt_h[i][cI]->SetMinimum(0.5);
      l1JetEt_h[i][cI]->Write("", TObject::kOverwrite);
      delete l1JetEt_h[i][cI];
      
      for(Int_t lI = 0; lI < nL1PtBins2D+1; ++lI){
	l1JetEta_h[i][cI][lI]->SetMinimum(0.5);
	l1JetPhi_h[i][cI][lI]->SetMinimum(0.5);    
	
	l1JetEta_h[i][cI][lI]->Write("", TObject::kOverwrite);
	delete l1JetEta_h[i][cI][lI];
	
	l1JetPhi_h[i][cI][lI]->Write("", TObject::kOverwrite);
	delete l1JetPhi_h[i][cI][lI];    
      }
    }
  }

  for(Int_t jI = 0; jI < nJetAlgos; ++jI){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      const std::string centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";

      jetPt_h[jI][cI]->Write("", TObject::kOverwrite);

      for(Int_t tI = 0; tI < nL1JetThresholds; ++tI){
	std::vector<TH1*> temp;
	std::vector<TH1*> tempMatch;
	for(Int_t lI = 0; lI < nL1Algo; ++lI){
	  jetPt_Trig_h[jI][lI][cI][tI]->Write("", TObject::kOverwrite);
	  jetPt_MatchTrig_h[jI][lI][cI][tI]->Write("", TObject::kOverwrite);
	  nJet_MissedTrig_h[jI][lI][cI][tI]->Write("", TObject::kOverwrite);
	  if(cI == 0) hiBin_MissedTrig_h[jI][lI][tI]->Write("", TObject::kOverwrite);
	  allJetEta_MissedTrig_h[jI][lI][cI][tI]->Write("", TObject::kOverwrite);
	  
	  temp.push_back(jetPt_Trig_h[jI][lI][cI][tI]);
	  tempMatch.push_back(jetPt_MatchTrig_h[jI][lI][cI][tI]);
	}
	
	std::string jetLabel = "Offline Jet " + jetAlgos[jI] + "; L1 Algo. p_{T} > " + prettyString(l1JetThresholds[tI], 0, false) + " (" + centStr2 + ")";
	std::string jetLabelMatch = "Offline Jet " + jetAlgos[jI] + "; L1 Algo. p_{T} > " + prettyString(l1JetThresholds[tI], 0, false) + ", #DeltaR Match (" + centStr2 + ")";
	turnOnToCanv(outFile_p, jetPt_h[jI][cI], temp, prettyString(l1JetThresholds[tI], 0, true), jetLabel, jetAlgos[jI]);
	turnOnToCanv(outFile_p, jetPt_h[jI][cI], tempMatch, prettyString(l1JetThresholds[tI], 0, true), jetLabelMatch, jetAlgos[jI]);
	temp.clear();
	tempMatch.clear();
	
	for(Int_t lI = 0; lI < nL1Algo; ++lI){
	  delete nJet_MissedTrig_h[jI][lI][cI][tI];
	  if(cI == 0) delete hiBin_MissedTrig_h[jI][lI][tI];
	  delete allJetEta_MissedTrig_h[jI][lI][cI][tI];
	  delete jetPt_Trig_h[jI][lI][cI][tI];
	  delete jetPt_MatchTrig_h[jI][lI][cI][tI];
	}
      }
      
      delete jetPt_h[jI][cI];
    }
  }

  outFile_p->Close();
  delete outFile_p;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  cppWatch timer;
  timer.start();

  if(argc < 3){
    std::cout << "Usage: ./l1Comp.exe <inForestName> <endlessInL1Algos>" << std::endl;
    return 1;
  }

  std::vector<std::string> inputs;
  for(Int_t i = 2; i < argc; ++i){
    std::string fileName = argv[i];
    if(fileName.find(".root") != std::string::npos) inputs.push_back(fileName);
    else if(fileName.find(".txt") != std::string::npos) textFileToVectorAppend(&inputs, fileName, ".root");
  }

  std::cout << "INPUTS: " << std::endl;
  for(unsigned int i = 0; i < inputs.size(); ++i){
    std::cout << " " << i << ": " << inputs.at(i) << std::endl;
  }

  if(inputs.size() == 0){
    std::cout << "WARNING: Given set of inputs contains no valid .root files. return 1" << std::endl;
    return 1;
  }

  int retVal = l1Comp(argv[1], inputs);

  timer.stop();
  std::cout << "Timer: " << timer.total() << std::endl;

  return retVal;
}

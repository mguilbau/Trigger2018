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

#include "include/checkMakeDir.h"
#include "include/doGlobalDebug.h"
#include "include/L1Tools.h"
#include "include/mntToXRootdFileString.h"
#include "include/L1AnalysisEventDataFormat.h"
#include "include/L1AnalysisL1CaloTowerDataFormat.h"
#include "include/L1AnalysisL1UpgradeDataFormat.h"

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

int l1OfflineSubtract(const std::string inFileName, const std::string subtractType)
{
  const Int_t nSubType = 5;
  const std::string subType[nSubType] = {"None", "PhiRingHITower", "PhiRingPPTower", "PhiRingHIRegion", "ChunkyDonut"};

  bool isFound = false;
  for(Int_t i = 0; i < nSubType; ++i){
    if(subtractType.size() == subType[i].size() && subtractType.find(subType[i]) != std::string::npos){
      isFound = true;
      break;
    }
  }

  if(!isFound){
    std::cout << "Given subtraction type \'" << subtractType << "\' is not found. Please pick one of: " << std::endl;
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

  const int minSeedThresh = 8;
  const int nIEta = 83;
  const int nIPhi = 72;

  TFile* outFile_p = new TFile(("output/l1OfflineSubtract_" + subtractType + "_" + dateStr + ".root").c_str(), "RECREATE");
  TH1F* etaMiss_p = new TH1F("etaMiss_h", ";HW #eta;Counts (Misses)", 103, -51.5, 51.5);
  TH1F* phiMiss_p = new TH1F("phiMiss_h", ";HW #phi;Counts (Misses)", 150, -0.5, 149.5);

  TH1F* etaMiss_CMSSWJet_p = new TH1F("etaMiss_CMSSWJet_h", ";HW #eta;Counts (Misses)", 103, -51.5, 51.5);
  TH1F* phiMiss_CMSSWJet_p = new TH1F("phiMiss_CMSSWJet_h", ";HW #phi;Counts (Misses)", 150, -0.5, 149.5);

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  TFile* inFile_p = TFile::Open(mntToXRootdFileString(inFileName).c_str(), "READ");
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

  const Int_t nEntries = TMath::Min(100000, (Int_t)l1CaloTree_p->GetEntries());

  std::map<Int_t, Int_t> max;
  for(Int_t i = 1; i <= 41; ++i){
    max[i] = 0;
    max[-i] = 0;
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  int foundCount = 0;
  int missCount = 0;

  std::map<int, std::map<int, int> > maxTowIEt;
  for(int etaI = -41; etaI <= 41; ++etaI){
    for(int phiI = 0; phiI <= nIPhi; ++phiI){
      (maxTowIEt[etaI])[phiI] = 0;
    }
  }

  bool isNone = subtractType.size() == std::string("None").size() && subtractType.find("None") != std::string::npos;
  bool isChunkyDonut = subtractType.size() == std::string("ChunkyDonut").size() && subtractType.find("ChunkyDonut") != std::string::npos;
  bool isPhiRingPPTower = subtractType.size() == std::string("PhiRingPPTower").size() && subtractType.find("PhiRingPPTower") != std::string::npos;
  bool isPhiRingHITower = subtractType.size() == std::string("PhiRingHITower").size() && subtractType.find("PhiRingHITower") != std::string::npos;
  bool isPhiRingHIRegion = subtractType.size() == std::string("PhiRingHIRegion").size() && subtractType.find("PhiRingHIRegion") != std::string::npos;

  std::cout << "Processing " << nEntries << "..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;

    l1CaloTree_p->GetEntry(entry);
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    inL1AlgoUpgradeTree_p->GetEntry(entry);
    inL1AlgoEvtTree_p->GetEntry(entry);
    
    //    if(entry != 22117) continue;
    /*
    if(evt->run != 1) continue;
    if(evt->lumi != 126) continue;
    if(evt->event < 12500 || evt->event > 12520) continue;
    */  

    //    std::cout << "  Event: " << evt->event << std::endl;

    //    if(entry != 190099) continue;
    //if(evt->run != 1) continue;
    //if(evt->lumi != 2541) continue;
    //if(evt->event < 254000) continue;
    //if(evt->event > 254100) continue;
    //if(evt->event != 254046) continue;

    int leadJtIEt = -999;
    int leadJtIEta = -999;
    int leadJtIPhi = -999;
    int leadJtSeedPt = -999;

    for(unsigned int jI = 0; jI < upgrade->jetIEt.size(); ++jI){
      //      if(TMath::Abs(upgrade->jetIEta.at(jI)) >= gtEta(25)) continue;
      if(upgrade->jetIEt.at(jI) > leadJtIEt){
	leadJtIEt = upgrade->jetIEt.at(jI);
	leadJtIEta = upgrade->jetIEta.at(jI);
	leadJtIPhi = upgrade->jetIPhi.at(jI);
	leadJtSeedPt = upgrade->jetSeedEt.at(jI);
      }
    }

    if(leadJtIEt == -999) leadJtIEt = 0;
  
    std::map<int, std::map<int, int> > newTowIEt;
    for(int etaI = -41; etaI <= 41; ++etaI){
      for(int phiI = 0; phiI <= nIPhi; ++phiI){
	(newTowIEt[etaI])[phiI] = 0;
      }
    }

    std::map<Int_t, Int_t> uePerEta;
    std::map<Int_t, Int_t> tempMax;
    for(Int_t i = 1; i <= 41; ++i){
      uePerEta[i] = 0;
      uePerEta[-i] = 0;

      tempMax[i] = 0;
      tempMax[-i] = 0;
    }


    for(unsigned int i = 0; i < towers_->iet.size(); ++i){
      //     std::cout << towers_->iet.at(i) << ", " << towers_->ieta.at(i) << ", " << towers_->iphi.at(i) << std::endl;

      if(towers_->ieta[i] == -25) continue;
      if(towers_->ieta[i] == -26) continue;
      if(towers_->ieta[i] == -27) continue;
      if(towers_->ieta[i] == -28) continue;

      if(towers_->ieta[i] == 25) continue;
      if(towers_->ieta[i] == 26) continue;
      if(towers_->ieta[i] == 27) continue;
      if(towers_->ieta[i] == 28) continue;

      uePerEta[towers_->ieta[i]] += towers_->iet[i];
      tempMax[towers_->ieta[i]] += 1;
    }
    
    for(Int_t etaI = 29; etaI < 41; ++etaI){
      uePerEta[-etaI] = uePerEta[-(etaI+1)];
      uePerEta[etaI] = uePerEta[etaI+1];
    }
    uePerEta[-41] = 0;
    uePerEta[41] = 0;
    

    for(Int_t phiI = 0; phiI <= nIPhi; ++phiI){
      uePerEta[-28] = 0;
      uePerEta[-27] = 0;
      uePerEta[-26] = 0;
      uePerEta[-25] = 0;
      
      uePerEta[28] = 0;
      uePerEta[27] = 0;
      uePerEta[26] = 0;
      uePerEta[25] = 0;	
    }
    
    if(isNone || isChunkyDonut){
      for(Int_t i = 1; i <= 41; ++i){
	uePerEta[i] = 0;
	uePerEta[-i] = 0;
      }
    }
    else if(isPhiRingHITower){
      for(unsigned int i = 1; i < 41; ++i){
	if(i%2 != 0) continue;
	Float_t tempVal = (uePerEta[i-1] + uePerEta[i])/2;
	uePerEta[i-1] = tempVal;
	uePerEta[i] = tempVal;
	
	tempVal = (uePerEta[-(i-1)] + uePerEta[-i])/2;
	uePerEta[-(i-1)] = tempVal;
	uePerEta[-i] = tempVal;
      }
    }
    else if(isPhiRingHIRegion){
      for(unsigned int i = 1; i < 41; ++i){
	if(i%4 == 3) continue;
	if(i%4 == 0) continue;
	if(i%4 == 2) continue;

	Float_t tempVal = (uePerEta[i] + uePerEta[i+1] + uePerEta[i+2] + uePerEta[i+3])/4;
	uePerEta[i] = tempVal;
	uePerEta[i+1] = tempVal;
	uePerEta[i+2] = tempVal;
	uePerEta[i+3] = tempVal;
	
	tempVal = (uePerEta[-i] + uePerEta[-i-1] + uePerEta[-i-2] + uePerEta[-i-3])/4;
	uePerEta[-i] = tempVal;
	uePerEta[-i-1] = tempVal;
	uePerEta[-i-2] = tempVal;
	uePerEta[-i-3] = tempVal;
      }
    }


    for(unsigned int i = 0; i < towers_->iet.size(); ++i){
      int puEtaPos = towers_->ieta[i];
      if(puEtaPos >= 29) puEtaPos--;
      else if(puEtaPos <= -29) puEtaPos++;

      int iEt = TMath::Max(0, towers_->iet.at(i) - ((int)(uePerEta[puEtaPos]/72)));

      //      int ietaPos = towers_->ieta[i] + 41;
      //      if(towers_->ieta[i] > 0) ietaPos -= 1;

      //      if(towers_->iphi[i] == 72) towers_->iphi[i] = 0;

      /*
      if(isPhiRingHITower || isPhiRingPPTower || isPhiRingHIRegion){
	if(towers_->ieta[i] <= -25 && towers_->ieta[i] >= -28) newTowIEt[ietaPos][towers_->iphi[i]] = 0;
	else if(towers_->ieta[i] <= 28 && towers_->ieta[i] >= 25) newTowIEt[ietaPos][towers_->iphi[i]] = 0;
	else newTowIEt[ietaPos][towers_->iphi[i]] = iEt;
      }
      */
     
      (newTowIEt[towers_->ieta[i]])[towers_->iphi[i]] = iEt;
      //      newTowIEt[ietaPos][towers_->iphi[i]] = iEt;
    }

        
    for(Int_t phiI = 0; phiI <= nIPhi; ++phiI){
      for(Int_t etaI = 29; etaI < 41; ++etaI){
	(newTowIEt[-etaI])[phiI] = (newTowIEt[-(etaI+1)])[phiI];
	(newTowIEt[etaI])[phiI] = (newTowIEt[etaI+1])[phiI];
      }
      (newTowIEt[-41])[phiI] = 0;
      (newTowIEt[41])[phiI] = 0;
    }
    for(Int_t phiI = 0; phiI <= nIPhi; ++phiI){
      (newTowIEt[-28])[phiI] = 0;
      (newTowIEt[-27])[phiI] = 0;
      (newTowIEt[-26])[phiI] = 0;
      (newTowIEt[-25])[phiI] = 0;
      
      (newTowIEt[28])[phiI] = 0;
      (newTowIEt[27])[phiI] = 0;
      (newTowIEt[26])[phiI] = 0;
      (newTowIEt[25])[phiI] = 0;	

    }


    for(Int_t phiI = 0; phiI <= nIPhi; ++phiI){
      for(Int_t etaI = -41; etaI <= 41; ++etaI){
	if((maxTowIEt[etaI])[phiI] < (newTowIEt[etaI])[phiI]) (maxTowIEt[etaI])[phiI] = (newTowIEt[etaI])[phiI];
      }
    }

    const std::string csvStr = "output/l1Towers_" + subtractType + "_Run" + std::to_string(evt->run) + "_Lumi" + std::to_string(evt->lumi) + "_Evt" + std::to_string(evt->event) + "_" + dateStr + ".csv";
    
    //    std::cout << csvStr << std::endl;
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
	disp = std::to_string((newTowIEt[etaI])[phiI]);
	if(disp.size() == 1) disp = " " + disp + ",";
	else if(disp.size() == 2) disp = disp + ",";

	file << disp;
	//std::cout << disp;
      }
      file << std::endl;
      //std::cout << std::endl;
    }    
    file.close();
  
    std::vector<int> jetPt;
    std::vector<int> jetPhi;
    std::vector<int> jetEta;
    std::vector<int> jetPhiPrev;
    std::vector<int> jetEtaPrev;
    std::vector<int> seedPt;
    std::vector< std::vector<int> > fullSeedPt;
    std::vector< std::vector<int> > fullSeedEta;
    std::vector< std::vector<int> > fullSeedPhi;

    for(int etaI = -41; etaI <= 41; ++etaI){

      if(etaI == 0) continue;

      if(etaI == -25) continue;
      if(etaI == -26) continue;
      if(etaI == -27) continue;
      if(etaI == -28) continue;

      if(etaI == 25) continue;
      if(etaI == 26) continue;
      if(etaI == 27) continue;
      if(etaI == 28) continue;
      
      for(int phiI = 0; phiI <= nIPhi; ++phiI){
	if((newTowIEt[etaI])[phiI] >= minSeedThresh){
	  int tempJetPt = 0;

	  std::vector<int> tempSeedPt;
	  std::vector<int> tempSeedPhi;
	  std::vector<int> tempSeedEta;

	  bool goodSeed = true;

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

	    for(Int_t phiI2 = phiI - 4; phiI2 <= phiI + 4; ++phiI2){ 	      
	      int phiPos = phiI2;
	      if(phiPos <= 0) phiPos += nIPhi;
	      if(phiPos > nIPhi) phiPos -= nIPhi;

	      tempJetPt += (newTowIEt[etaPos2])[phiPos];

	      tempSeedPt.push_back((newTowIEt[etaPos2])[phiPos]);
	      tempSeedPhi.push_back(phiPos);
	      tempSeedEta.push_back(etaPos2);
	      
	      if(etaPos2 == etaI && phiPos == phiI) continue;
	      else{
		int dPhi = phiI2 - phiI + 4;
		int dEta = etaI2 - etaI + 4;

		if(mask_[8-dPhi][dEta] == 1 && (newTowIEt[etaI])[phiI] < newTowIEt[etaPos2][phiPos]) goodSeed = false;
		else if(mask_[8-dPhi][dEta] == 2 && (newTowIEt[etaI])[phiI] <= (newTowIEt[etaPos2])[phiPos]) goodSeed = false;

		if(!goodSeed){
		  break;
		}
	      }
	    }
	    if(!goodSeed) break;
	  }

	  if(!goodSeed) continue;

	  if(isChunkyDonut){
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

		ptFlapPhiUp += (newTowIEt[etaPos])[phiPos];
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

		ptFlapPhiDown += (newTowIEt[etaPos])[phiPos];
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
		
		ptFlapEtaDown += (newTowIEt[etaPos])[phiPos];
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

		ptFlapEtaUp += (newTowIEt[etaPos])[phiPos];
	      }
	    }
	  
	    std::vector<int> flaps = {ptFlapEtaUp, ptFlapEtaDown, ptFlapPhiUp, ptFlapPhiDown};
	    if((newTowIEt[etaI])[phiI] == 17){
	      //	      std::cout << " Seedpt, phi eta: " << (newTowIEt[etaI])[phiI] << ", " << phiI << ", " << etaI << std::endl;
	      //	      std::cout << " Flaps etaup, etadown, phiup, phidown: " << flaps.at(0) << ", " << flaps.at(1) << ", " << flaps.at(2) << ", " << flaps.at(3) << std::endl;
	    }

	    std::sort(std::begin(flaps), std::end(flaps));

	    tempJetPt -= (flaps.at(0) + flaps.at(1) + flaps.at(2));
	    if(tempJetPt < 0) tempJetPt = 0;
	    if(tempJetPt > 0){
	      //	      std::cout << " Flaps: " << flaps.at(0) << ", " << flaps.at(1) << ", " << flaps.at(2) << ", " << flaps.at(3) << std::endl;
	    }
	  }

	  if(newTowIEt[etaI][phiI] == 510 || newTowIEt[etaI][phiI] == 509 || newTowIEt[etaI][phiI] == 511) tempJetPt = 2047;

	  jetPt.push_back(tempJetPt);
	  jetPhi.push_back(gtPhi(phiI));
	  jetEta.push_back(gtEta(etaI));
	  jetPhiPrev.push_back(phiI);
	  jetEtaPrev.push_back(etaI);
	  seedPt.push_back((newTowIEt[etaI])[phiI]);
	  fullSeedPt.push_back(tempSeedPt);
	  fullSeedEta.push_back(tempSeedEta);
	  fullSeedPhi.push_back(tempSeedPhi);
	  tempSeedPt.clear();
	  tempSeedEta.clear();
	  tempSeedPhi.clear();
	}
      }
    }

    if(jetPt.size() == 0) continue;

    for(unsigned int jI = 0; jI < jetPt.size() - 1; ++jI){
      for(unsigned int jI2 = jI+1; jI2 < jetPt.size(); ++jI2){
	if(jetPt.at(jI) < jetPt.at(jI2)){
	  int tempPt = jetPt.at(jI);
	  int tempPhi = jetPhi.at(jI);
	  int tempEta = jetEta.at(jI);
	  int tempPhiPrev = jetPhiPrev.at(jI);
	  int tempEtaPrev = jetEtaPrev.at(jI);
	  int tempSeed = seedPt.at(jI);
	  std::vector<int> tempSeedPt = fullSeedPt.at(jI);
	  std::vector<int> tempSeedPhi = fullSeedPhi.at(jI);
	  std::vector<int> tempSeedEta = fullSeedEta.at(jI);

	  jetPt.at(jI) = jetPt.at(jI2);
	  jetPhi.at(jI) = jetPhi.at(jI2);
	  jetEta.at(jI) = jetEta.at(jI2);
	  jetPhiPrev.at(jI) = jetPhiPrev.at(jI2);
	  jetEtaPrev.at(jI) = jetEtaPrev.at(jI2);
	  seedPt.at(jI) = seedPt.at(jI2);
	  fullSeedPt.at(jI) = fullSeedPt.at(jI2);
	  fullSeedPhi.at(jI) = fullSeedPhi.at(jI2);
	  fullSeedEta.at(jI) = fullSeedEta.at(jI2);

	  jetPt.at(jI2) = tempPt;
	  jetPhi.at(jI2) = tempPhi;
	  jetEta.at(jI2) = tempEta;
	  jetPhiPrev.at(jI2) = tempPhiPrev;
	  jetEtaPrev.at(jI2) = tempEtaPrev;
	  seedPt.at(jI2) = tempSeed;
	  fullSeedPt.at(jI2) = tempSeedPt;
	  fullSeedPhi.at(jI2) = tempSeedPhi;
	  fullSeedEta.at(jI2) = tempSeedEta;
	} 
      }
    }

    for(unsigned int jI = 0; jI < upgrade->jetIEt.size(); ++jI){
      bool isJtFound = false;

      for(unsigned int jI2 = 0; jI2 < jetPt.size(); ++jI2){
	if(jetPt.at(jI2) != upgrade->jetIEt.at(jI)) continue;
	if(jetPhi.at(jI2) != upgrade->jetIPhi.at(jI)) continue;
	if(jetEta.at(jI2) != upgrade->jetIEta.at(jI)) continue;

	isJtFound = true;
	break;
      }
     
      if(!isJtFound){
	std::cout << "Warning in event: " << evt->event << ", entry " << entry << std::endl;
	std::cout << " Cannot find jet (pt, phi, eta): " << upgrade->jetIEt.at(jI) << ", " << upgrade->jetIPhi.at(jI) << ", " << upgrade->jetIEta.at(jI) << std::endl;
	std::cout << " Candidates were: " << std::endl;
	for(unsigned int jI2 = 0; jI2 < jetPt.size(); ++jI2){
	  std::cout << "  " << jetPt.at(jI2) << ", " << jetPhi.at(jI2) << ", " << jetEta.at(jI2) << std::endl;
	}
      }
    }
  
    if(leadJtIEt != jetPt.at(0)){
      missCount++;
      //      std::cout << "Lead iet, iphi, ieta, seedEt: " << leadJtIEt << ", " << leadJtIPhi << ", " << leadJtIEta << ", " << leadJtSeedPt << std::endl;
      //      std::cout << "NJets: " << jetPt.size() << std::endl;
      for(int i = 0; i < TMath::Min(8, (int)jetPt.size()); ++i){
	//	std::cout << " " << i << "/" << jetPt.size() << ", pt, phi, eta, seed): " << jetPt.at(i) << ", " << jetPhi.at(i) << ", " << jetEta.at(i) << ", " << seedPt.at(i) << " (" << jetPhiPrev.at(i) << ", " << jetEtaPrev.at(i) << ")" << std::endl;
		
	for(unsigned int j = 0; j < fullSeedPt.at(i).size(); ++j){
	  //	  if(fullSeedPt.at(i).at(j) > 0) std::cout << "  " << fullSeedPt.at(i).at(j) << "," << fullSeedPhi.at(i).at(j) << "," << fullSeedEta.at(i).at(j) << std::endl;
	}
	
      }
    }
    else foundCount++;

    if(leadJtIEt - jetPt.at(0) != 0 || jetPhi.at(0) != leadJtIPhi || jetEta.at(0) != leadJtIEta){
      etaMiss_p->Fill(jetEta.at(0));
      phiMiss_p->Fill(jetPhi.at(0));

      etaMiss_CMSSWJet_p->Fill(leadJtIEta);
      phiMiss_CMSSWJet_p->Fill(leadJtIPhi);
    }

    for(std::map<Int_t,Int_t>::iterator it=max.begin(); it!=max.end(); ++it){
      if(it->second < tempMax[it->first]) max[it->first] = tempMax[it->first];
    }        

    fullSeedPt.clear();
  }

  /*
  std::cout << "Tower display: " << std::endl;
  
  for(int etaI = 0; etaI < nIEta; ++etaI){
    std::string disp = std::to_string(TMath::Abs(etaI-41)) + ",";
    if(disp.size() == 2) disp = " " + disp;
    std::cout << disp;
  }
  std::cout << std::endl;
  
  for(int phiI = 0; phiI < nIPhi; ++phiI){	
    for(int etaI = 0; etaI < nIEta; ++etaI){
      std::string disp = std::to_string(maxTowIEt[etaI][phiI]);
      if(disp.size() == 1) disp = " " + disp + ",";
      else if(disp.size() == 2) disp = disp + ",";
      //	else if(disp.size() == 3) disp =  disp + ",";
      std::cout << disp;
    }
    std::cout << std::endl;
  }
  */

  std::cout << "End tower display" << std::endl;

  std::cout << "Found: " << foundCount << std::endl;
  std::cout << "Miss: " << missCount << std::endl;

  std::cout << "Processing complete!" << std::endl;
  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();
  etaMiss_p->Write("", TObject::kOverwrite);
  phiMiss_p->Write("", TObject::kOverwrite);

  etaMiss_CMSSWJet_p->Write("", TObject::kOverwrite);
  phiMiss_CMSSWJet_p->Write("", TObject::kOverwrite);
  
  delete etaMiss_p;
  delete phiMiss_p;

  delete etaMiss_CMSSWJet_p;
  delete phiMiss_CMSSWJet_p;

  outFile_p->Close();
  delete outFile_p;;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage ./bin/l1OfflineSubtract.exe <inFileName> <subTypeStr>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += l1OfflineSubtract(argv[1], argv[2]);
  return retVal;
}

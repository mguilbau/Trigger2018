#include <string>
#include <iostream>
#include <vector>
#include <map>
//#include <pair>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TMath.h"

#include "include/doGlobalDebug.h"
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



const float kGTEtaLSB = 0.0435;
const float kGTPhiLSB = 0.0435;
const int kHFBegin=29;
const int kHFEnd=41;
const int kNPhi=72;


std::pair<float,float> towerEtaBounds(int ieta)
{
  if(ieta==0) ieta = 1;
  if(ieta>kHFEnd) ieta = kHFEnd;
  if(ieta<(-1*kHFEnd)) ieta = -1*kHFEnd;
  const float towerEtas[42] = {0,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.870,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172,2.322,2.5,2.650,2.853,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.191,5.191};
  return std::make_pair( towerEtas[abs(ieta)-1],towerEtas[abs(ieta)] );
}

float towerEta(int ieta)
{
  std::pair<float,float> bounds = towerEtaBounds(ieta);
  float eta = (bounds.second+bounds.first)/2.;
  float sign = ieta>0 ? 1. : -1.;
  return sign*eta;
}

float towerPhiSize()
{
  return 2.*M_PI/kNPhi;
}


float towerPhi(int iphi)
{
  float phi = (float(iphi)-0.5)*towerPhiSize();
  if (phi > M_PI) phi = phi - (2*M_PI);
  return phi;
}

int gtEta(int ieta) {

  double eta = towerEta(ieta);
  return round ( eta / kGTEtaLSB );

}

int gtPhi(int iphi) {

  double phi = towerPhi(iphi);
  if (phi<0) phi = phi + 2*M_PI;
  return round ( phi / kGTPhiLSB );

}

int l1OfflineSubtract(const std::string inFileName)
{
  const int minSeedThresh = 8;
  const int nIEta = 82;
  const int nIPhi = 72;

  TFile* outFile_p = new TFile("output/l1OfflineSubtract.root", "RECREATE");
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

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  inL1AlgoUpgradeTree_p->SetBranchAddress("L1Upgrade", &upgrade);


  inL1AlgoEvtTree_p->SetBranchStatus("*", 0);
  inL1AlgoEvtTree_p->SetBranchStatus("Event", 1);
  inL1AlgoEvtTree_p->SetBranchStatus("run", 1);
  inL1AlgoEvtTree_p->SetBranchStatus("lumi", 1);
  inL1AlgoEvtTree_p->SetBranchStatus("event", 1);

  inL1AlgoEvtTree_p->SetBranchAddress("Event", &evt);

  const Int_t nEntries = TMath::Min(200000, (Int_t)l1CaloTree_p->GetEntries());

  std::map<Int_t, Int_t> max;
  for(Int_t i = 1; i <= 41; ++i){
    max[i] = 0;
    max[-i] = 0;
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << "Processing " << nEntries << "..." << std::endl;
  for(Int_t entry = 190000; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;

    l1CaloTree_p->GetEntry(entry);
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    inL1AlgoUpgradeTree_p->GetEntry(entry);
    inL1AlgoEvtTree_p->GetEntry(entry);

    
    //    if(entry != 190099) continue;
    if(evt->run != 1) continue;
    if(evt->lumi != 2541) continue;
    if(evt->event < 254000) continue;
    if(evt->event > 254100) continue;
    if(evt->event != 254046) continue;

    std::cout << "Event: " << evt->event << std::endl;
    std::cout << "Entry: " << entry << std::endl;

    int leadJtIEt = -999;
    int leadJtIEta = -999;
    int leadJtIPhi = -999;
    
    for(unsigned int jI = 0; jI < upgrade->jetIEt.size(); ++jI){
      if(upgrade->jetIEt.at(jI) > leadJtIEt){
	leadJtIEt = upgrade->jetIEt.at(jI);
	leadJtIEta = upgrade->jetIEta.at(jI);
	leadJtIPhi = upgrade->jetIPhi.at(jI);
      }
    }

    int newTowIEt[nIEta][nIPhi];
    for(int etaI = 0; etaI < nIEta; ++etaI){
      for(int phiI = 0; phiI < nIPhi; ++phiI){
	newTowIEt[etaI][phiI] = 0;
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

      uePerEta[towers_->ieta[i]] += towers_->iet[i];
      tempMax[towers_->ieta[i]] += 1;
    }

    Float_t tempVal = (uePerEta[1] + uePerEta[-1])/2;
    uePerEta[1] = tempVal;
    uePerEta[-1] = tempVal;

    for(unsigned int i = 2; i <= 41; ++i){
      if(i%2 == 1) continue;
      tempVal = (uePerEta[i] + uePerEta[i+1])/2;
      uePerEta[i] = tempVal;
      uePerEta[i+1] = tempVal;

      tempVal = (uePerEta[-i] + uePerEta[-i-1])/2;
      uePerEta[-i] = tempVal;
      uePerEta[-i-1] = tempVal;
    }

    for(unsigned int i = 0; i < towers_->iet.size(); ++i){
      int iEt = TMath::Max(0, towers_->iet.at(i) - ((int)(uePerEta[towers_->ieta[i]]/72)));
      int ietaPos = towers_->ieta[i] + 41;
      if(towers_->ieta[i] > 0) ietaPos -= 1;
      newTowIEt[ietaPos][towers_->iphi[i]] = iEt;
    }

    std::vector<int> jetPt;
    std::vector<int> jetPhi;
    std::vector<int> jetEta;
    std::vector<int> jetPhiPrev;
    std::vector<int> jetEtaPrev;
    std::vector<int> seedPt;
    std::vector< std::vector<int> > fullSeedPt;
    std::vector< std::vector<int> > fullSeedEta;
    std::vector< std::vector<int> > fullSeedPhi;

    for(int etaI = 4; etaI < nIEta-4; ++etaI){
      int etaPos = etaI - nIEta/2;
      if(etaPos >= 0) etaPos += 1;
      
      if(TMath::Abs(etaPos) >= 25) continue;

      for(int phiI = 0; phiI < nIPhi; ++phiI){
	if(newTowIEt[etaI][phiI] >= minSeedThresh){	  

	  int tempJetPt = 0;

	  std::vector<int> tempSeedPt;
	  std::vector<int> tempSeedPhi;
	  std::vector<int> tempSeedEta;

	  bool goodSeed = true;
	  for(Int_t etaI2 = etaI - 4; etaI2 <= etaI + 4; ++etaI2){	
	    for(Int_t phiI2 = phiI - 4; phiI2 <= phiI + 4; ++phiI2){ 	      
	      int phiPos = phiI2;
	      if(phiPos < 0) phiPos += nIPhi;
	      if(phiPos >= nIPhi) phiPos -= nIPhi;

	      int etaPos2 = etaI2 - nIEta/2;
	      if(etaPos2 >= 0) etaPos2 += 1;

	      if(TMath::Abs(etaPos2) >= 25) continue;

	      tempJetPt += newTowIEt[etaI2][phiPos];

	      if(phiI == 71 && etaPos == -23) std::cout << " " << newTowIEt[etaI2][phiPos] << ", " << etaI2 << ", " << phiPos << ", (" << etaPos2 << ", " << phiI2 << ")" << std::endl;   
	      tempSeedPt.push_back(newTowIEt[etaI2][phiPos]);
	      tempSeedPhi.push_back(phiPos);
	      tempSeedEta.push_back(etaPos2);

	      if(newTowIEt[etaI2][phiPos] > 4){
		std::cout << "  " << newTowIEt[etaI2][phiPos] << ", " << etaI2 << ", " << phiPos << ", (" << etaPos2 << ", " << phiI2 << ")" << std::endl;
	      }
	      
	      if(etaI2 == 0 && phiI2 == 0) continue;
	      else{
		int dPhi = phiI2 - phiI + 4;
		int dEta = etaI2 - etaI + 4;

		if(mask_[dPhi][dEta] == 1 && newTowIEt[etaI][phiI] < newTowIEt[etaI2][phiPos]) goodSeed = false;
		else if(mask_[dPhi][dEta] == 2 && newTowIEt[etaI][phiI] <= newTowIEt[etaI2][phiPos]) goodSeed = false;

		if(!goodSeed) break;
	      }
	    }
	    if(!goodSeed) break;
	  }

	  if(!goodSeed) continue;
	
	  if(newTowIEt[etaI][phiI] == 510 || newTowIEt[etaI][phiI] == 509 || newTowIEt[etaI][phiI] == 511) tempJetPt = 2047;

	  jetPt.push_back(tempJetPt);
	  jetPhi.push_back(gtPhi(phiI));
	  jetEta.push_back(gtEta(etaPos));
	  jetPhiPrev.push_back(phiI);
	  jetEtaPrev.push_back(etaPos);
	  seedPt.push_back(newTowIEt[etaI][phiI]);
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

    std::cout << "Lead: " << leadJtIEt << ", " << leadJtIPhi << ", " << leadJtIEta << std::endl;
    std::cout << "NJets: " << jetPt.size() << std::endl;
    for(int i = 0; i < TMath::Min(4, (int)jetPt.size()); ++i){
      std::cout << " " << i << "/" << jetPt.size() << ": " << jetPt.at(i) << ", " << jetPhi.at(i) << ", " << jetEta.at(i) << ", " << seedPt.at(i) << " (" << jetPhiPrev.at(i) << ", " << jetEtaPrev.at(i) << ")" << std::endl;

      
      for(unsigned int j = 0; j < fullSeedPt.at(i).size(); ++j){
	if(fullSeedPt.at(i).at(j) > 0) std::cout << "  " << fullSeedPt.at(i).at(j) << "," << fullSeedPhi.at(i).at(j) << "," << fullSeedEta.at(i).at(j) << std::endl;
      }
      
    }
 
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

  std::cout << "Processing complete!" << std::endl;
  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  /*
  std::cout << "Maxes..." << std::endl;
  for(std::map<Int_t,Int_t>::iterator it=max.begin(); it!=max.end(); ++it){
    std::cout << " " << it->first << ", " << it->second << std::endl;
  }
  */

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
  if(argc != 2){
    std::cout << "Usage ./bin/l1OfflineSubtract.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += l1OfflineSubtract(argv[1]);
  return retVal;
}

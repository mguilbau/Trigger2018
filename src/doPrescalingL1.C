//cpp dependencies
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TH1F.h"
#include "TDatime.h"

//Local dependencies
#include "include/checkMakeDir.h"
#include "include/listOfPrimes.h"
#include "include/returnRootFileContentsList.h"
#include "include/stringUtil.h"

int doPrescalingL1(const std::string inFileName, const std::string prescaleConfigName, const std::string l1XmlFileName)
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

  if(!checkFile(l1XmlFileName)){
    std::cout << "Warning: Given l1XmlFileName \'" << l1XmlFileName << "\' is not a valid file. return 1" << std::endl;
    return 1;
  }
  else if(l1XmlFileName.find(".xml") == std::string::npos){
    std::cout << "Warning: Given l1XmlFileName \'" << l1XmlFileName << "\' is not a valid .xml file. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  std::vector<int> colRateTable;
  colRateTable.push_back(10);
  colRateTable.push_back(20);
  colRateTable.push_back(30);
  colRateTable.push_back(40);
  colRateTable.push_back(50);
  colRateTable.push_back(60);

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
    if(branchName.find("HLT_") != std::string::npos) continue;
    if(branchName.find("L1_") == std::string::npos) continue;

    finalListOfBranches.push_back(branchName);
  }

  //for(unsigned int i = 0; i< finalListOfBranches.size(); ++i)
  //   std::cout << i << " ~~> " << finalListOfBranches[i] << std::endl;

  std::vector<std::string> vZB;
  std::vector<std::string> vMB;
  std::vector<std::string> vMBnot;
  std::vector<std::string> vMu;
  std::vector<std::string> vJet;
  std::vector<std::string> vEG;
  std::vector<std::string> vETT;

  //Extract prescales
  std::ifstream prescaleConfig(prescaleConfigName.c_str());

  std::vector<std::string> vl1bitname;
  std::string l1bitname;
  std::vector< std::vector<int> > vl1ps(colRateTable.size(), std::vector<int>());
  int l1ps;

  while(!prescaleConfig.eof()){
    prescaleConfig >> l1bitname;
    vl1bitname.push_back(l1bitname);

    for(unsigned int ips = 0; ips < colRateTable.size(); ++ips)
    {
       prescaleConfig >> l1ps;
       vl1ps[ips].push_back(getNearestPrime(l1ps));
    }
  }
  vl1bitname.pop_back();
  for(unsigned int ips = 0; ips < colRateTable.size(); ++ips) vl1ps[ips].pop_back();

  for(unsigned int i = 0; i< vl1bitname.size(); ++i){
  //   std::cout << i << " ## " << vl1bitname[i];
  //   for(unsigned int j = 0; j<colRateTable.size(); ++j)
  //      std::cout << vl1ps[j][i] << " "; 
  //   std::cout << std::endl;
       if( vl1bitname[i].find("ZeroBias") != std::string::npos ){
         bool decision = false;
         for(unsigned int j = 0; j < colRateTable.size(); ++j){
            if(vl1ps[j][i] > 0) decision = true;
         }
         if(decision) vZB.push_back(vl1bitname[i]);
       }
       else if( vl1bitname[i].find("NotMinimumBias") != std::string::npos ){
         bool decision = false;
         for(unsigned int j = 0; j < colRateTable.size(); ++j){
            if(vl1ps[j][i] > 0) decision = true;
         } 
         if(decision) vMBnot.push_back(vl1bitname[i]);
       }
       else if( vl1bitname[i].find("MinimumBias") != std::string::npos && (vl1bitname[i].find("ETT") == std::string::npos && vl1bitname[i].find("Mu") == std::string::npos)){
         bool decision = false;
         for(unsigned int j = 0; j < colRateTable.size(); ++j){
            if(vl1ps[j][i] > 0) decision = true;
         } 
         if(decision) vMB.push_back(vl1bitname[i]);
       }
       else if( vl1bitname[i].find("SingleMu") != std::string::npos || vl1bitname[i].find("DoubleMu") != std::string::npos ){
         bool decision = false;
         for(unsigned int j = 0; j < colRateTable.size(); ++j){
            if(vl1ps[j][i] > 0) decision = true;
         } 
         if(decision) vMu.push_back(vl1bitname[i]);
       }
       else if( vl1bitname[i].find("SingleJet") != std::string::npos || vl1bitname[i].find("DoubleJet") != std::string::npos ){
         bool decision = false;
         for(unsigned int j = 0; j < colRateTable.size(); ++j){
            if(vl1ps[j][i] > 0) decision = true;
         } 
         if(decision) vJet.push_back(vl1bitname[i]);
       }
       else if( vl1bitname[i].find("SingleEG") != std::string::npos || vl1bitname[i].find("DoubleEG") != std::string::npos ){
         bool decision = false;
         for(unsigned int j = 0; j < colRateTable.size(); ++j){
            if(vl1ps[j][i] > 0) decision = true;
         } 
         if(decision) vEG.push_back(vl1bitname[i]);
       }
       else if( vl1bitname[i].find("ETT") != std::string::npos ){
         bool decision = false;
         for(unsigned int j = 0; j < colRateTable.size(); ++j){
            if(vl1ps[j][i] > 0) decision = true;
         } 
         if(decision) vETT.push_back(vl1bitname[i]);
       }
  }

  //Display for check
  std::cout << "~~~> ZeroBias" << std::endl;
  for(unsigned int itrg=0; itrg < vZB.size(); ++itrg)
  {
     std::cout << "    -" << vZB[itrg].c_str() << std::endl;
  }
  std::cout << "~~~> MinimumBias NOT" << std::endl;
  for(unsigned int itrg=0; itrg < vMBnot.size(); ++itrg)
  {
     std::cout << "    -" << vMBnot[itrg].c_str() << std::endl;
  }
  std::cout << "~~~> MinimumBias" << std::endl;
  for(unsigned int itrg=0; itrg < vMB.size(); ++itrg)
  {
     std::cout << "    -" << vMB[itrg].c_str() << std::endl;
  }
  std::cout << "~~~> sMu + dMu" << std::endl;
  for(unsigned int itrg=0; itrg < vMu.size(); ++itrg)
  {
     std::cout << "    -" << vMu[itrg].c_str() << std::endl;
  }
  std::cout << "~~~> sJet + dJet" << std::endl;
  for(unsigned int itrg=0; itrg < vJet.size(); ++itrg)
  {
     std::cout << "    -" << vJet[itrg].c_str() << std::endl;
  }
  std::cout << "~~~> sEG + dEG" << std::endl;
  for(unsigned int itrg=0; itrg < vEG.size(); ++itrg)
  {
     std::cout << "    -" << vEG[itrg].c_str() << std::endl;
  }
  std::cout << "~~~> ETT + ETTAsym + ETTNot" << std::endl;
  for(unsigned int itrg=0; itrg < vETT.size(); ++itrg)
  {
     std::cout << "    -" << vETT[itrg].c_str() << std::endl;
  }

  std::vector<std::string>::iterator it;
  //Size checking
  if( finalListOfBranches.size() != vl1bitname.size()){
     std::cerr << "Size inconsistancy between HLT Tree and L1 Menu" << std::endl
               << "   - HLT Tree number of branches: " << finalListOfBranches.size() << std::endl
               << "   - L1 Menu number of bits: " << vl1bitname.size() << std::endl;

     for(unsigned int i = 0; i < vl1bitname.size(); ++i){
        std::cout << "Checking item: " << i << " -- ";
        it = std::find(finalListOfBranches.begin(), finalListOfBranches.end(), vl1bitname[i]);
        if(it == finalListOfBranches.end()){
          std::cerr << "Cannot find matching between HLT Tree and L1 Menu for bit: " << vl1bitname[i].c_str() << std::endl;
        }
        else{
          std::cout << *it << " | " << vl1bitname[i] << std::endl;
        }
     }

     return 0;
  }
  else{
     std::cout << "HLT Tree vs. L1bits size check passed" << std::endl;
  }

  //Branch matching
  for(unsigned int i = 0; i < vl1bitname.size(); ++i){
     //std::cout << "Checking item: " << i << " -- ";
     it = std::find(finalListOfBranches.begin(), finalListOfBranches.end(), vl1bitname[i]);
     if(it == finalListOfBranches.end()){
       std::cerr << "Cannot find matching between HLT Tree and L1 Menu for bit: " << *it << std::endl;
       return 0;
     }
     //std::cout << std::endl;
  }
  std::cout << "HLT Tree vs. L1bits branch matching passed" << std::endl;
   

  const Int_t nTrig = vl1bitname.size();
  const Int_t nPS = colRateTable.size();
  Int_t trigVal[nTrig] = {0};
  Int_t trigFires[nTrig] = {0};
  Int_t trigPrescaledFires[nTrig][nPS] = {{0}};
  Int_t totalFires[nPS] = {0};

  Int_t trigFires_ZB[nPS]    = {0};  
  Int_t trigFires_MBnot[nPS] = {0}; 
  Int_t trigFires_MB[nPS]    = {0}; 
  Int_t trigFires_Mu[nPS]    = {0}; 
  Int_t trigFires_Jet[nPS]   = {0}; 
  Int_t trigFires_EG[nPS]    = {0}; 
  Int_t trigFires_ETT[nPS]   = {0}; 

  hltTree_p->SetBranchStatus("*", 0);
  for(Int_t bI = 0; bI < nTrig; ++bI){
    hltTree_p->SetBranchStatus(vl1bitname.at(bI).c_str(), 1);
    hltTree_p->SetBranchAddress(vl1bitname.at(bI).c_str(), &(trigVal[bI]));
  }

  const Int_t nEntries = hltTree_p->GetEntries();
  const Int_t nDiv = TMath::Max(1, nEntries/20);

  std::cout << "Processing " << nEntries << "..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;

    hltTree_p->GetEntry(entry);

    std::vector<bool> doesGlobalFire(colRateTable.size(),false);
    std::vector<bool> doesZBFire(colRateTable.size(),false);
    std::vector<bool> doesMBnotFire(colRateTable.size(),false);
    std::vector<bool> doesMBFire(colRateTable.size(),false);
    std::vector<bool> doesMuFire(colRateTable.size(),false);
    std::vector<bool> doesJetFire(colRateTable.size(),false);
    std::vector<bool> doesEGFire(colRateTable.size(),false);
    std::vector<bool> doesETTFire(colRateTable.size(),false);

    for(Int_t bI = 0; bI < nTrig; ++bI){
      if(trigVal[bI] == 1){
        for(unsigned int ips = 0; ips < colRateTable.size(); ++ips){
	   if(trigFires[bI]%vl1ps.at(ips).at(bI) == 0 && vl1ps.at(ips).at(bI)>0){
	     doesGlobalFire.at(ips) = true;
	     ++trigPrescaledFires[bI][ips];

             for(unsigned int ipd = 0; ipd < vZB.size(); ++ipd)
             {
                if( vl1bitname[bI] == vZB[ipd] ) doesZBFire[ips] = true;
             }
             for(unsigned int ipd = 0; ipd < vMBnot.size(); ++ipd)
             {
                if( vl1bitname[bI] == vMBnot[ipd] ) doesMBnotFire[ips] = true;
             }
             for(unsigned int ipd = 0; ipd < vMB.size(); ++ipd)
             {
                if( vl1bitname[bI] == vMB[ipd] ) doesMBFire[ips] = true;
             }
             for(unsigned int ipd = 0; ipd < vMu.size(); ++ipd)
             {
                if( vl1bitname[bI] == vMu[ipd] ) doesMuFire[ips] = true;
             }
             for(unsigned int ipd = 0; ipd < vJet.size(); ++ipd)
             {
                if( vl1bitname[bI] == vJet[ipd] ) doesJetFire[ips] = true;
             }
             for(unsigned int ipd = 0; ipd < vEG.size(); ++ipd)
             {
                if( vl1bitname[bI] == vEG[ipd] ) doesEGFire[ips] = true;
             }
             for(unsigned int ipd = 0; ipd < vETT.size(); ++ipd)
             {
                if( vl1bitname[bI] == vETT[ipd] ) doesETTFire[ips] = true;
             }
	   }
        }
	++trigFires[bI];
      }
    }
    for(unsigned int ips = 0; ips < colRateTable.size(); ++ips){
       if(doesGlobalFire.at(ips)) ++totalFires[ips];
       if(doesZBFire.at(ips)) ++trigFires_ZB[ips];
       if(doesMBnotFire.at(ips)) ++trigFires_MBnot[ips];
       if(doesMBFire.at(ips)) ++trigFires_MB[ips];
       if(doesMuFire.at(ips)) ++trigFires_Mu[ips];
       if(doesJetFire.at(ips)) ++trigFires_Jet[ips];
       if(doesEGFire.at(ips)) ++trigFires_EG[ips];
       if(doesETTFire.at(ips)) ++trigFires_ETT[ips];
    }
  }
    
  inFile_p->Close();
  delete inFile_p;


  std::vector< std::vector<double> > trgRateBPS;
  std::vector< std::vector<double> > trgRateBPS_err;
  trgRateBPS.resize(nTrig);
  trgRateBPS_err.resize(nTrig);

  for(Int_t bI = 0; bI < nTrig; ++bI){
    for(Int_t ips = 0; ips < nPS; ++ips){
      trgRateBPS[bI].push_back(((Double_t)trigFires[bI])*colRateTable.at(ips)*1000./((Double_t)nEntries));
      trgRateBPS_err[bI].push_back(((Double_t)trigFires[bI] + TMath::Sqrt(trigFires[bI]))*colRateTable.at(ips)*1000./((Double_t)nEntries) - trgRateBPS[bI][ips]);
    }
  }

  std::vector< std::vector<double> > trgRateAPS;
  std::vector< std::vector<double> > trgRateAPS_err;
  trgRateAPS.resize(nTrig);
  trgRateAPS_err.resize(nTrig);

  for(Int_t bI = 0; bI < nTrig; ++bI){
    for(Int_t ips = 0; ips < nPS; ++ips){
      trgRateAPS[bI].push_back(((Double_t)trigPrescaledFires[bI][ips])*colRateTable.at(ips)*1000./((Double_t)nEntries));
      trgRateAPS_err[bI].push_back(((Double_t)trigPrescaledFires[bI][ips] + TMath::Sqrt(trigPrescaledFires[bI][ips]))*colRateTable.at(ips)*1000./((Double_t)nEntries) - trgRateAPS[bI][ips]);
    }
  }

  std::cout.precision(5);
  for(Int_t bI = 0; bI < nTrig; ++bI){
    std::cout << std::setw(60) << vl1bitname.at(bI).c_str(); 
    for(Int_t ips = 0; ips < nPS; ++ips){
       if(vl1ps.at(ips).at(bI) < 0) vl1ps[ips][bI] = 0;
       std::cout << std::setw(3) << " | " << std::setw(8) << trgRateBPS[bI][ips] << std::setw(1) << " " << std::setw(3) << vl1ps.at(ips).at(bI) << std::setw(1) << " " << std::setw(7) << trgRateAPS[bI][ips]; 
    }
    std::cout << std::endl;
  }



  for(Int_t ips = 0; ips < nPS; ++ips){
     std::cout << std::endl << std::endl;
     std::cout << "SUMMARY: Total fires, Rate at " << (Int_t)colRateTable.at(ips) << " kHz (Hz): " 
               << totalFires[ips] << ", " << ((Double_t)totalFires[ips])*colRateTable.at(ips)*1000./((Double_t)nEntries) << " Hz" << std::endl;

     std::cout << std::setw(21) << "            - ZB: " 
               << std::setw(21) << trigFires_ZB[ips] << std::setw(2) << ", " << std::setw(8) << ((Double_t)trigFires_ZB[ips])*colRateTable.at(ips)*1000./((Double_t)nEntries) << std::setw(3) << " Hz" << std::endl;
     std::cout << std::setw(21) << "            - MBnot: " 
               << std::setw(21) << trigFires_MBnot[ips] << std::setw(2) << ", " << std::setw(8) << ((Double_t)trigFires_MBnot[ips])*colRateTable.at(ips)*1000./((Double_t)nEntries) << std::setw(3) << " Hz" << std::endl;
     std::cout << std::setw(21) << "            - MB: " 
               << std::setw(21) << trigFires_MB[ips] << std::setw(2) << ", " << std::setw(8) << ((Double_t)trigFires_MB[ips])*colRateTable.at(ips)*1000./((Double_t)nEntries) << std::setw(3) << " Hz" << std::endl;
     std::cout << std::setw(21) << "            - Mu: " 
               << std::setw(21) << trigFires_Mu[ips] << std::setw(2) << ", " << std::setw(8) << ((Double_t)trigFires_Mu[ips])*colRateTable.at(ips)*1000./((Double_t)nEntries) << std::setw(3) << " Hz" << std::endl;
     std::cout << std::setw(21) << "            - Jet: " 
               << std::setw(21) << trigFires_EG[ips] << std::setw(2) << ", " << std::setw(8) << ((Double_t)trigFires_Jet[ips])*colRateTable.at(ips)*1000./((Double_t)nEntries) << std::setw(3) << " Hz" << std::endl;
     std::cout << std::setw(21) << "            - EG: " 
               << std::setw(21) << trigFires_EG[ips] << std::setw(2) << ", " << std::setw(8) << ((Double_t)trigFires_EG[ips])*colRateTable.at(ips)*1000./((Double_t)nEntries) << std::setw(3) << " Hz" << std::endl;
     std::cout << std::setw(21) << "            - ETT: " 
               << std::setw(21) << trigFires_ETT[ips] << std::setw(2) << ", " << std::setw(8) << ((Double_t)trigFires_ETT[ips])*colRateTable.at(ips)*1000./((Double_t)nEntries) << std::setw(3) << " Hz" << std::endl;

  }
////
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 3){
    std::cout << "Usage: ./bin/doPrescalingL1.exe <inFileName> <prescaleConfigName> <l1XmlFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += doPrescalingL1(argv[1], argv[2], argv[3]);
  return retVal;
}

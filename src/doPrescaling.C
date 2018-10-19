//cpp dependencies
#include <iostream>
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

int doPrescaling(const std::string inFileName, const std::string prescaleConfigName, const std::string l1XmlFileName, const double inCollisionRatekHz)
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
  std::vector<int> trigThreshold;

  std::map<std::string, bool> pdMapIsFired;
  std::map<std::string, int> pdMapToFires;

  std::map<std::string, bool> subPDMapIsFired;
  std::map<std::string, int> subPDMapToFires;

  std::vector<std::string> uniqueSubPD;
  std::vector<std::vector<double> > subPDThresholds;

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
    tempStr.replace(0, tempStr.find(",")+1, "");
    int threshold = std::stoi(tempStr.substr(0, tempStr.find(",")));

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
    trigThreshold.push_back(threshold);

    if(pdMapToFires.count(pd) == 0) pdMapToFires[pd] = 0;
    if(subPDMapToFires.count(subPD) == 0) subPDMapToFires[subPD] = 0;

    
    int idPos = -1;
    bool isUniqueSubPD = true;
    for(unsigned int pI = 0; pI < uniqueSubPD.size(); ++pI){
      if(isStrSame(uniqueSubPD.at(pI), subPD)){
	isUniqueSubPD = false;
	idPos = pI;
	break;
      }
    }

    if(isUniqueSubPD){
      uniqueSubPD.push_back(subPD);
      subPDThresholds.push_back({(Double_t)threshold});
    }
    else{
      subPDThresholds.at(idPos).push_back((Double_t)threshold);
    }
  }
  prescaleConfig.close();

  const Int_t nTrig = trigNames.size();
  Int_t trigVal[nTrig];
  Int_t trigFires[nTrig];
  Int_t trigPrescaledFires[nTrig];

  std::string matchingL1FromXML[nTrig];

  std::ifstream xmlFile(l1XmlFileName.c_str());
  std::vector<std::string> fullXmlList;
  std::vector<bool> fullXmlListMatched;
  while(std::getline(xmlFile, tempStr)){
    if(tempStr.find("<name>L1_") == std::string::npos) continue;
    tempStr.replace(tempStr.find("<name>"), std::string("<name>").size(), "");
    tempStr.replace(tempStr.find("</name>"), std::string("</name>").size(), "");

    while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, "");}

    fullXmlList.push_back(tempStr);
    fullXmlListMatched.push_back(false);
  }
  xmlFile.close();

  for(Int_t tI = 0; tI < nTrig; ++tI){
    std::string tempTrigName = trigNames.at(tI);
    tempTrigName.replace(0,4,"");
    tempTrigName.replace(tempTrigName.rfind("_v"), tempTrigName.size(), "");

    bool trigIsFound = false;
    for(unsigned int xI = 0; xI < fullXmlList.size(); ++xI){
      if(fullXmlListMatched.at(xI)) continue;
      std::string tempL1Name = fullXmlList.at(xI);
      while(tempL1Name.find("_") != std::string::npos){tempL1Name.replace(tempL1Name.find("_"), 1, "");}
      //      std::cout << "tempL1Name: \'" << tempL1Name << "\', \'" << tempTrigName << "\'" << std::endl;

      if(isStrSame(tempL1Name, tempTrigName)){
	fullXmlListMatched.at(xI) = true;
	matchingL1FromXML[tI] = fullXmlList.at(xI);
	trigIsFound = true;
	break;
      }
    }
    
    if(!trigIsFound){
      std::cout << "WARNING: \'" << trigNames.at(tI) << "\' is not found in xml \'" << l1XmlFileName << "\'. Setting blank" << std::endl;
      matchingL1FromXML[tI] = "MISSING";
    }
  }

  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){
    outFileName.replace(0, outFileName.find("/")+1, "");
  }
  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), outFileName.size(), "");

  outFileName = "output/" + dateStr + "/" + outFileName + "_RatePlots_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  const Int_t nSubPD = uniqueSubPD.size();
  TH1F* subPDFiresAtKHz_h[nSubPD];
  for(Int_t sI = 0; sI < nSubPD; ++sI){
    std::vector<double> temp = subPDThresholds.at(sI);
    std::sort(std::begin(temp), std::end(temp));
    subPDThresholds.at(sI) = temp;

    const Int_t nBins = subPDThresholds.at(sI).size();
    Double_t bins[nBins+1];

    if(nBins != 1){
      bins[0] = subPDThresholds.at(sI).at(0) - (subPDThresholds.at(sI).at(1) - subPDThresholds.at(sI).at(0))/2.;
      
      for(Int_t bIX = 0; bIX < nBins-1; ++bIX){
	bins[bIX+1] = subPDThresholds.at(sI).at(bIX) + (subPDThresholds.at(sI).at(bIX+1) - subPDThresholds.at(sI).at(bIX))/2.;
      }
      
      bins[nBins] = subPDThresholds.at(sI).at(nBins-1) + (subPDThresholds.at(sI).at(nBins-1) - subPDThresholds.at(sI).at(nBins-2))/2.;
    }
    else{
      std::cout << "Warning: subPD \'" << uniqueSubPD.at(sI) << "\' has only 1 threshold" << std::endl;
      bins[0] = subPDThresholds.at(sI).at(0) - 1;
      bins[1] = subPDThresholds.at(sI).at(0) + 1;
    }
      
    subPDFiresAtKHz_h[sI] = new TH1F(("subPDFires" + std::to_string((Int_t)inCollisionRatekHz) + "KHz_" + uniqueSubPD.at(sI) + "_h").c_str(), (";Threshold;Rate at " + std::to_string((Int_t)inCollisionRatekHz) + " kHz (Hz)").c_str(), nBins, bins);
  }

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

  outFile_p->cd();

  for(Int_t bI = 0; bI < nTrig; ++bI){
    Int_t subPDPos = -1;

    for(Int_t sI = 0; sI < nSubPD; ++sI){
      if(isStrSame(uniqueSubPD.at(sI), trigSubPD[sI])){
	subPDPos = sI;
	break;
      }
    }

    Int_t binPos = subPDFiresAtKHz_h[subPDPos]->FindBin(trigThreshold.at(bI));

    Double_t val = ((Double_t)trigFires[bI])*inCollisionRatekHz*1000./((Double_t)nEntries);
    Double_t valErr = ((Double_t)trigFires[bI] + TMath::Sqrt(trigFires[bI]))*inCollisionRatekHz*1000./((Double_t)nEntries) - val;

    subPDFiresAtKHz_h[subPDPos]->SetBinContent(binPos, val);
    subPDFiresAtKHz_h[subPDPos]->SetBinError(binPos, valErr);
  }
  
  for(Int_t sI = 0; sI < nSubPD; ++sI){
    subPDFiresAtKHz_h[sI]->Write("", TObject::kOverwrite);
    delete subPDFiresAtKHz_h[sI];
  }

  outFile_p->Close();
  delete outFile_p;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  const std::string outPrescaleName = "output/" + dateStr + "/prescales_" + std::to_string((Int_t)inCollisionRatekHz) + "kHz_" + dateStr + ".csv";
  std::ofstream outPrescaleFile(outPrescaleName.c_str());
  outPrescaleFile << "L1 Trigger,Unprescaled Rate in Hz (" << std::to_string((Int_t)inCollisionRatekHz) << " kHz),Prescale (" << std::to_string((Int_t)inCollisionRatekHz) << " kHz),Prescaled Rate in Hz (" << std::to_string((Int_t)inCollisionRatekHz) << " kHz)," << std::endl;
  for(Int_t bI = 0; bI < nTrig; ++bI){
    outPrescaleFile << matchingL1FromXML[bI] << "," << ((Double_t)trigFires[bI])*inCollisionRatekHz*1000./((Double_t)nEntries) << "," << trigPrescale[bI] << "," << ((Double_t)trigPrescaledFires[bI])*inCollisionRatekHz*1000./((Double_t)nEntries) << "," << std::endl;
  }
  outPrescaleFile.close();
  
  std::cout << "#: Name, PD, SubPD, Final prescale, Fires, Prescaled Fires, Rate at " << (Int_t)inCollisionRatekHz << "kHz (Hz), Prescaled Rate at " << (Int_t)inCollisionRatekHz << "kHz (Hz)" << std::endl;
  for(Int_t bI = 0; bI < nTrig; ++bI){
    std::cout << " " << bI << "/" << nTrig << ": " << trigNames[bI] << ", " << trigPD[bI] << ", " << trigSubPD[bI] << ", " << trigPrescale[bI] << ", " << trigFires[bI] << ", " << trigPrescaledFires[bI] << ", " << ((Double_t)trigFires[bI])*inCollisionRatekHz*1000./((Double_t)nEntries) << ", " << ((Double_t)trigPrescaledFires[bI])*inCollisionRatekHz*1000./((Double_t)nEntries) << std::endl;
  }

  std::cout << std::endl;

  std::cout << "PD, Total Prescaled Fires, Rate at " << (Int_t)inCollisionRatekHz << "kHz (Hz)" << std::endl;
  for(auto const &iter : pdMapToFires){
    std::cout << iter.first << ", " << iter.second << ", " << ((Double_t)iter.second)*inCollisionRatekHz*1000./((Double_t)nEntries) << std::endl;
  }

  std::cout << std::endl;

  std::cout << "SubPD, Total Prescaled Fires, Rate at " << (Int_t)inCollisionRatekHz << "kHz (Hz)" << std::endl;
  for(auto const &iter : subPDMapToFires){
    std::cout << iter.first << ", " << iter.second << ", " << ((Double_t)iter.second)*inCollisionRatekHz*1000./((Double_t)nEntries) << std::endl;
  }

  std::cout << std::endl;

  std::cout << "Total fires, Rate at " << (Int_t)inCollisionRatekHz << " kHz (Hz): " << totalFires << ", " << ((Double_t)totalFires)*inCollisionRatekHz*1000./((Double_t)nEntries) << std::endl;

  std::cout << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 5){
    std::cout << "Usage: ./bin/doPrescaling.exe <inFileName> <prescaleConfigName> <l1XmlFileName> <inCollisionRatekHz>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += doPrescaling(argv[1], argv[2], argv[3], std::stod(argv[4]));
  return retVal;
}

//cpp dependencies
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"

//local dependencies
#include "include/mntToXRootdFileString.h"
#include "include/plotUtilities.h"

std::string to_string_with_precision(double a_value, const int n)
{
  std::ostringstream out;
  out << std::setprecision(n) << a_value;
  return out.str();
}

int quickHiBin(const std::string inFileName, const Int_t nBins, const std::string tagStr = "")
{
  TFile* inFile_p = TFile::Open(mntToXRootdFileString(inFileName).c_str(), "READ");
  TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
  Float_t hiHF_;

  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("hiHF", 1);

  hiTree_p->SetBranchAddress("hiHF", &hiHF_);

  std::vector<double> hiHFVect;

  for(Int_t entry = 0; entry < hiTree_p->GetEntries(); ++entry){
    hiTree_p->GetEntry(entry);

    hiHFVect.push_back(((double)hiHF_));
  }  

  inFile_p->Close();
  delete inFile_p;

  std::sort(std::begin(hiHFVect), std::end(hiHFVect));

  std::cout << "hiHF min, max: " << hiHFVect.at(0) << ", " << hiHFVect.at(hiHFVect.size()-1) << std::endl;

  std::cout << "In percentiles of 5" << std::endl;

  const int precisionDec = 6;
  Double_t bins[nBins+1];
  Int_t precision[nBins+1];

  for(int i = 0; i < nBins; ++i){
    Int_t posLow = i*hiHFVect.size()/nBins;
    Int_t posHi = (i+1)*hiHFVect.size()/nBins;
    if(i == nBins-1) posHi -= 1;

    double valLow = hiHFVect.at(posLow);
    double valHi = hiHFVect.at(posHi);
    int div = 10;
    precision[i] = precisionDec+1;
    precision[i+1] = precisionDec+1;
    while(((int)valLow/div) > 0){
      precision[i] += 1;
      div *= 10;
    }
    div = 10;
    while(((int)valHi/div) > 0){
      precision[i+1] += 1;
      div *= 10;
    }

    std::cout << prettyString(100. - (i+1)*100./(Double_t)nBins, 1, false) << "-" << prettyString(100. - i*100./(Double_t)nBins, 1, false) << "%: " << to_string_with_precision(valLow, precision[i]) << "-" << to_string_with_precision(valHi, precision[i+1]) << std::endl;
    bins[i] = hiHFVect.at(posLow);
  }

  bins[0] = 0.0;
  bins[nBins] = 50000.0;

  precision[0] = precisionDec+1;
  precision[nBins] = precisionDec+1;
  int div = 10;
  while(((int)bins[0]/div) > 0){
    precision[0] += 1;
    div *= 10;
  }
  div = 10;
  while(((int)bins[nBins]/div) > 0){
    precision[nBins] += 1;
    div *= 10;
  }

  //Create a quick header
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  const std::string dateStr2 = std::to_string(date->GetYear()) + "." + std::to_string(date->GetMonth()) + "." + std::to_string(date->GetDay()) ;
  delete date;

  std::string finalTagStr = "";
  std::string finalTagStr2 = "";
  if(tagStr.size() != 0){
    finalTagStr = tagStr + "_";
    finalTagStr2 = "_" + tagStr;
  }

  const std::string outFileName = "quickCentralityTable_nBins" + std::to_string(nBins) + "_" + finalTagStr + dateStr + ".h";
  std::ofstream outFile(outFileName.c_str());

  outFile << "#ifndef QUICKCENTRALITYTABLE_" << finalTagStr << dateStr << "_H" << std::endl;
  outFile << "#define QUICKCENTRALITYTABLE_" << finalTagStr << dateStr << "_H" << std::endl;

  outFile << std::endl;
  outFile << "//Created with file \'" << inFileName << "\'" << std::endl;
  outFile << "//File tagStr give \'" << tagStr << "\'" << std::endl;
  outFile << "//Created on " << dateStr2 << std::endl;
  outFile << "//Number of events=" << hiHFVect.size() << std::endl;
  outFile << "//nBins=" << nBins << ", meaning percentiles of " << prettyString(100./(Double_t)nBins,2,false) << "%" << std::endl;
  outFile << "//Reworked to map on typical hiBin distribution of 200" << std::endl;
  outFile << "//hiBin interval of " << prettyString(200./(Double_t)nBins, 2, false) << std::endl;
  outFile << std::endl;

  outFile << "const Int_t nBins" << finalTagStr2 << " = " << nBins << ";" << std::endl;
  outFile << "const Double_t bins" << finalTagStr2 << "[nBins" << finalTagStr2 << "+1] = {";
  for(Int_t i = 0; i < nBins; ++i){
    outFile << to_string_with_precision(bins[i], precision[i]) << ", ";
  }
  outFile << to_string_with_precision(bins[nBins], precision[nBins]) << "};" << std::endl;;

  outFile << std::endl;

  outFile << "//Returns lower bound value of hibin at the granularity of nBins=" << nBins << std::endl;
  outFile << "Int_t getHiBinFromHiHF" << finalTagStr2 << "(const Double_t hiHF)" << std::endl;
  outFile << "{" << std::endl;
  outFile << "  Int_t binPos = -1;" << std::endl;
  outFile << "  for(int i = 0; i < nBins" << finalTagStr2 << "; ++i){" << std::endl;
  outFile << "    if(hiHF >= bins" << finalTagStr2 << "[i] && hiHF < bins" << finalTagStr2 << "[i+1]){" << std::endl;
  outFile << "      binPos = i;" << std::endl;
  outFile << "      break;" << std::endl;
  outFile << "    }" << std::endl;
  outFile << "  }" << std::endl;
  outFile << std::endl;
  outFile << "  binPos = nBins" << finalTagStr2 << " - 1 - binPos;" << std::endl;
  outFile <<   std::endl;  
  outFile << "  return (Int_t)(200*((Double_t)binPos)/((Double_t)nBins" << finalTagStr2 << "));" << std::endl;
  outFile << "}" << std::endl;

  outFile << std::endl;
  outFile << "#endif" << std::endl;

  outFile.close();


  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 3 || argc > 4){
    std::cout << "Usage: ./bin/quickHiBin.exe <inFileName> <nBins> <tagStr-Opt>" << std::endl;
    return 1;     
  }
  
  int retVal = 0;
  if(argc == 3) retVal += quickHiBin(argv[1], std::stoi(argv[2]));
  else if(argc == 4) retVal += quickHiBin(argv[1], std::stoi(argv[2]), argv[3]);
  return retVal;
}

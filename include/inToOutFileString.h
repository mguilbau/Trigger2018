#ifndef INTOOUTFILESTRING_H
#define INTOOUTFILESTRING_H

#include <string>

#include "TDatime.h"

std::string inToOutFileString(std::string inFileName, std::string appStr = "")
{
  TDatime* date = new TDatime();
  while(inFileName.find("/") != std::string::npos){inFileName.replace(0, inFileName.find("/")+1, "");}

  const Int_t nExt = 2;
  const std::string ext[nExt] = {".root", ".txt"};
  for(Int_t i = 0; i < nExt; ++i){
    if(inFileName.find(ext[i]) != std::string::npos) inFileName.replace(inFileName.find(ext[i]), ext[i].size(), "");
  }
  
  if(appStr.size() != 0) inFileName = inFileName + "_" + appStr;
  inFileName = inFileName + "_" + std::to_string(date->GetDate()) + ".root";

  delete date;

  return inFileName;
}

#endif


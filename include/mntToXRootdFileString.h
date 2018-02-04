#ifndef MNTTOXROOTDFILESTRING_H
#define MNTTOXROOTDFILESTRING_H

#include <string>

std::string mntToXRootdFileString(std::string inFileName)
{
  const std::string mntStr = "/mnt/hadoop/cms";
  if(inFileName.size() > mntStr.size()){
    if(inFileName.substr(0, mntStr.size()).find(mntStr) != std::string::npos){
      inFileName.replace(0, mntStr.size(), "");
      inFileName = "root://xrootd.cmsaf.mit.edu/" + inFileName;
    }
  }

  return inFileName;
}

#endif

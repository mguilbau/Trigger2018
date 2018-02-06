#ifndef TEXTFILETOVECTOR_H
#define TEXTFILETOVECTOR_H

#include <fstream>
#include <vector>
#include <string>

std::vector<std::string> textFileToVector(const std::string inFileName, const std::string filterStr)
{
  std::vector<std::string> retVect;

  std::ifstream file(inFileName.c_str());
  std::string tempStr;
  while(std::getline(file, tempStr)){
    if(tempStr.size() == 0) continue;
    if(tempStr.find(filterStr) == std::string::npos) continue;

    retVect.push_back(tempStr);
  }

  file.close();

  return retVect;
}

bool textFileToVectorAppend(std::vector<std::string>* inVect_p, const std::string inFileName, const std::string filterStr)
{
  const unsigned int startSize = inVect_p->size();
  std::vector<std::string> tempVect = textFileToVector(inFileName, filterStr);
  inVect_p->insert(inVect_p->end(), tempVect.begin(), tempVect.end());
  return (inVect_p->size() != startSize);
}

#endif

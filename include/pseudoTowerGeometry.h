#ifndef PSEUDOTOWERGEOMETRY_H
#define PSEUDOTOWERGEOMETRY_H

#include <vector>
#include "TMath.h"

class pseudoTowGeo{
 public:
  pseudoTowGeo();

  const static int nEtaTow = 82;
  const double etaTowBounds[nEtaTow+1] = {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -3, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0.000, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.870, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.740, 1.830, 1.930, 2.043, 2.172, 2.322, 2.500, 2.650, 2.853, 3.000, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};
  const double etaAbsTowBounds[nEtaTow/2+1] = {0.000, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.870, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.740, 1.830, 1.930, 2.043, 2.172, 2.322, 2.500, 2.650, 2.853, 3.000, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

  const int nTowInPhi[nEtaTow] = {18, 18, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 72, 72, 72, 72, 72, 72, 72,72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 18, 18};

  const int nAbsTowInPhi[nEtaTow/2] = {72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 18, 18};
  
  std::vector<double> etaTowBoundsVect;
  std::vector<int> nTowInPhiVect;
  std::vector<double> etaAbsTowBoundsVect;
  std::vector<int> nAbsTowInPhiVect;

  std::vector<double> getEtaTowBounds();
  std::vector<int> getNTowInPhi();

  std::vector<double> getAbsEtaTowBounds();

  std::vector<double> getPhiBoundsForEta(const double inEta);
};

pseudoTowGeo::pseudoTowGeo()
{
  etaTowBoundsVect.reserve(nEtaTow+1);
  for(unsigned int iter = 0; iter < nEtaTow+1; ++iter){
    etaTowBoundsVect.push_back(etaTowBounds[iter]);
  }

  etaAbsTowBoundsVect.reserve(nEtaTow/2+1);
  for(unsigned int iter = 0; iter < nEtaTow/2+1; ++iter){
    etaAbsTowBoundsVect.push_back(etaAbsTowBounds[iter]);
  }

  nTowInPhiVect.reserve(nEtaTow);
  for(unsigned int iter = 0; iter < nEtaTow; ++iter){
    nTowInPhiVect.push_back(nTowInPhi[iter]);
  }

  nTowInPhiVect.reserve(nEtaTow/2);
  for(unsigned int iter = 0; iter < nEtaTow/2; ++iter){
    nAbsTowInPhiVect.push_back(nAbsTowInPhi[iter]);
  }
  return;
}

std::vector<double> pseudoTowGeo::getEtaTowBounds(){return etaTowBoundsVect;}
std::vector<int> pseudoTowGeo::getNTowInPhi(){return nTowInPhiVect;}
std::vector<double> pseudoTowGeo::getAbsEtaTowBounds(){return etaAbsTowBoundsVect;}

std::vector<double> pseudoTowGeo::getPhiBoundsForEta(const double inEta)
{
  std::vector<double> phiBounds;

  int etaPos = -1;
  for(int etaIter = 0; etaIter < nEtaTow; ++etaIter){
    if(inEta >= etaTowBounds[etaIter] && inEta < etaTowBounds[etaIter+1]){
      etaPos = etaIter;
      break;
    }
  }
  
  const int nTow = nTowInPhi[etaPos];
  phiBounds.reserve(nTow+1);
  for(int phiIter = 0; phiIter < nTow+1; ++phiIter){
    phiBounds.push_back(-TMath::Pi() + phiIter*2.*TMath::Pi()/nTow);
  }
  
  return phiBounds;
}

#endif

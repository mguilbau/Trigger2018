#ifndef L1TOOLS_H
#define L1TOOLS_H

//#include "L1AnalysisEventDataFormat.h"
//#include "L1AnalysisL1CaloTowerDataFormat.h"
//#include "L1AnalysisL1UpgradeDataFormat_OLD_AsOf20180907.h"

//const float M_PI = 3.14159265358979323846;
const float kGTEtaLSB = 0.0435;
const float kGTPhiLSB = 0.0435;
const int kHFBegin=29;
const int kHFEnd=41;
const int kNPhi=72;

//For Jet -> Tower conversions
const Int_t nIEtaJetNew = 239;
const Int_t iEtaTowerNew[nIEtaJetNew] = {-41,-999,-999,-40,-999,-999,-999,-999,-999,-39,-999,-999,-999,-38,-999,-999,-999,-37,-999,-999,-999,-36,-999,-999,-999,-35,-999,-999,-999,-34,-999,-999,-999,-33,-999,-999,-999,-32,-999,-999,-999,-31,-999,-999,-999,-30,-999,-999,-999,-999,-29,-999,-999,-999,-999,-999,-28,-999,-999,-999,-27,-999,-999,-999,-26,-999,-999,-25,-999,-999,-999,-24,-999,-23,-999,-999,-22,-999,-21,-999,-20,-999,-19,-999,-18,-999,-17,-999,-16,-999,-15,-999,-14,-999,-13,-999,-12,-999,-11,-999,-10,-999,-9,-999,-8,-999,-7,-999,-6,-999,-5,-999,-4,-999,-3,-999,-2,-999,-1,-999,1,-999,2,-999,3,-999,4,-999,5,-999,6,-999,7,-999,8,-999,9,-999,10,-999,11,-999,12,-999,13,-999,14,-999,15,-999,16,-999,17,-999,18,-999,19,-999,20,-999,21,-999,22,-999,-999,23,-999,24,-999,-999,-999,25,-999,-999,26,-999,-999,-999,27,-999,-999,-999,28,-999,-999,-999,-999,-999,29,-999,-999,-999,-999,30,-999,-999,-999,31,-999,-999,-999,32,-999,-999,-999,33,-999,-999,-999,34,-999,-999,-999,35,-999,-999,-999,36,-999,-999,-999,37,-999,-999,-999,38,-999,-999,-999,39,-999,-999,-999,-999,-999,40,-999,-999,41};
const Int_t iEtaJetNew[nIEtaJetNew] = {-119,-118,-117,-116,-115,-114,-113,-112,-111,-110,-109,-108,-107,-106,-105,-104,-103,-102,-101,-100,-99,-98,-97,-96,-95,-94,-93,-92,-91,-90,-89,-88,-87,-86,-85,-84,-83,-82,-81,-80,-79,-78,-77,-76,-75,-74,-73,-72,-71,-70,-69,-68,-67,-66,-65,-64,-63,-62,-61,-60,-59,-58,-57,-56,-55,-54,-53,-52,-51,-50,-49,-48,-47,-46,-45,-44,-43,-42,-41,-40,-39,-38,-37,-36,-35,-34,-33,-32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119};

const Int_t nIPhiJetNew = 143;
const Int_t iPhiTowerNew[nIPhiJetNew] = {1,-999,2,-999,3,-999,4,-999,5,-999,6,-999,7,-999,8,-999,9,-999,10,-999,11,-999,12,-999,13,-999,14,-999,15,-999,16,-999,17,-999,18,-999,19,-999,20,-999,21,-999,22,-999,23,-999,24,-999,25,-999,26,-999,27,-999,28,-999,29,-999,30,-999,31,-999,32,-999,33,-999,34,-999,35,-999,36,-999,37,-999,38,-999,39,-999,40,-999,41,-999,42,-999,43,-999,44,-999,45,-999,46,-999,47,-999,48,-999,49,-999,50,-999,51,-999,52,-999,53,-999,54,-999,55,-999,56,-999,57,-999,58,-999,59,-999,60,-999,61,-999,62,-999,63,-999,64,-999,65,-999,66,-999,67,-999,68,-999,69,-999,70,-999,71,-999,72};
const Int_t iPhiJetNew[nIPhiJetNew] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143};


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

int convertGtEtaToTowerIEta(int gtIEta)
{
  if(TMath::Abs(gtIEta) > iEtaJetNew[nIEtaJetNew-1]) return -999;
  return iEtaTowerNew[gtIEta-iEtaJetNew[0]];
}

int convertGtPhiToTowerIPhi(int gtIPhi)
{
  if(gtIPhi < iPhiJetNew[0]) return -999;
  else if(gtIPhi > iPhiJetNew[nIPhiJetNew-1]) return -999;
  return iPhiTowerNew[gtIPhi-iPhiJetNew[0]];
}

#endif

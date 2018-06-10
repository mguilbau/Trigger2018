#ifndef L1TOOLS_H
#define L1TOOLS_H

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

#endif

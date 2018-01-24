#ifndef RUNLUMIEVTKEY_H
#define RUNLUMIEVTKEY_H

#include <map>

class runLumiEvtKey{
 public:
  std::map<unsigned long long, int> runLumiEvtKeyToEntry;

  runLumiEvtKey(){};
  bool addKey(unsigned int run, unsigned int lumi, unsigned long evt, int entry);
  int getEntryFromKey(unsigned int run, unsigned int lumi, unsigned long evt);

 private:
  unsigned long long makeKey(unsigned int run, unsigned int lumi, unsigned long evt);

  static const unsigned long long runMult = 10000000000000;
  static const unsigned long long lumiMult = 1000000000;

  static const unsigned int runMax = 999999;
  static const unsigned int lumiMax = 9999;
  static const unsigned long evtMax = 999999999;
};


unsigned long long runLumiEvtKey::makeKey(unsigned int run, unsigned int lumi, unsigned long evt)
{
  unsigned long long outNum = ((unsigned long long)run)*runMult;
  outNum += ((unsigned long long)lumi)*lumiMult;
  outNum += evt;

  //  std::cout << run << ", " << lumi << ", " << evt << ", " << outNum << std::endl;
  
  return outNum;
}


bool runLumiEvtKey::addKey(unsigned int run, unsigned int lumi, unsigned long evt, int entry)
{
  if(run > runMax){
    std::cout << "Input run \'" << run << "\' is greater than maxVal for key, \'" << runMax << "\'. In runLumiEvtKey::addKey, return false" << std::endl;
    return false;
  }
  else if(lumi > lumiMax){
    std::cout << "Input lumi \'" << lumi << "\' is greater than maxVal for key, \'" << lumiMax << "\'. In runLumiEvtKey::addKey, return false" << std::endl;
    return false;
  }
  else if(evt > evtMax){
    std::cout << "Input evt \'" << evt << "\' is greater than maxVal for key, \'" << evtMax << "\'. In runLumiEvtKey::addKey, return false" << std::endl;
    return false;
  } 
  unsigned long long key = makeKey(run, lumi, evt);
  int keyCount = runLumiEvtKeyToEntry.count(key);
  if(keyCount != 0){
    std::cout << "Input run, lumi, evt (" << run << ", " << lumi << ", " << evt << ") gives repeat key! In runLumiEvtKey::addKey, return false" << std::endl;
    return false;
  }
  runLumiEvtKeyToEntry[key] = entry;

  return true;
}

int runLumiEvtKey::getEntryFromKey(unsigned int run, unsigned int lumi, unsigned long evt)
{
  if(run > runMax){
    std::cout << "Input run \'" << run << "\' is greater than maxVal for key, \'" << runMax << "\'. In runLumiEvtKey::getEntryFromKey, return false" << std::endl;
    return false;
  }
  else if(lumi > lumiMax){
    std::cout << "Input lumi \'" << lumi << "\' is greater than maxVal for key, \'" << lumiMax << "\'. In runLumiEvtKey::getEntryFromKey, return false" << std::endl;
    return false;
  }
  else if(evt > evtMax){
    std::cout << "Input evt \'" << evt << "\' is greater than maxVal for key, \'" << evtMax << "\'. In runLumiEvtKey::getEntryFromKey, return false" << std::endl;
    return false;
  }

  unsigned long long key = makeKey(run, lumi, evt);
  int outEntry = -1;
  if(runLumiEvtKeyToEntry.count(key) == 1) outEntry = runLumiEvtKeyToEntry[key];

  return outEntry;
}

#endif

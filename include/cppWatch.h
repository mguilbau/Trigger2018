#ifndef CPPWATCH_H
#define CPPWATCH_H

#include <ctime>

class cppWatch{
 public:
  cppWatch();
  ~cppWatch(){}

  double totalInt;
  double currentInt;
  std::clock_t c_start;
  
  void start(){c_start = std::clock(); return;};
  void stop();
  double total(){return totalInt;}
  double current(){return currentInt;}
  void clear(){totalInt = 0; currentInt = 0; return;}
};

cppWatch::cppWatch(){totalInt = 0; currentInt = 0; return;}

void cppWatch::stop()
{
  std::clock_t c_end = std::clock();
  
  totalInt += double((c_end - c_start)/(double)CLOCKS_PER_SEC); 
  currentInt = double((c_end - c_start)/(double)CLOCKS_PER_SEC);
  return;
}

#endif

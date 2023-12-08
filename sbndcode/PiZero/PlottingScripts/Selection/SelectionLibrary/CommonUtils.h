#include "TString.h"
#include "TChain.h"
#include <iostream>

namespace CommonUtils
{
  const double goalPOT     = 10e20;
  const double potPerSpill = 5e12;
  const double goalSpills  = goalPOT / potPerSpill;

  TString POTString();

  double GetPOT(TChain *subruns);

  int GetGenEvents(TChain *subruns);

  void GetScaling(TChain *rockboxSubruns, TChain *intimeSubruns, double &rockboxScaling, double &intimeScaling);
}

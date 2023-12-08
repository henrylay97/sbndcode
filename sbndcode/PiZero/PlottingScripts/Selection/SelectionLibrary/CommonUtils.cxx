#include "CommonUtils.h"

TString CommonUtils::POTString()
{
  TString potString = Form(" (%g POT)", goalPOT);
  potString.ReplaceAll("e+", "x10^{");
  potString.ReplaceAll(" POT", "} POT");

  return potString;
}

double CommonUtils::GetPOT(TChain *subruns)
{
  double sum = 0., pot = 0;

  subruns->SetBranchAddress("pot", &pot);

  for(int i = 0; i < subruns->GetEntries(); ++i)
    {
      subruns->GetEntry(i);
      sum += pot;
    }

  return sum;
}

int CommonUtils::GetGenEvents(TChain *subruns)
{
  int sum = 0, ngenevts = 0;

  subruns->SetBranchAddress("ngenevts", &ngenevts);

  for(int i = 0; i < subruns->GetEntries(); ++i)
    {
      subruns->GetEntry(i);
      sum += ngenevts;
    }

  return sum;
}

void CommonUtils::GetScaling(TChain *rockboxSubruns, TChain *intimeSubruns, double &rockboxScaling, double &intimeScaling)
{
  const double rockboxPOT = GetPOT(rockboxSubruns);
  const int rockboxSpills = GetGenEvents(rockboxSubruns);
  const int intimeSpills  = GetGenEvents(intimeSubruns);

  rockboxScaling = goalPOT / rockboxPOT;

  const double scaledRockboxSpills = rockboxScaling * rockboxSpills;
  
  intimeScaling = (goalSpills - scaledRockboxSpills) / intimeSpills;
}

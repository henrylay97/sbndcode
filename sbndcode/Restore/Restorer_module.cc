////////////////////////////////////////////////////////////////////////
// Class:       Restorer
// Plugin Type: analyzer (art v3_06_03)
// File:        Restorer_module.cc
//
// Generated at Fri Apr  9 05:26:45 2021 by Henry Lay using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace restore {
  class Restorer;
}


class restore::Restorer : public art::EDAnalyzer {
public:
  explicit Restorer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Restorer(Restorer const&) = delete;
  Restorer(Restorer&&) = delete;
  Restorer& operator=(Restorer const&) = delete;
  Restorer& operator=(Restorer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

};


restore::Restorer::Restorer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void restore::Restorer::analyze(art::Event const& e)
{
  std::cout << "Viewing event: " << e.id() << std::endl;
}

DEFINE_ART_MODULE(restore::Restorer)

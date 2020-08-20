////////////////////////////////////////////////////////////////////////
// Class:       PhotonValidation
// Module Type: analyzer
// File:        PhotonValidation_module.cc
//
// Generated at Tue Nov  5 08:54:30 2019 by Dominic Barker using artmod
// from cetpkgsupport v1_14_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"

//Shower Finder includes 
#include "larreco/RecoAlg/MCRecoUtils/RecoUtils.h"
#include "larreco/RecoAlg/MCRecoUtils/ShowerUtils.h"


//C++ Includes
#include <map>
#include <vector>
#include <iostream>

class PhotonValidation;

class PhotonValidation : public art::EDAnalyzer {
public:
  explicit PhotonValidation(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotonValidation(PhotonValidation const &) = delete;
  PhotonValidation(PhotonValidation &&) = delete;
  PhotonValidation & operator = (PhotonValidation const &) = delete;
  PhotonValidation & operator = (PhotonValidation &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  //Function to get the shower mothers
  std::map<int,std::vector<int> > GetShowerMothers(std::vector<art::Ptr<recob::Hit> >& hits, std::vector<art::Ptr<sim::SimChannel> >& simchannels);

private:

  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  float fSimEnergyCut;
  float fDensityCut;
  std::string fHitsModuleLabel;
  std::string fLArGeantModuleLabel;
  // Declare member data here.

};


PhotonValidation::PhotonValidation(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset),
  fSimEnergyCut(pset.get<float>("SimEnergyCut")), //~.1.2
  fDensityCut(pset.get<float>("DensityCut")), //~1.1
  fHitsModuleLabel(pset.get<std::string>("HitsModuleLabel")), //gaushit
  fLArGeantModuleLabel(pset.get<std::string>("LArGeantModuleLabel")) //largeant
 // More initializers here.
{}

void PhotonValidation::analyze(art::Event const & evt)
{

  //Getting the Hit Information
  art::Handle<std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if(evt.getByLabel(fHitsModuleLabel,hitListHandle))
    {art::fill_ptr_vector(hits, hitListHandle);}

  //Getting the SimWire Information
  //Get the SimChannels so that we can find the IDEs deposited on them.
  art::Handle<std::vector<sim::SimChannel> > simChannelHandle;
  std::vector<art::Ptr<sim::SimChannel> > simchannels;
  if(evt.getByLabel(fLArGeantModuleLabel,simChannelHandle))
    {art::fill_ptr_vector(simchannels, simChannelHandle);}
  

  //Save the shower chain. The int is the track id of the mother the shower. The rest of the particles in the em shower have there track ids saved in the vector.
  std::map<int,std::vector<int> > ShowerMothers = GetShowerMothers(hits,simchannels);

  std::cout << "You have " << ShowerMothers.size() << " showers left in the event" << std::endl;
  //Now you have your final showers. You might want to use some more logic now to only keep showers come from a neutral pion.


}

std::map<int,std::vector<int> > PhotonValidation::GetShowerMothers(std::vector<art::Ptr<recob::Hit> >& hits, std::vector<art::Ptr<sim::SimChannel> >& simchannels){
  
  //List the monte carlo in the particles in the event
  const sim::ParticleList& particles = particleInventory->ParticleList();

  //Loop over the monte carlo particles and save them in something more understandable
  std::map<int,const simb::MCParticle*> trueParticles;

  //Save the shower chain. The int is the track id of the mother the shower. The rest of the particles in the em shower have there track ids saved in the vector.
  std::map<int,std::vector<int> > ShowersMothers; 

  //This function returns the energy deposited in MeV in the liquid argon. The energy has to bee seen somehow on the wires. The int corresponds to the track id and the flot to the energy that was deposited.  
  std::map<int,float> MCTrack_Energy_map = RecoUtils::TrueEnergyDepositedFromMCTracks(simchannels);

  //This function returns the number of wires per plane that each track id deposited.
  std::map<int,std::map<geo::PlaneID,int> > MCWires_Track_Map = RecoUtils::NumberofMCWiresHitMap(simchannels);

    //Make a Loop over the mc particles and add to the map.
    for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
      const simb::MCParticle *particle = particleIt->second;
      trueParticles[particle->TrackId()] = particle;

      //example output statment to show how to use MCParticles. Look at the doxygen.
      std::cout << "True Particle with track ID: " << particle->TrackId() << " Has code of: " << particle->PdgCode() << " and Energy of: " << particle->E() << " With Mother: " << particle->Mother() << " Proccess: " << particle->Process() << " End Process: "  << particle->EndProcess() << std::endl;
      
    }

    //Get the shower mothers. The key of the map is the mother and the vector is the daughter chain. 
    ShowersMothers = ShowerUtils::GetShowerMothersCandidates(trueParticles);
    
    //Cut on the shower mothers with energy. If the total true energy of particle is less than the fSimEnergyCut the candidate is removed. 
    ShowerUtils::CutShowerMothersByE(ShowersMothers, trueParticles, fSimEnergyCut);

    //Cut on shower mothers with density. If number of hits / the number of wires hit > fDensity cut the candidate is removed 
    ShowerUtils::CutShowerMothersByDensity(ShowersMothers, trueParticles, hits, fDensityCut);
    
    //If more than 90% of the energy of the particle is deposited in the detector the shower is kept.
    ShowerUtils::RemoveNoneContainedParticles(ShowersMothers,trueParticles,MCTrack_Energy_map);

    return ShowersMothers;
  
}

DEFINE_ART_MODULE(PhotonValidation)
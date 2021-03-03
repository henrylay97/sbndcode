////////////////////////////////////////////////////////////////////////
// Class:       vertexEVD
// Plugin Type: analyzer (art v3_05_01)
// File:        vertexEVD_module.cc
//
// Generated at Tue Sep 29 05:57:52 2020 by Henry Lay using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

//Standard
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Sim Base
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimChannel.h"

//Reco Base
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

//Tools
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

//LArSoft
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larcorealg/Geometry/ChannelMapStandardAlg.h"
#include "larcore/Geometry/Geometry.h"

//Root
#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "TVector3.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TList.h"

class vertexEVD;


class vertexEVD : public art::EDAnalyzer {
public:
  explicit vertexEVD(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  vertexEVD(vertexEVD const&) = delete;
  vertexEVD(vertexEVD&&) = delete;
  vertexEVD& operator=(vertexEVD const&) = delete;
  vertexEVD& operator=(vertexEVD&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;


private:

  // Declare member data here.
  void ClearData();
  void SetupMaps(art::Event const &e);

  void SimGraphs(art::Event const &e);
  void RecoGraphs(art::Event const &e);

  const art::Ptr<simb::MCParticle> GetMCParticle(int const &trackID);
  const art::Ptr<recob::PFParticle> GetPrimaryPFP();
  const art::Ptr<recob::PFParticle> GetPFP(unsigned int const &index);
  const std::vector<art::Ptr<recob::PFParticle> > GetPrimaryNeutrinoPFPs();
  const std::vector<art::Ptr<recob::PFParticle> > GetSlicePrimary(std::vector<art::Ptr<recob::PFParticle> > const &slicePFPs);
  const art::Ptr<simb::MCParticle> GetTrueParticle(std::vector<art::Ptr<recob::Hit> > const &trackHits);
  const art::Ptr<simb::MCParticle> GetTrueShowerParticle(std::vector<art::Ptr<recob::Hit> > const &showerHits);
  unsigned int HierarchyPrimary(art::Ptr<recob::PFParticle> const &pfp);
  float TrueTrackLength(art::Ptr<simb::MCParticle> const &particle);


  art::Handle<std::vector<simb::MCTruth> > eHandleNeutrinos;
  art::Handle<std::vector<simb::MCTruth> > eHandleCosmics;
  art::Handle<std::vector<simb::MCParticle> > eHandleParticles;
  art::Handle<std::vector<recob::Hit> > eHandleHits;
  art::Handle<std::vector<recob::Track> > eHandleTracks;
  art::Handle<std::vector<recob::Shower> > eHandleShowers;
  art::Handle<std::vector<recob::PFParticle> > eHandlePFPs;
  art::Handle<std::vector<sim::SimChannel> > eHandleSimChannels;
  art::Handle<std::vector<recob::Vertex> > eHandleVertices;
  art::Handle<std::vector<recob::Slice> > eHandleSlices;
  art::Handle<std::vector<recob::SpacePoint> > eHandleSpacePoints;

  detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
  detinfo::DetectorPropertiesData propData = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob();
  art::ServiceHandle<geo::Geometry> GeoHandler;
  geo::GeometryCore const* geom = art::ServiceHandle<geo::Geometry>()->provider();

  std::string fNuGenModuleLabel, fCosmicGenModuleLabel, fLArGeantModuleLabel,
    fPFParticleModuleLabel, fTrackModuleLabel, fShowerModuleLabel,
    fVertexModuleLabel, fHitsModuleLabel, fClusterModuleLabel,
    fSpacePointModuleLabel, fParticleIDModuleLabel, fSliceModuleLabel;
  bool fProcessCosmics;

  TTree *fEventTree;

  int eventID, subRunID, runID;

  TList *mcGraphs = new TList();
  TList *trackGraphs = new TList();
  TList *showerGraphs = new TList();
  TList *mcVertices = new TList();
  TList *recoVertices = new TList();

  std::vector<int> mcPDGList, trackPDGList, showerPDGList;
  std::vector<float> vertexErrorList, xCorrection;

  std::map<int,int> MCPDGMap;
  TVector3 true_vertex;
  bool negTPC;
  
  
};

vertexEVD::vertexEVD(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fNuGenModuleLabel (p.get<std::string>("NuGenModuleLabel")),
  fCosmicGenModuleLabel (p.get<std::string>("CosmicGenModuleLabel")),
  fLArGeantModuleLabel (p.get<std::string>("LArGeantModuleLabel")),
  fPFParticleModuleLabel (p.get<std::string>("PFParticleModuleLabel")),
  fTrackModuleLabel (p.get<std::string>("TrackModuleLabel")),
  fShowerModuleLabel (p.get<std::string>("ShowerModuleLabel")),
  fVertexModuleLabel (p.get<std::string>("VertexModuleLabel")),
  fHitsModuleLabel (p.get<std::string>("HitsModuleLabel")),
  fClusterModuleLabel (p.get<std::string>("ClusterModuleLabel")),
  fSpacePointModuleLabel (p.get<std::string>("SpacePointModuleLabel")),
  fParticleIDModuleLabel (p.get<std::string>("ParticleIDModuleLabel")),
  fSliceModuleLabel (p.get<std::string>("SliceModuleLabel")),
  fProcessCosmics (p.get<bool>("ProcessCosmics"))
  // More initializers here.
  {
    // Call appropriate consumes<>() for any products to be retrieved by this module.
    art::ServiceHandle<art::TFileService> tfs;

    mcGraphs->SetName("mcGraphs");
    trackGraphs->SetName("trackGraphs");
    showerGraphs->SetName("showerGraphs");
    mcVertices->SetName("mcVertices");
    recoVertices->SetName("recoVertices");
    
    fEventTree = tfs->make<TTree>("EventTree","Data tree");
    fEventTree->Branch("eventID",&eventID);
    fEventTree->Branch("subRunID",&subRunID);
    fEventTree->Branch("runID",&runID);
    fEventTree->Branch("mcGraphs",&mcGraphs);
    fEventTree->Branch("mcPDGList",&mcPDGList);
    fEventTree->Branch("mcVertices",&mcVertices);
    fEventTree->Branch("trackGraphs",&trackGraphs);
    fEventTree->Branch("trackPDGList",&trackPDGList);
    fEventTree->Branch("showerGraphs",&showerGraphs);
    fEventTree->Branch("showerPDGList",&showerPDGList);
    fEventTree->Branch("recoVertices",&recoVertices);    
    fEventTree->Branch("vertexErrorList",&vertexErrorList);    
  }

void vertexEVD::SimGraphs(art::Event const &e)
{
  TGraph2D *vertices = new TGraph2D();

  for(unsigned int nu_i = 0; nu_i < eHandleNeutrinos->size(); ++nu_i){
    const art::Ptr<simb::MCTruth> truthNeutrino(eHandleNeutrinos,nu_i);
    if(truthNeutrino.isNull()) continue;

    if(truthNeutrino->Origin() == 1) {
      const simb::MCNeutrino neutrino = truthNeutrino->GetNeutrino();
      const simb::MCParticle neutrinoParticle = neutrino.Nu();

      if(neutrinoParticle.Vx() < 0) negTPC = true;

      art::FindManyP<simb::MCParticle> nuParticleAssn(eHandleNeutrinos,e,fLArGeantModuleLabel);
      std::vector<art::Ptr<simb::MCParticle> > particles = nuParticleAssn.at(truthNeutrino.key());

      for(unsigned int part_i = 0; part_i < particles.size(); ++part_i) {
	const art::Ptr<simb::MCParticle> particle = particles[part_i];
	if(particle->Mother() != 0) continue;
	if(particle->StatusCode() != 1) continue;
	if(particle->PdgCode() == 2112) continue;
	if(particle->NumberTrajectoryPoints() < 3) continue;

	TGraph2D *g = new TGraph2D();
	
	for(unsigned int i = 0; i < particle->NumberTrajectoryPoints(); ++i){
	  if(particle->Vx(i) == -999 || particle->Vy(i) == -999 || particle->Vz(i) == -999) break;
	  if(particle->Vx(i) > 200 || particle->Vx(i) < -200 || particle->Vy(i) > 200 || particle->Vy(i) < -200 ||
	     particle->Vz(i) > 500 || particle->Vz(i) < 0) break;
	  g->SetPoint(i,particle->Vx(i),particle->Vy(i),particle->Vz(i));
	}
	mcPDGList.push_back(particle->PdgCode());
	mcGraphs->Add(g);
      }
      detinfo::DetectorPropertiesData propD = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
      util::GeometryUtilities geomU = util::GeometryUtilities(*geom,clockData,propD);
      
      for(auto pID : geom->IteratePlaneIDs()){
        xCorrection.push_back(propD.ConvertTicksToX(clockData.TPCG4Time2Tick(neutrinoParticle.T()),pID) - propD.ConvertTicksToX(0,pID));
	break;
      }
      
      vertices->SetPoint(nu_i,neutrinoParticle.Vx(),neutrinoParticle.Vy(),neutrinoParticle.Vz());
      true_vertex = TVector3(neutrinoParticle.Vx(),neutrinoParticle.Vy(),neutrinoParticle.Vz());
    }
  }
  mcVertices->Add(vertices);
}

void vertexEVD::RecoGraphs(art::Event const &e)
{
  art::FindManyP<recob::Vertex> pfpVertexAssoc(eHandlePFPs,e,fVertexModuleLabel);
  art::FindManyP<recob::Track> pfpTrackAssoc(eHandlePFPs,e,fTrackModuleLabel);
  art::FindManyP<recob::Shower> pfpShowerAssoc(eHandlePFPs,e,fShowerModuleLabel);
  art::FindManyP<recob::PFParticle> slicePFPAssoc(eHandleSlices,e,fSliceModuleLabel);
  art::FindManyP<recob::SpacePoint> pfpSPAssoc(eHandlePFPs,e,fPFParticleModuleLabel);
  art::FindManyP<recob::Hit> trackHitsAssoc(eHandleTracks,e,fTrackModuleLabel);
  art::FindManyP<recob::Hit> showerHitsAssoc(eHandleShowers,e,fShowerModuleLabel);

  TGraph2D *vertices = new TGraph2D();
  
  for(unsigned int slice_i = 0; slice_i < eHandleSlices->size(); ++slice_i){
    const art::Ptr<recob::Slice> slice(eHandleSlices,slice_i);
    std::vector<art::Ptr<recob::PFParticle> > slicePFPs = slicePFPAssoc.at(slice.key());

    std::vector<art::Ptr<recob::PFParticle> > primaryPFPs = GetSlicePrimary(slicePFPs);
    if(primaryPFPs.size() == 0) continue;
    else if(primaryPFPs.size() > 1) {std::cout << "Multiple primaries!!!" << std::endl; continue;}

    art::Ptr<recob::PFParticle> nuPrimary = primaryPFPs[0];

    std::vector<art::Ptr<recob::Vertex> > pfpVertices = pfpVertexAssoc.at(nuPrimary.key());
    if(pfpVertices.size() != 1) {std::cout << "Multiple vertices!!!" << std::endl; return;}

    float x = pfpVertices[0]->position().X();
    if(negTPC) x -= xCorrection[0];
    else x += xCorrection[0];
        
    vertices->SetPoint(slice_i,x,pfpVertices[0]->position().Y(),pfpVertices[0]->position().Z());
    TVector3 reco_vertex = TVector3(x,pfpVertices[0]->position().Y(),pfpVertices[0]->position().Z());
    
    vertexErrorList.push_back((reco_vertex - true_vertex).Mag());

    std::vector<long unsigned int> nuDaughters = nuPrimary->Daughters();

    for(auto id : nuDaughters){
      const art::Ptr<recob::PFParticle> daughter = GetPFP(id);
      std::vector<art::Ptr<recob::SpacePoint> > spacePoints = pfpSPAssoc.at(daughter.key());
      if(spacePoints.size() < 3) continue;

      TGraph2D *g = new TGraph2D(spacePoints.size());
      
      for(unsigned int i = 0; i < spacePoints.size(); ++i){
	art::Ptr<recob::SpacePoint> sp = spacePoints[i];
	if(sp->position().X() == -999 || sp->position().Y() == -999 || sp->position().Z() == -999) break;
	float corrX = sp->position().X();
	if(corrX < 0) corrX -= xCorrection[0];
	else corrX += xCorrection[0];

	g->SetPoint(i,corrX,sp->position().Y(),sp->position().Z());
      }
      
      if(daughter->PdgCode() == 13) {
	std::vector<art::Ptr<recob::Track> > daughterTracks = pfpTrackAssoc.at(daughter.key());
        if(daughterTracks.size() != 1) continue;
        else{
	  art::Ptr<recob::Track> track = daughterTracks[0];
	  std::vector<art::Ptr<recob::Hit> > trackHits = trackHitsAssoc.at(track.key());
	  art::Ptr<simb::MCParticle> mc = GetTrueParticle(trackHits);
	  if(mc.isNull()) return;
	  
	  trackPDGList.push_back(mc->PdgCode());
	  trackGraphs->Add(g);
	}
      }
      else if(daughter->PdgCode() == 11) {
	std::vector<art::Ptr<recob::Shower> > daughterShowers = pfpShowerAssoc.at(daughter.key());
	if(daughterShowers.size() != 1) continue;
	else {
	  art::Ptr<recob::Shower> shower = daughterShowers[0];
	  std::vector<art::Ptr<recob::Hit> > showerHits = showerHitsAssoc.at(shower.key());
	  art::Ptr<simb::MCParticle> mc = GetTrueShowerParticle(showerHits);
	  if(mc.isNull()) return;
	  
	  showerPDGList.push_back(mc->PdgCode());
	  showerGraphs->Add(g);
	}
      }
    }
  }
  recoVertices->Add(vertices);
}

void vertexEVD::SetupMaps(art::Event const &e) 
{
  MCPDGMap.clear();

  for(unsigned int part_i = 0; part_i < eHandleParticles->size(); ++part_i) {
    const art::Ptr<simb::MCParticle> particle(eHandleParticles,part_i);
    MCPDGMap[particle->TrackId()] = particle->PdgCode();
  }

}

const art::Ptr<simb::MCParticle> vertexEVD::GetMCParticle(int const &trackID)
{
  for(unsigned int part_i = 0; part_i < eHandleParticles->size(); ++part_i) {
    const art::Ptr<simb::MCParticle> particle(eHandleParticles,part_i);
    if(particle->TrackId() == trackID) return particle;
  }

  const art::Ptr<simb::MCParticle> nullReturn;
  return nullReturn;
}

const art::Ptr<recob::PFParticle> vertexEVD::GetPrimaryPFP() 
{
  for(unsigned int pfp_i = 0; pfp_i < eHandlePFPs->size(); ++pfp_i) {
    const art::Ptr<recob::PFParticle> pfp(eHandlePFPs,pfp_i);
    if(pfp->IsPrimary()) return pfp;
  }
  
  art::Ptr<recob::PFParticle> nullReturn;
  return nullReturn;
}

const art::Ptr<recob::PFParticle> vertexEVD::GetPFP(unsigned int const &index) 
{
  for(unsigned int pfp_i = 0; pfp_i < eHandlePFPs->size(); ++pfp_i) {
    const art::Ptr<recob::PFParticle> pfp(eHandlePFPs,pfp_i);
    if(pfp->Self()==index) return pfp;
  }
  
  art::Ptr<recob::PFParticle> nullReturn;
  return nullReturn;
}

const std::vector<art::Ptr<recob::PFParticle> > vertexEVD::GetPrimaryNeutrinoPFPs() 
{
  std::vector<art::Ptr<recob::PFParticle> > primaries;
  for(unsigned int pfp_i = 0; pfp_i < eHandlePFPs->size(); ++pfp_i) {
    const art::Ptr<recob::PFParticle> pfp(eHandlePFPs,pfp_i);
    if(pfp->IsPrimary()){
      if(pfp->PdgCode() == 12 || pfp->PdgCode() == -12 
	 || pfp->PdgCode() == 14 || pfp->PdgCode() == -14){
	primaries.push_back(pfp);
      }
    }
  }
  return primaries;
}

const std::vector<art::Ptr<recob::PFParticle> > vertexEVD::GetSlicePrimary(std::vector<art::Ptr<recob::PFParticle> > const &slicePFPs)
{
  std::vector<art::Ptr<recob::PFParticle> > primaries;
  for(auto pfp : slicePFPs) {
    if(pfp->IsPrimary()){
      primaries.push_back(pfp);
    }
  }
  return primaries;
}

const art::Ptr<simb::MCParticle> vertexEVD::GetTrueParticle(std::vector<art::Ptr<recob::Hit> > const &trackHits)
{
  int particleID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,trackHits,true);

  for(unsigned int i = 0; i < eHandleParticles->size(); i++) {
      const::art::Ptr<simb::MCParticle> particle(eHandleParticles,i);
      if(particle->TrackId() == particleID) return particle;
    }

  art::Ptr<simb::MCParticle> nullReturn;
  return nullReturn;
}

const art::Ptr<simb::MCParticle> vertexEVD::GetTrueShowerParticle(std::vector<art::Ptr<recob::Hit> > const &showerHits)
{
  art::Ptr<simb::MCParticle> original = GetTrueParticle(showerHits);
  if(original.isNull()) return original;

  int trackID = original->Mother();

  while(MCPDGMap[trackID] == 11 || 
	MCPDGMap[trackID] == -11 || 
	MCPDGMap[trackID] == 22) {
    original = GetMCParticle(trackID);
    trackID = original->Mother();
  }
  return original;
}

unsigned int vertexEVD::HierarchyPrimary(art::Ptr<recob::PFParticle> const &pfp)
{
  bool regress = true;
  art::Ptr<recob::PFParticle> parent = pfp;
  
  while(regress){
    if(parent->IsPrimary()) return parent->Self();
    else parent = GetPFP(parent->Parent());
  }
  return 999999;  
}

float vertexEVD::TrueTrackLength(art::Ptr<simb::MCParticle> const &particle)
{
  float length = 0;
  unsigned int nTrajPoints = particle->NumberTrajectoryPoints();

  if(nTrajPoints < 2) return length;

  for(unsigned int point = 1; point < nTrajPoints; ++point) {
    TVector3 l = particle->Position(point).Vect();
    if(l.X() > 200 || l.X() < -200 || l.Y() > 200 || l.Y() < -200 ||
       l.Z() > 500 || l.Z() < 0) break;

    TVector3 diff = particle->Position(point).Vect() - particle->Position(point-1).Vect();
    length += TMath::Sqrt(diff.Mag2());
  }

  return length;
}

void vertexEVD::ClearData()
{
  mcGraphs->Clear(); trackGraphs->Clear(); showerGraphs->Clear();
  mcPDGList.clear(); trackPDGList.clear(); showerPDGList.clear();
  mcVertices->Clear(); recoVertices->Clear();
  vertexErrorList.clear(); xCorrection.clear();
  true_vertex = TVector3(-999,-999,-999);
  negTPC = false;
}


void vertexEVD::analyze(art::Event const& e)
{
  e.getByLabel(fNuGenModuleLabel,eHandleNeutrinos);
  if(fProcessCosmics) e.getByLabel(fCosmicGenModuleLabel,eHandleCosmics);
  e.getByLabel(fLArGeantModuleLabel,eHandleParticles);
  e.getByLabel(fHitsModuleLabel,eHandleHits);
  e.getByLabel(fTrackModuleLabel,eHandleTracks);
  e.getByLabel(fShowerModuleLabel,eHandleShowers);
  e.getByLabel(fPFParticleModuleLabel,eHandlePFPs);
  e.getByLabel(fLArGeantModuleLabel,eHandleSimChannels);
  e.getByLabel(fVertexModuleLabel,eHandleVertices);
  e.getByLabel(fSliceModuleLabel,eHandleSlices);
  e.getByLabel(fSpacePointModuleLabel,eHandleSpacePoints);

  clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  eventID = e.event();
  subRunID = e.subRun();
  runID = e.run();

  ClearData();
  SimGraphs(e);
  RecoGraphs(e);
  
  fEventTree->Fill();
}

DEFINE_ART_MODULE(vertexEVD)

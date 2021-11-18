#include "CRTHitRecoAlg.h"

#include "lardataalg/DetectorInfo/DetectorClocksData.h"

namespace sbnd{

CRTHitRecoAlg::CRTHitRecoAlg(const Config& config){

  this->reconfigure(config);
}


CRTHitRecoAlg::CRTHitRecoAlg(){
}


CRTHitRecoAlg::~CRTHitRecoAlg(){

}


void CRTHitRecoAlg::reconfigure(const Config& config){

  fUseReadoutWindow = config.UseReadoutWindow(); 
  fQPed = config.QPed();
  fQSlope = config.QSlope();
  fNpeScaleShift = config.NpeScaleShift();
  fTimeCoincidenceLimit = config.TimeCoincidenceLimit();
  fClockSpeedCRT = config.ClockSpeedCRT();

  return;
}

std::map<std::pair<std::string, unsigned>, std::vector<CRTStrip>> CRTHitRecoAlg::CreateTaggerStrips(detinfo::DetectorClocksData const& clockData,
                                                                                                    detinfo::DetectorPropertiesData const& detProp,
                                                                                                    std::vector<art::Ptr<sbnd::crt::CRTData>> crtList,
												    art::Handle<std::vector<sbnd::crt::CRTData>> crtListHandle,
												    const art::Event &event,
												    TTree *fHitTree){

  double readoutWindowMuS  = clockData.TPCTick2Time((double)detProp.ReadOutWindowSize()); // [us]
  double driftTimeMuS = fTpcGeo.MaxX()/detProp.DriftVelocity(); // [us]

  std::map<std::pair<std::string, unsigned>, std::vector<CRTStrip>> taggerStrips;

  art::FindManyP<sim::AuxDetIDE> dataToIDEsAssn(crtListHandle, event, "crt");
  
  std::vector<std::string> *tStripName = 0;
  std::vector<unsigned> *tChannel = 0, *tPulseT0 = 0, *tPulseT1 = 0, *tPulseADC = 0;
  double tPE1, tPE2;
  std::string *tTaggerName;
  std::vector<double> *tStripTrueMin = 0, *tStripTrueMax = 0;
  unsigned tPlane;
  double tHitT0, tHitX, tHitXErr, tHitPEs;
  std::vector<double> *tLimits = 0;
  unsigned nIDEs;
  std::vector<double> *tEntryX = 0, *tEntryY = 0, *tEntryZ = 0, *tEntryT = 0, 
    *tExitX = 0, *tExitY = 0, *tExitZ = 0, *tExitT = 0;

  fHitTree->SetBranchAddress("StripName",&tStripName);
  fHitTree->SetBranchAddress("Channel",&tChannel);
  fHitTree->SetBranchAddress("PulseT0",&tPulseT0);
  fHitTree->SetBranchAddress("PulseT1",&tPulseT1);
  fHitTree->SetBranchAddress("PulseADC",&tPulseADC);
  fHitTree->SetBranchAddress("PE1",&tPE1);
  fHitTree->SetBranchAddress("PE2",&tPE2);
  fHitTree->SetBranchAddress("StripTrueMin",&tStripTrueMin);
  fHitTree->SetBranchAddress("StripTrueMax",&tStripTrueMax);
  fHitTree->SetBranchAddress("TaggerName",&tTaggerName);
  fHitTree->SetBranchAddress("Plane",&tPlane);
  fHitTree->SetBranchAddress("HitT0",&tHitT0);
  fHitTree->SetBranchAddress("HitX",&tHitX);
  fHitTree->SetBranchAddress("HitXErr",&tHitXErr);
  fHitTree->SetBranchAddress("HitPEs",&tHitPEs);
  fHitTree->SetBranchAddress("Limits",&tLimits);
  fHitTree->SetBranchAddress("nIDEs",&nIDEs);
  fHitTree->SetBranchAddress("EntryX",&tEntryX);
  fHitTree->SetBranchAddress("EntryY",&tEntryY);
  fHitTree->SetBranchAddress("EntryZ",&tEntryZ);
  fHitTree->SetBranchAddress("EntryT",&tEntryT);
  fHitTree->SetBranchAddress("ExitX",&tExitX);
  fHitTree->SetBranchAddress("ExitY",&tExitY);
  fHitTree->SetBranchAddress("ExitZ",&tExitZ);
  fHitTree->SetBranchAddress("ExitT",&tExitT);
  
  for (size_t i = 0; i < crtList.size(); i+=2){

    // Get the time
    //fTrigClock.SetTime(crtList[i]->T0());
    //double t1 = fTrigClock.Time(); // [us]
    double t1 = (double)(int)crtList[i]->T0()/fClockSpeedCRT; // [tick -> us]
    if(fUseReadoutWindow){
      if(!(t1 >= -driftTimeMuS && t1 <= readoutWindowMuS)) continue;
    }

    tStripName->clear(); tChannel->clear(); tPulseT0->clear(); tPulseT1->clear(); tPulseADC->clear(); tStripTrueMin->clear(); tStripTrueMax->clear();
    tTaggerName = new std::string(""); tPlane = 99999; tHitT0 = -99999.; tHitX = -99999.; tHitXErr = -99999.; tHitPEs = -99999.;
    tPE1 = -99999; tPE2 = -99999;
    tLimits->clear(); nIDEs = 99999; tEntryX->clear(); tEntryY->clear(); tEntryZ->clear(); tEntryT->clear(); tExitX->clear(); tExitY->clear(); tExitZ->clear(); tExitT->clear();

    CRTStrip strip = CreateCRTStrip(crtList[i], crtList[i+1], i);

    tStripName->push_back(fCrtGeo.ChannelToStripName(crtList[i]->Channel()));
    tStripName->push_back(fCrtGeo.ChannelToStripName(crtList[i+1]->Channel()));
    tChannel->push_back(crtList[i]->Channel());
    tChannel->push_back(crtList[i+1]->Channel());
    tPulseT0->push_back(crtList[i]->T0());
    tPulseT0->push_back(crtList[i+1]->T0());
    tPulseT1->push_back(crtList[i]->T1());
    tPulseT1->push_back(crtList[i+1]->T1());
    tPulseADC->push_back(crtList[i]->ADC());
    tPulseADC->push_back(crtList[i+1]->ADC());
    tPE1 = ((double)crtList[i]->ADC() - fQPed)/fQSlope;
    tPE2 = ((double)crtList[i+1]->ADC() - fQPed)/fQSlope;
    
    sbnd::CRTStripGeo stripGeo = fCrtGeo.GetStrip(tStripName->at(0));
    tStripTrueMin = new std::vector<double>({stripGeo.minX, stripGeo.minY, stripGeo.minZ});
    tStripTrueMax = new std::vector<double>({stripGeo.maxX, stripGeo.maxY, stripGeo.maxZ});
    tTaggerName = new std::string(strip.tagger.first);
    tPlane = strip.tagger.second;
    tHitT0 = strip.t0;
    tHitX = strip.x;
    tHitXErr = strip.ex;
    tHitPEs = strip.pes;

    tLimits = new std::vector<double>(fCrtGeo.StripLimitsWithChargeSharing(tStripName->at(0), strip.x, strip.ex));

    std::vector<art::Ptr<sim::AuxDetIDE> > ides = dataToIDEsAssn.at(crtList[i].key());
    nIDEs = ides.size();
    for (auto ide : ides) {
      tEntryX->push_back(ide->entryX);
      tEntryY->push_back(ide->entryY);
      tEntryZ->push_back(ide->entryZ);
      tEntryT->push_back(ide->entryT);
      tExitX->push_back(ide->exitX);
      tExitY->push_back(ide->exitY);
      tExitZ->push_back(ide->exitZ);
      tExitT->push_back(ide->exitT);
    }
    taggerStrips[strip.tagger].push_back(strip);
    fHitTree->Fill();

  }

  return taggerStrips;

}


CRTStrip CRTHitRecoAlg::CreateCRTStrip(art::Ptr<sbnd::crt::CRTData> sipm1, art::Ptr<sbnd::crt::CRTData> sipm2, size_t ind){

  // Get the time, channel, center and width
  //fTrigClock.SetTime(sipm1->T0());
  //double t1 = fTrigClock.Time(); // [us]
  double t1 = (double)(int)sipm1->T0()/fClockSpeedCRT; // [tick -> us]

  // Get strip info from the geometry service
  uint32_t channel = sipm1->Channel();

  std::pair<std::string,unsigned> tagger = ChannelToTagger(channel);

  // Get the time of hit on the second SiPM
  //fTrigClock.SetTime(sipm2->T0());
  //double t2 = fTrigClock.Time(); // [us]
  double t2 = (double)(int)sipm2->T0()/fClockSpeedCRT; // [tick -> us]

  // Calculate the number of photoelectrons at each SiPM
  double npe1 = ((double)sipm1->ADC() - fQPed)/fQSlope;
  double npe2 = ((double)sipm2->ADC() - fQPed)/fQSlope;

  // Calculate the distance between the SiPMs
  std::pair<double, double> sipmDist = DistanceBetweenSipms(sipm1, sipm2);

  double time = (t1 + t2)/2.;

  CRTStrip stripHit = {time, channel, sipmDist.first, sipmDist.second, npe1+npe2, tagger, ind};
  return stripHit;

}

std::pair<double, double> CRTHitRecoAlg::DistanceBetweenSipms(art::Ptr<sbnd::crt::CRTData> sipm1, art::Ptr<sbnd::crt::CRTData> sipm2){
  
  uint32_t channel = sipm1->Channel();
  std::string stripName = fCrtGeo.ChannelToStripName(channel);
  double width = fCrtGeo.GetStrip(stripName).width;

  // Calculate the number of photoelectrons at each SiPM
  double npe1 = ((double)sipm1->ADC() - fQPed)/fQSlope;
  double npe2 = ((double)sipm2->ADC() - fQPed)/fQSlope;

  // Calculate the distance between the SiPMs
  double x = (width/2.)*atan(log(1.*npe2/npe1)) + (width/2.);

  // Calculate the error
  double normx = x + 0.344677*x - 1.92045;
  double ex = 1.92380e+00+1.47186e-02*normx-5.29446e-03*normx*normx;

  return std::make_pair(x, ex);

}


std::vector<std::pair<sbn::crt::CRTHit, std::vector<int>>> CRTHitRecoAlg::CreateCRTHits(std::map<std::pair<std::string, unsigned>, std::vector<CRTStrip>> taggerStrips,
											const std::vector<art::Ptr<sbnd::crt::CRTData>> &crtList, 
											const art::Handle< std::vector<sbnd::crt::CRTData>> &crtListHandle,
											const art::Event &event,
											TTree *fThreeDHitTree)
{

  std::vector<std::pair<sbn::crt::CRTHit, std::vector<int>>> returnHits;

  std::vector<uint8_t> tfeb_id = {0};
  std::map<uint8_t, std::vector<std::pair<int,float>>> tpesmap;
  tpesmap[0] = {std::make_pair(0,0)};
  
  // Remove any duplicate (same channel and time) hit strips
  for(auto &tagStrip : taggerStrips){
    std::sort(tagStrip.second.begin(), tagStrip.second.end(),
              [](const CRTStrip & a, const CRTStrip & b) -> bool{
                return (a.t0 < b.t0) || 
                       ((a.t0 == b.t0) && (a.channel < b.channel));
              });
    // Remove hits with the same time and channel
    tagStrip.second.erase(std::unique(tagStrip.second.begin(), tagStrip.second.end(),
                                         [](const CRTStrip & a, const CRTStrip & b) -> bool{
                                           return a.t0 == b.t0 && a.channel == b.channel;
                                          }), tagStrip.second.end());
  }

  std::vector<std::string> usedTaggers;

  for (auto &tagStrip : taggerStrips){
    if (std::find(usedTaggers.begin(),usedTaggers.end(),tagStrip.first.first)!=usedTaggers.end()) continue;
    usedTaggers.push_back(tagStrip.first.first);
    unsigned planeID = 0;
    if(tagStrip.first.second==0) planeID = 1;
    std::pair<std::string,unsigned> otherPlane = std::make_pair(tagStrip.first.first, planeID);

    for (size_t hit_i = 0; hit_i < tagStrip.second.size(); hit_i++){
      // Get the position (in real space) of the 4 corners of the hit, taking charge sharing into account
      std::vector<double> limits1 =  ChannelToLimits(tagStrip.second[hit_i]);

      // Check for overlaps on the first plane
      if(CheckModuleOverlap(tagStrip.second[hit_i].channel)){

        // Loop over all the hits on the parallel (odd) plane
        for (size_t hit_j = 0; hit_j < taggerStrips[otherPlane].size(); hit_j++){
          // Get the limits in the two variable directions
          std::vector<double> limits2 = ChannelToLimits(taggerStrips[otherPlane][hit_j]);

          // If the time and position match then record the pair of hits
          std::vector<double> overlap = CrtOverlap(limits1, limits2);
          double t0_1 = tagStrip.second[hit_i].t0;
          double t0_2 = taggerStrips[otherPlane][hit_j].t0;
          if (overlap[0] != -99999 && std::abs(t0_1 - t0_2) < fTimeCoincidenceLimit){
            // Calculate the mean and error in x, y, z
            TVector3 mean((overlap[0] + overlap[1])/2., 
                          (overlap[2] + overlap[3])/2., 
                          (overlap[4] + overlap[5])/2.);
            TVector3 error(std::abs((overlap[1] - overlap[0])/2.), 
                           std::abs((overlap[3] - overlap[2])/2.), 
                           std::abs((overlap[5] - overlap[4])/2.));

            // Average the time
            double time = (t0_1 + t0_2)/2;
            //double pes = tagStrip.second[hit_i].pes + taggerStrips[otherPlane][hit_j].pes;
            double pes = CorrectNpe(tagStrip.second[hit_i], taggerStrips[otherPlane][hit_j], mean);

            // Create a CRT hit
            sbn::crt::CRTHit crtHit = FillCrtHit(tfeb_id, tpesmap, pes, time, 0, mean.X(), error.X(), 
                                            mean.Y(), error.Y(), mean.Z(), error.Z(), tagStrip.first.first);
            std::vector<int> dataIds;
            dataIds.push_back(tagStrip.second[hit_i].dataID);
            dataIds.push_back(tagStrip.second[hit_i].dataID+1);
            dataIds.push_back(taggerStrips[otherPlane][hit_j].dataID);
            dataIds.push_back(taggerStrips[otherPlane][hit_j].dataID+1);
            returnHits.push_back(std::make_pair(crtHit, dataIds));

	    FillThreeDTree(event, crtHit, dataIds, crtList, crtListHandle, fThreeDHitTree, t0_1, t0_2);
	   }

        }

      }
      // If module doesn't overlap with a perpendicular one create 1D hits
      else{
        TVector3 mean((limits1[0] + limits1[1])/2., 
                      (limits1[2] + limits1[3])/2., 
                      (limits1[4] + limits1[5])/2.);
        TVector3 error(std::abs((limits1[1] - limits1[0])/2.), 
                       std::abs((limits1[3] - limits1[2])/2.), 
                       std::abs((limits1[5] - limits1[4])/2.));

        double time = tagStrip.second[hit_i].t0;
        double pes = tagStrip.second[hit_i].pes;

        // Just use the single plane limits as the crt hit
        sbn::crt::CRTHit crtHit = FillCrtHit(tfeb_id, tpesmap, pes, time, 0, mean.X(), error.X(), 
                                        mean.Y(), error.Y(), mean.Z(), error.Z(), tagStrip.first.first);
        std::vector<int> dataIds;
        dataIds.push_back(tagStrip.second[hit_i].dataID);
        dataIds.push_back(tagStrip.second[hit_i].dataID+1);
        returnHits.push_back(std::make_pair(crtHit, dataIds));
	
	FillThreeDTree(event, crtHit, dataIds, crtList, crtListHandle, fThreeDHitTree, time, time);
      }

    }
    // Loop over tagger modules on the perpendicular plane to look for 1D hits
    for (size_t hit_j = 0; hit_j < taggerStrips[otherPlane].size(); hit_j++){
      // Get the limits in the two variable directions
      std::vector<double> limits1 = ChannelToLimits(taggerStrips[otherPlane][hit_j]);

      // Check if module overlaps with a perpendicular one
      if(!CheckModuleOverlap(taggerStrips[otherPlane][hit_j].channel)){
        TVector3 mean((limits1[0] + limits1[1])/2., 
                      (limits1[2] + limits1[3])/2., 
                      (limits1[4] + limits1[5])/2.);
        TVector3 error(std::abs((limits1[1] - limits1[0])/2.), 
                       std::abs((limits1[3] - limits1[2])/2.), 
                       std::abs((limits1[5] - limits1[4])/2.));

        double time = taggerStrips[otherPlane][hit_j].t0;
        double pes = taggerStrips[otherPlane][hit_j].pes;

        // Just use the single plane limits as the crt hit
        sbn::crt::CRTHit crtHit = FillCrtHit(tfeb_id, tpesmap, pes, time, 0, mean.X(), error.X(), 
                                        mean.Y(), error.Y(), mean.Z(), error.Z(), otherPlane.first);
        std::vector<int> dataIds;
        dataIds.push_back(taggerStrips[otherPlane][hit_j].dataID);
        dataIds.push_back(taggerStrips[otherPlane][hit_j].dataID+1);
        returnHits.push_back(std::make_pair(crtHit, dataIds));

	FillThreeDTree(event, crtHit, dataIds, crtList, crtListHandle, fThreeDHitTree, time, time);
      }

    }    
  }

  return returnHits;

}

void CRTHitRecoAlg::FillThreeDTree(const art::Event &event, const sbn::crt::CRTHit &crtHit, const std::vector<int> dataIDs,
				   const std::vector<art::Ptr<sbnd::crt::CRTData>> &crtList, const art::Handle< std::vector<sbnd::crt::CRTData>> &crtListHandle, 
				   TTree *fThreeDHitTree, double &t0_1, double &t0_2)
  {
    float tPEsHit = -99999.f, tPosX = -99999.f, tPosY = -99999.f, tPosZ = -99999.f, tErrX = -99999.f, tErrY = -99999.f, tErrZ = -99999.f;
    unsigned tTime = 99999, nIDEs = 0;
    std::string *tTagger = 0;
    std::vector<double> *tEntryX = 0, *tEntryY = 0, *tEntryZ = 0, *tEntryT = 0, 
      *tExitX = 0, *tExitY = 0, *tExitZ = 0, *tExitT = 0;
    double tTimeHit0 = -99999.f, tTimeHit1 = -99999.f, tDistHit0 = -99999.f, tDistHit1 = -99999.f;

    fThreeDHitTree->SetBranchAddress("PEsHit",&tPEsHit);
    fThreeDHitTree->SetBranchAddress("PosX",&tPosX);
    fThreeDHitTree->SetBranchAddress("PosY",&tPosY);
    fThreeDHitTree->SetBranchAddress("PosZ",&tPosZ);
    fThreeDHitTree->SetBranchAddress("ErrX",&tErrX);
    fThreeDHitTree->SetBranchAddress("ErrY",&tErrY);
    fThreeDHitTree->SetBranchAddress("ErrZ",&tErrZ);
    fThreeDHitTree->SetBranchAddress("Time",&tTime);
    fThreeDHitTree->SetBranchAddress("tTagger",&tTagger);
    fThreeDHitTree->SetBranchAddress("nIDEs",&nIDEs);
    fThreeDHitTree->SetBranchAddress("EntryX",&tEntryX);
    fThreeDHitTree->SetBranchAddress("EntryY",&tEntryY);
    fThreeDHitTree->SetBranchAddress("EntryZ",&tEntryZ);
    fThreeDHitTree->SetBranchAddress("EntryT",&tEntryT);
    fThreeDHitTree->SetBranchAddress("ExitX",&tExitX);
    fThreeDHitTree->SetBranchAddress("ExitY",&tExitY);
    fThreeDHitTree->SetBranchAddress("ExitZ",&tExitZ);
    fThreeDHitTree->SetBranchAddress("ExitT",&tExitT);
    fThreeDHitTree->SetBranchAddress("TimeHit0",&tTimeHit0);
    fThreeDHitTree->SetBranchAddress("TimeHit1",&tTimeHit1);
    fThreeDHitTree->SetBranchAddress("DistHit0",&tDistHit0);
    fThreeDHitTree->SetBranchAddress("DistHit1",&tDistHit1);

    tPEsHit = crtHit.peshit;
    tPosX = crtHit.x_pos;
    tPosY = crtHit.y_pos;
    tPosZ = crtHit.z_pos;
    tErrX = crtHit.x_err;
    tErrY = crtHit.y_err;
    tErrZ = crtHit.z_err;
    tTime = crtHit.ts0_ns;
    tTagger = new std::string(crtHit.tagger);

    if(dataIDs.size() == 4)
      {
	std::string name0 = fCrtGeo.ChannelToStripName(crtList[dataIDs[0]]->Channel());
	std::string name1 = fCrtGeo.ChannelToStripName(crtList[dataIDs[3]]->Channel());
	tDistHit0 = fCrtGeo.DistanceDownStrip({tPosX, tPosY, tPosZ}, name0);
	tDistHit1 = fCrtGeo.DistanceDownStrip({tPosX, tPosY, tPosZ}, name1);
	tTimeHit0 = t0_1;
	tTimeHit1 = t0_2;
      }
    else if(dataIDs.size() == 2)
      {
	std::string name0 = fCrtGeo.ChannelToStripName(crtList[dataIDs[0]]->Channel());
	tDistHit0 = fCrtGeo.DistanceDownStrip({tPosX, tPosY, tPosZ}, name0);
	tTimeHit0 = t0_1;
      }

    art::FindManyP<sim::AuxDetIDE> dataToIDEsAssn(crtListHandle, event, "crt");
    for(auto const &dataID : dataIDs)
      {
	std::vector<art::Ptr<sim::AuxDetIDE>> ides = dataToIDEsAssn.at(crtList[dataID].key());

	for(auto const &ide : ides)
	  {
	    ++nIDEs;
	    tEntryX->push_back(ide->entryX);
	    tEntryY->push_back(ide->entryY);
	    tEntryZ->push_back(ide->entryZ);
	    tEntryT->push_back(ide->entryT);
	    tExitX->push_back(ide->exitX);
	    tExitY->push_back(ide->exitY);
	    tExitZ->push_back(ide->exitZ);
	    tExitT->push_back(ide->exitT);
	  }
      }
    fThreeDHitTree->Fill();
  }

// Function to calculate the strip position limits in real space from channel
std::vector<double> CRTHitRecoAlg::ChannelToLimits(CRTStrip stripHit){

  std::string stripName = fCrtGeo.ChannelToStripName(stripHit.channel);
  return fCrtGeo.StripLimitsWithChargeSharing(stripName, stripHit.x, stripHit.ex);

} // CRTHitRecoAlg::ChannelToLimits()


// Function to calculate the overlap between two crt strips
std::vector<double> CRTHitRecoAlg::CrtOverlap(std::vector<double> strip1, std::vector<double> strip2){

  // Get the minimum and maximum X, Y, Z coordinates
  double minX = std::max(strip1[0], strip2[0]);
  double maxX = std::min(strip1[1], strip2[1]);
  double minY = std::max(strip1[2], strip2[2]);
  double maxY = std::min(strip1[3], strip2[3]);
  double minZ = std::max(strip1[4], strip2[4]);
  double maxZ = std::min(strip1[5], strip2[5]);

  std::vector<double> null = {-99999, -99999, -99999, -99999, -99999, -99999};
  std::vector<double> overlap = {minX, maxX, minY, maxY, minZ, maxZ};

  // If the two strips overlap in 2 dimensions then return the overlap
  if ((minX<maxX && minY<maxY) || (minX<maxX && minZ<maxZ) || (minY<maxY && minZ<maxZ)) return overlap;
  // Otherwise return a "null" value
  return null;

} // CRTHitRecoAlg::CRTOverlap()


// Function to return the CRT tagger name and module position from the channel ID
std::pair<std::string,unsigned> CRTHitRecoAlg::ChannelToTagger(uint32_t channel){

  std::string stripName = fCrtGeo.ChannelToStripName(channel);
  size_t planeID = fCrtGeo.GetModule(fCrtGeo.GetStrip(stripName).module).planeID;
  std::string tagName = fCrtGeo.GetModule(fCrtGeo.GetStrip(stripName).module).tagger;
  
  std::pair<std::string, unsigned> output = std::make_pair(tagName, planeID);

  return output;

} // CRTHitRecoAlg::ChannelToTagger()


// Function to check if a CRT strip overlaps with a perpendicular module
bool CRTHitRecoAlg::CheckModuleOverlap(uint32_t channel){

  std::string stripName = fCrtGeo.ChannelToStripName(channel);
  return fCrtGeo.StripHasOverlap(stripName);

} // CRTHitRecoAlg::CheckModuleOverlap


// Function to make filling a CRTHit a bit faster
sbn::crt::CRTHit CRTHitRecoAlg::FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, 
                              std::vector<std::pair<int,float>>> tpesmap, float peshit, double time, int plane, 
                              double x, double ex, double y, double ey, double z, double ez, std::string tagger){

  sbn::crt::CRTHit crtHit;

  crtHit.feb_id      = tfeb_id;
  crtHit.pesmap      = tpesmap;
  crtHit.peshit      = peshit;
  crtHit.ts0_s_corr  = 0; 
  crtHit.ts0_ns      = time * 1e3;
  crtHit.ts0_ns_corr = 0; 
  crtHit.ts1_ns      = time * 1e3;
  crtHit.ts0_s       = time * 1e-6;
  crtHit.plane       = plane;
  crtHit.x_pos       = x;
  crtHit.x_err       = ex;
  crtHit.y_pos       = y; 
  crtHit.y_err       = ey;
  crtHit.z_pos       = z;
  crtHit.z_err       = ez;
  crtHit.tagger      = tagger;

  return crtHit;

} // CRTHitRecoAlg::FillCrtHit()


// Function to correct number of photoelectrons by distance down strip
double CRTHitRecoAlg::CorrectNpe(CRTStrip strip1, CRTStrip strip2, TVector3 position){
  geo::Point_t pos {position.X(), position.Y(), position.Z()};

  // Get the strip name from the channel ID
  std::string name1 = fCrtGeo.ChannelToStripName(strip1.channel);
  std::string name2 = fCrtGeo.ChannelToStripName(strip2.channel);

  // Get the distance from the CRT hit to the sipm end
  double stripDist1 = fCrtGeo.DistanceDownStrip(pos, name1);
  double stripDist2 = fCrtGeo.DistanceDownStrip(pos, name2);

  // Correct the measured pe
  double pesCorr1 = strip1.pes * pow(stripDist1 - fNpeScaleShift, 2) / pow(fNpeScaleShift, 2);
  double pesCorr2 = strip2.pes * pow(stripDist2 - fNpeScaleShift, 2) / pow(fNpeScaleShift, 2);

  // Add the two strips together
  return pesCorr1 + pesCorr2;
}

}

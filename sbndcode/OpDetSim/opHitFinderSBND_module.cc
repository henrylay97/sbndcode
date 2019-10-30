////////////////////////////////////////////////////////////////////////
// Class:       opHitFinderSBND
// Module Type: producer
// File:        opHitFinderSBND_module.cc
//
// This module produces an OpHit object for light analysis
// Created by L. Paulucci and F. Marinho
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
//#include "larana/OpticalDetector/OpHitFinder/OpHitAlg.h"
//#include "lardataobj/Simulation/BeamGateInfo.h"

#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "sbndcode/Utilities/SignalShapingServiceSBND.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
//#include "larsim/MCCheater/PhotonBackTracker.h"

#include <memory>
#include <algorithm>
#include "TMath.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TF1.h"

#include "sbndPDMapAlg.h" 

namespace opdet{

  class opHitFinderSBND;

  class opHitFinderSBND : public art::EDProducer {
  public:
    explicit opHitFinderSBND(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
    opHitFinderSBND(opHitFinderSBND const &) = delete;
    opHitFinderSBND(opHitFinderSBND &&) = delete;
    opHitFinderSBND & operator = (opHitFinderSBND const &) = delete;
    opHitFinderSBND & operator = (opHitFinderSBND &&) = delete;

  // Required functions.
    void produce(art::Event & e) override;
    opdet::sbndPDMapAlg map; //map for photon detector types

  private:

  // Declare member data here.
    std::string fInputModuleName;
  //  art::ServiceHandle<cheat::PhotonBackTracker> pbt;
    double fSampling; //in MHz
    double fBaselineSample; //in ticks
    double fUseDenoising; 
    double fPulsePolarityPMT; 
    double fPulsePolarityArapuca; 
    double fSaturation; //in number of p.e.
    double fArea1pePMT; //area of 1 pe in ADC*ns for PMTs
    double fArea1peSiPM; //area of 1 pe in ADC*ns for Arapucas
    int fThresholdPMT; //in ADC
    int fThresholdArapuca; //in ADC
    int fEvNumber;
    int fChNumber;
    std::vector<short unsigned int> fwaveform;
    //int fSize;
    //int fTimePMT;         //Start time of PMT signal
    //int fTimeMax;         //Time of maximum (minimum) PMT signal
    void subtractBaseline(std::vector<short unsigned int>& waveform, std::string pdtype, double& rms);
    bool findPeak(std::vector<short unsigned int>& waveform, size_t& time, double& Area, double rms, double& amplitude, std::string type);
    void denoise(std::vector<short unsigned int>& waveform);
    void TV1D_denoise(float* input, float*& output, const int width, const float lambda);
//    std::stringstream histname;
  };

  opHitFinderSBND::opHitFinderSBND(fhicl::ParameterSet const & p)
   : EDProducer{p}
// Initialize member data here.
  {
    fInputModuleName = p.get< std::string >("InputModule" );
    fBaselineSample  = p.get< int    >("BaselineSample"); //in ticks
    fSaturation      = p.get< double >("Saturation"   ); //in number of p.e.
    fArea1pePMT      = p.get< double >("Area1pePMT"   ); //in ADC*ns for PMTs
    fArea1peSiPM     = p.get< double >("Area1peSiPM"  ); //in ADC*ns for SiPMs
    fThresholdPMT    = p.get< double >("ThresholdPMT" ); //in ADC
    fThresholdArapuca = p.get< double >("ThresholdArapuca"); //in ADC
    fUseDenoising     = p.get< int   >("UseDenoising"); 
    fPulsePolarityPMT = p.get< int   >("PulsePolarityPMT");
    fPulsePolarityArapuca = p.get< int   >("PulsePolarityArapuca");

    auto const *timeService = lar::providerFrom< detinfo::DetectorClocksService >();
    fSampling = (timeService->OpticalClock().Frequency()); // MHz

  // Call appropriate produces<>() functions here.
    produces<std::vector<recob::OpHit>>();
  }

  void opHitFinderSBND::produce(art::Event & e)
  {
    // Implementation of required member function here.
    fEvNumber = e.id().event();
    std::cout << "Event #" << fEvNumber << std::endl;

    std::unique_ptr< std::vector< recob::OpHit > > pulseVecPtr(std::make_unique< std::vector< recob::OpHit > > ());
    fwaveform.reserve(20000); // TODO: no hardcoded value

    art::ServiceHandle<art::TFileService> tfs;
    //   histname.str(std::string());
    // histname << "event_" << fEvNumber;
    art::Handle< std::vector< raw::OpDetWaveform > > wvfHandle;
    e.getByLabel(fInputModuleName, wvfHandle);

    if(!wvfHandle.isValid()){
      std::cout <<Form("Did not find any waveform") << std::endl;
    }

    size_t timebin=0;
    double FWHM=1, Area=0, phelec, fasttotal=3./4., rms=0, amplitude=0, time=0;
    unsigned short frame=1;
//    int histogram_number = 0;
    for(auto const& wvf : (*wvfHandle)){
      if (wvf.size() == 0 ) {
        std::cout << "Empty waveform, continue." << std::endl;
        continue;
      }

      fChNumber = wvf.ChannelNumber();
      fwaveform.resize(wvf.size());
      for(unsigned int i=0;i<wvf.size();i++){
        fwaveform[i]=wvf[i];
      }
      subtractBaseline(fwaveform, map.pdName(fChNumber), rms);
      if((map.pdName(fChNumber)=="pmt") || (map.pdName(fChNumber)== "barepmt")){
      }else{
        if(fUseDenoising==1) denoise(fwaveform);
      }
      int i=1;
      while(findPeak(fwaveform,timebin,Area,rms,amplitude,map.pdName(fChNumber))){
        time = wvf.TimeStamp() + (double)timebin/fSampling;

        if(map.pdName(fChNumber)=="pmt" || map.pdName(fChNumber) == "barepmt"){
          phelec=Area/fArea1pePMT;
          //   std::cout << 0 << " " << time << " " << Area << " " << phelec << std::endl;
        }else{
          phelec=Area/fArea1peSiPM;
          //  std::cout << 1 << " " << time << " " << Area << " " << phelec << std::endl;
        }
        i++;
        recob::OpHit opHit(fChNumber, time, time, frame, FWHM, Area, amplitude, phelec, fasttotal);//including hit info: OpChannel, PeakTime, PeakTimeAbs, Frame, Width, Area, PeakHeight, PE, FastToTotal
        pulseVecPtr->emplace_back(opHit);
      }
      //     histogram_number += 1;
      // fwaveform.clear();
    } // for(auto const& wvf : (*wvfHandle)){
    e.put(std::move(pulseVecPtr));
  } // void opHitFinderSBND::produce(art::Event & e)

  DEFINE_ART_MODULE(opHitFinderSBND)

  void opHitFinderSBND::subtractBaseline(std::vector<short unsigned int>& waveform, std::string pdtype, double& rms){
    double baseline = 0.0; 
    rms=0.0;
    int cnt = 0;
    for(int i=0; i<fBaselineSample; i++){
	baseline+=waveform[i];
	rms+=pow(waveform[i],2.0);
        cnt++;
    }
    baseline=baseline/cnt;
    rms=sqrt(rms/cnt-baseline*baseline);
    rms=rms/sqrt(cnt-1);

    if(pdtype=="pmt" || pdtype == "barepmt"){
        for(unsigned int i=0; i<waveform.size(); i++) waveform[i]=fPulsePolarityPMT*(waveform[i]-baseline);
    }else{
        for(unsigned int i=0; i<waveform.size(); i++) waveform[i]=fPulsePolarityArapuca*(waveform[i]-baseline);
    }
  }

  bool opHitFinderSBND::findPeak(std::vector<short unsigned int>& waveform, size_t& time, double& Area, double rms, double& amplitude, std::string type){

    //Gets info from highest peak and suppress it 
    double aux = *max_element(waveform.begin(), waveform.end());
    double max;
    size_t time_end, bin, binmax = distance(waveform.begin(),max_element(waveform.begin(), waveform.end()));
    int threshold;
    Area=0;

    if(type=="pmt" || type == "barepmt"){
      threshold=fThresholdPMT;
    }else{
      threshold=fThresholdArapuca;
    }

    bin=binmax;
    amplitude=aux;
    max=aux;

    if(aux<threshold)return false;

    while(aux>=rms){
        bin++;	
        aux = waveform[bin];
    }
    time_end=bin-1; //looking for the length of the peak

    aux=max;
    while(aux>=rms){
        bin--;
        aux = waveform[bin];
    }
    time = bin+1; //for rise time

    for(unsigned int j=time; j<=time_end; j++) Area+=waveform[j];
    Area=Area/fSampling;

    bin = time;
    aux = waveform[time];  

   while(aux>=rms){
        waveform[bin]=0.0;
        bin++;
        aux = waveform[bin];
    }
    //std::cout << time << " " << time_end << " " << (time_end - time);
    time=binmax; //returning the peak time
    return true;		
  }

  void opHitFinderSBND::denoise(std::vector<short unsigned int>& waveform){

    int wavelength = waveform.size();
    float lambda = 10.0;
    float* input;
    float* output;

    //std::cout <<"denoise wavl:" << wavelength << std::endl;

    input = (float*)malloc(sizeof(*input)*wavelength);
    output = (float*)malloc(sizeof(*output)*wavelength);

    for(int i=0; i<wavelength; i++){
      input[i] = waveform[i];
      output[i] = input[i];
    }

    //std::cout <<"denoise full array" << std::endl;

    TV1D_denoise(input, output, (const int)wavelength, (const float)lambda);

    //std::cout <<"denoise denoise" << std::endl;

    for(int i=0; i<wavelength; i++){
      if(output[i]) waveform[i]=output[i];
    }

 }

  void opHitFinderSBND::TV1D_denoise(float* input, float*& output, const int width, const float lambda) {
	if (width>0) {				/*to avoid invalid memory access to input[0]*/
		int k=0, k0=0;			/*k: current sample location, k0: beginning of current segment*/
		float umin=lambda, umax=-lambda;	/*u is the dual variable*/
		float vmin=input[0]-lambda, vmax=input[0]+lambda;	/*bounds for the segment's value*/
		int kplus=0, kminus=0; 	/*last positions where umax=-lambda, umin=lambda, respectively*/
		const float twolambda=2.0*lambda;	/*auxiliary variable*/
		const float minlambda=-lambda;		/*auxiliary variable*/
		for (;;) {				/*simple loop, the exit test is inside*/
			while (k==width-1) {	/*we use the right boundary condition*/
				if (umin<0.0) {			/*vmin is too high -> negative jump necessary*/
					do output[k0++]=vmin; while (k0<=kminus);
					umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
				} else if (umax>0.0) {	/*vmax is too low -> positive jump necessary*/
					do output[k0++]=vmax; while (k0<=kplus);
					umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
				} else {
					vmin+=umin/(k-k0+1); 
					do output[k0++]=vmin; while(k0<=k); 
					return;
				}
			}
			if ((umin+=input[k+1]-vmin)<minlambda) {		/*negative jump necessary*/
				do output[k0++]=vmin; while (k0<=kminus);
				vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
				umin=lambda; umax=minlambda;
			} else if ((umax+=input[k+1]-vmax)>lambda) {	/*positive jump necessary*/
				do output[k0++]=vmax; while (k0<=kplus);
				vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
				umin=lambda; umax=minlambda;
			} else { 	/*no jump necessary, we continue*/
				k++;
				if (umin>=lambda) {		/*update of vmin*/
					vmin+=(umin-lambda)/((kminus=k)-k0+1);
					umin=lambda;
				} 
				if (umax<=minlambda) {	/*update of vmax*/
					vmax+=(umax+lambda)/((kplus=k)-k0+1);
					umax=minlambda;
				} 	
			}
		}
	}
  }

}

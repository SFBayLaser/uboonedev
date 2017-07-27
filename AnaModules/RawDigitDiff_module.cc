////////////////////////////////////////////////////////////////////////
//
// RawDigitDiff class - This attempts to compare RawDigits by taking their
//                      difference and outputting a new RawDigit of this diff
// usher@slac.stanford.edu
//
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>
#include <vector>
#include <utility> // std::pair<>
#include <memory> // std::unique_ptr<>
#include <iomanip>

// ROOT libraries
#include "TH1D.h"
#include "TProfile.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"

///creation of calibrated signals on wires
namespace caldata {

class RawDigitDiff : public art::EDProducer
{
public:
    
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit RawDigitDiff(fhicl::ParameterSet const& pset); 
    virtual ~RawDigitDiff();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
    
    void reconfFFT(int temp_fftsize);
private:
    
    art::InputTag fTopRawDigitLabel;     ///< module that made digits
    art::InputTag fBotRawDigitLabel;     ///< module that made digits
    
    std::vector<TH1D*>     fAveHistVec;
    std::vector<TH1D*>     fRMSHistVec;
    std::vector<TProfile*> fAveVsChanHistVec;
    std::vector<TProfile*> fRMSVsChanHistVec;
    
    const geo::GeometryCore*           fGeometry = lar::providerFrom<geo::Geometry>();
    const lariov::DetPedestalProvider& fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
}; // class RawDigitDiff

DEFINE_ART_MODULE(RawDigitDiff)
  
//-------------------------------------------------
RawDigitDiff::RawDigitDiff(fhicl::ParameterSet const& pset) :
    fPedestalRetrievalAlg(*lar::providerFrom<lariov::DetPedestalService>())
{
    this->reconfigure(pset);

    produces< std::vector<raw::RawDigit> >();
}
  
//-------------------------------------------------
RawDigitDiff::~RawDigitDiff()
{
}

//////////////////////////////////////////////////////
void RawDigitDiff::reconfigure(fhicl::ParameterSet const& p)
{
    fTopRawDigitLabel = p.get< std::string >("TopRawDigitLabel", "wcNoiseFilter");
    fBotRawDigitLabel = p.get< std::string >("BotRawDigitLabel", "digitfilter"  );

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;
    
    fAveHistVec.resize(3);
    fRMSHistVec.resize(3);
    fAveVsChanHistVec.resize(3);
    fRMSVsChanHistVec.resize(3);
    
    for(size_t planeIdx = 0; planeIdx < 3; planeIdx++)
    {
        fAveHistVec[planeIdx]       = tfs->make<TH1D>(Form("Ave_%02zu",planeIdx),";Average;", 200, -20., 20.);
        fRMSHistVec[planeIdx]       = tfs->make<TH1D>(Form("RMS_%02zu",planeIdx),";rms;",     100,   0., 10.);
        fAveVsChanHistVec[planeIdx] = tfs->make<TProfile>(Form("AveVsChan_%2zu",planeIdx), ";Wire #", fGeometry->Nwires(planeIdx), 0., fGeometry->Nwires(planeIdx), -20., 20.);
        fRMSVsChanHistVec[planeIdx] = tfs->make<TProfile>(Form("RMSVsChan_%2zu",planeIdx), ";Wire #", fGeometry->Nwires(planeIdx), 0., fGeometry->Nwires(planeIdx),   0., 10.);
    }

    return;
}

//-------------------------------------------------
void RawDigitDiff::beginJob()
{
} // beginJob

//////////////////////////////////////////////////////
void RawDigitDiff::endJob()
{  
}

//////////////////////////////////////////////////////
void RawDigitDiff::produce(art::Event& event)
{
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    std::unique_ptr<std::vector<raw::RawDigit> > filteredRawDigit(new std::vector<raw::RawDigit>);
    
    // Read in the top RawDigit
    art::Handle< std::vector<raw::RawDigit> > topDigitVecHandle;
    event.getByLabel(fTopRawDigitLabel, topDigitVecHandle);
    
    // Read in the top RawDigit
    art::Handle< std::vector<raw::RawDigit> > botDigitVecHandle;
    event.getByLabel(fBotRawDigitLabel, botDigitVecHandle);

    if (topDigitVecHandle->size() > 0 && botDigitVecHandle->size())
    {
        raw::RawDigit::ADCvector_t rawDigitVec;
        
        filteredRawDigit->reserve(topDigitVecHandle->size());
        
        // First task is to get channel to vector mapping for both collections since they can get out of step
        using ChannelToVectorMap = std::map<raw::ChannelID_t,const raw::RawDigit*>;
        
        ChannelToVectorMap topMap;
        
        for(size_t rdIter = 0; rdIter < topDigitVecHandle->size(); rdIter++)
        {
            const raw::RawDigit* rawDigit = &topDigitVecHandle->at(rdIter);
            
            topMap[rawDigit->Channel()] = rawDigit;
        }
        
        ChannelToVectorMap botMap;
        
        for(size_t rdIter = 0; rdIter < botDigitVecHandle->size(); rdIter++)
        {
            const raw::RawDigit* rawDigit = &botDigitVecHandle->at(rdIter);
            
            botMap[rawDigit->Channel()] = rawDigit;
        }
        
        // Now iterate through the top map
        for(const auto& topPair : topMap)
        {
            // get the reference to the current raw::RawDigit
            const raw::RawDigit* topDigitVec = topPair.second;
            raw::ChannelID_t     channel     = topPair.first;
            
            unsigned int dataSize = topDigitVec->Samples();
            
            // Resize output buffer
            rawDigitVec.resize(dataSize, 0);
            
            // vector holding uncompressed adc values
            std::vector<short> toprawadc(dataSize);

            // uncompress the data
            raw::Uncompress(topDigitVec->ADCs(), toprawadc, topDigitVec->Compression());
            
            // Can we find the corresponding channel in the bottom raw digit?
            ChannelToVectorMap::iterator botMapItr = botMap.find(channel);
            
            if (botMapItr != botMap.end())
            {
                const raw::RawDigit* botDigitVec = botMapItr->second;
                
                // vector holding uncompressed adc values
                std::vector<short> botrawadc(dataSize);
                
                // uncompress the data
                raw::Uncompress(botDigitVec->ADCs(), botrawadc, botDigitVec->Compression());
                
                // Recover the database version of the pedestal
                float pedestal = fPedestalRetrievalAlg.PedMean(channel);
                
                std::transform(toprawadc.begin(),toprawadc.end(),botrawadc.begin(),rawDigitVec.begin(),[pedestal](const auto& left, const auto& right){return std::round(left - right + pedestal);});
                
                raw::RawDigit::ADCvector_t locRawDigit = rawDigitVec;
                
                int meanSum = std::accumulate(locRawDigit.begin(),locRawDigit.end(),0);
                
                float mean    = float(meanSum) / float(locRawDigit.size());
                int   intMean = std::round(mean);
                
                std::transform(locRawDigit.begin(),locRawDigit.end(),locRawDigit.begin(),std::bind2nd(std::minus<short>(),intMean));
                
                // And now the local rms calculation...
                int   rmsSum = std::inner_product(locRawDigit.begin(), locRawDigit.end(), locRawDigit.begin(), 0);
                float rmsVal = std::sqrt(std::max(float(0.),float(rmsSum) / float(locRawDigit.size())));
                
                // First up, determine what kind of wire we have
                std::vector<geo::WireID> wids    = fGeometry->ChannelToWire(channel);
                const geo::PlaneID&      planeID = wids[0].planeID();
                
                fAveHistVec.at(planeID.Plane)->Fill(mean-pedestal, 1.);
                fRMSHistVec.at(planeID.Plane)->Fill(rmsVal, 1.);
                fAveVsChanHistVec.at(planeID.Plane)->Fill(wids[0].Wire, mean-pedestal, 1.);
                fRMSVsChanHistVec.at(planeID.Plane)->Fill(wids[0].Wire, rmsVal, 1.);
            }
            else
            {
                // Copy
                std::copy(toprawadc.begin(),toprawadc.end(),rawDigitVec.begin());
            }
            
            filteredRawDigit->emplace_back(raw::RawDigit(channel, rawDigitVec.size(), rawDigitVec, raw::kNone));
            filteredRawDigit->back().SetPedestal(topDigitVec->GetPedestal(),topDigitVec->GetSigma());
        }
    }
    
    // Add tracks and associations to event.
    event.put(std::move(filteredRawDigit));

    return;
} // produce


} // end namespace caldata

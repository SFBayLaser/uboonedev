////////////////////////////////////////////////////////////////////////
/// \file   DrawWire2D.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "IDraw2DObjects.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

#include "TEveManager.h"
#include "TEveEventManager.h"
#include "TEveScene.h"
#include "TEveViewer.h"
#include <TEveBrowser.h>
#include "TEveRGBAPaletteOverlay.h"
#include "TEveFrameBox.h"
#include "TEveQuadSet.h"
#include "TEveTrans.h"
#include "TEveProjectionAxes.h"
#include "TEveProjectionManager.h"
#include "TEveWindow.h"
#include "TGLViewer.h"
#include "TGLCameraOverlay.h"
#include "TGTab.h"

#include "TRandom.h"
#include "TStyle.h"

#include <cmath>
#include <fstream>

namespace display
{

class DrawWire2D : IDraw2DObjects
{
public:
    explicit DrawWire2D(const fhicl::ParameterSet& pset);
    
    ~DrawWire2D();
    
    void drawWire2D(const art::Event&) const override;
    
private:
    // Member variables from the fhicl file
    int      fPlane;   // The plane being drawn
    
    // Keep track of our stuff
    TEveViewer*            fXYView;
    TEveProjectionManager* fXYMgr;
    TEveScene*             fDetXYScene;
    TEveScene*             fEvtXYScene;
    
    TEveEventManager*      fEventManager;
    TEveFrameBox*          fFrameBox;
    TEveRGBAPalette*       fPalette;
    mutable TEveQuadSet*   fQuadSet;
    
    const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>();
};
    
//----------------------------------------------------------------------
// Constructor.
DrawWire2D::DrawWire2D(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fPlane = pset.get<int>("Plane", 2);

    // Create our own event manager
    fEventManager = new TEveEventManager("Event2D", "MicroBooNE Detector Event");
    
    gEve->AddEvent(fEventManager);

    // What follows is initially taken from a root tutorial (qaudset.C)
    gStyle->SetPalette(1, 0);
    
    fPalette  = new TEveRGBAPalette(-50, 200); // (0, 130);
    fFrameBox = new TEveFrameBox();
    
//    fFrameBox->SetAAQuadXY(-10, -10, 0, 20, 20); //-100., -10., 0., 3420., 6410.);
    fFrameBox->SetAAQuadXY(0., 0., 0., 3400., 6400.);
    fFrameBox->SetFrameColor(kGreen);

    fQuadSet = new TEveQuadSet(TEveQuadSet::kQT_RectangleXY, kFALSE, 1000, "RectangleXY");
    
    fQuadSet->SetOwnIds(kTRUE);
    fQuadSet->SetPalette(fPalette);
    fQuadSet->SetFrame(fFrameBox);
    
    // Set up a new tab and viewer for the 2D drawing
    // Create detector and event scenes for the 2D drawing
    fDetXYScene = gEve->SpawnNewScene("Detector 2D Scene", "");
    fEvtXYScene = gEve->SpawnNewScene("Detector 2D Scene", "");
    
//    gEve->AddElement(fEventManager, fDetXYScene);
    
    // Create XY projection mgr, draw projected axes, & add them to scenes
    fXYMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
    
    gEve->AddToListTree(fXYMgr, kFALSE);
    
    TEveProjectionAxes* axes_xy = new TEveProjectionAxes(fXYMgr);
    axes_xy->SetMainColor(kWhite);
    fDetXYScene->AddElement(axes_xy);
    
    // Create XY view in new tab & add det/evt scenes
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TEveWindowSlot* slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    TEveWindowPack* pack = slot->MakePack();
    
    pack->SetElementName("Wire Views");
    pack->SetHorizontal();
    pack->SetShowTitleBar(kTRUE);
    pack->NewSlot()->MakeCurrent();
    
    fXYView = gEve->SpawnNewViewer("Wire View", "");
    fXYView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
    fXYView->AddScene(fDetXYScene);
    fXYView->AddScene(fEvtXYScene);
    
    gEve->GetBrowser()->GetTabRight()->SetTab(0);

    TGLCameraOverlay* co = fXYView->GetGLViewer()->GetCameraOverlay();
    
    co->SetShowOrthographic(kTRUE);
    co->SetOrthographicMode(TGLCameraOverlay::kGridFront);
    
    // Uncomment these two lines to get internal highlight / selection.
    // q->SetPickable(1);
    // q->SetAlwaysSecSelect(1);
    
    TEveRGBAPaletteOverlay *po = new TEveRGBAPaletteOverlay(fPalette, 0.55, 0.1, 0.4, 0.05);
//    TEveRGBAPaletteOverlay *po = new TEveRGBAPaletteOverlay(fPalette, 100, 100, 400, 200);
    
    fXYView->GetGLViewer()->AddOverlayElement(po);
    
    
    
    gEve->GetDefaultGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    
    TGLCameraOverlay* coEve = gEve->GetDefaultGLViewer()->GetCameraOverlay();
    
    coEve->SetShowOrthographic(kTRUE);
    coEve->SetOrthographicMode(TGLCameraOverlay::kGridFront);
    
    // Uncomment these two lines to get internal highlight / selection.
    // q->SetPickable(1);
    // q->SetAlwaysSecSelect(1);
    
    TEveRGBAPaletteOverlay *poEve = new TEveRGBAPaletteOverlay(fPalette, 0.55, 0.1, 0.4, 0.05);
//    TEveRGBAPaletteOverlay *poEve = new TEveRGBAPaletteOverlay(fPalette, 100, 100, 400, 200);
    
    gEve->GetDefaultGLViewer()->AddOverlayElement(poEve);

    // To set user-interface (GUI + overlay) to display real values
    // mapped with a linear function: r = 0.1 * i + 0;
    // pal->SetUIDoubleRep(kTRUE, 0.1, 0);
    
    // Register
//    gEve->AddElement(fQuadSet);
//    gEve->Redraw3D(kTRUE);
    
    return;
}
    
DrawWire2D::~DrawWire2D()
{
}
    
void DrawWire2D::drawWire2D(const art::Event& event) const
{
    // Start by recovering the data
    art::Handle< std::vector<recob::Wire> > wireVecHandle;
    event.getByLabel("caldata",wireVecHandle);
    
//    fEvtXYScene->DestroyElements();
//    fEventManager->DestroyElements();
    
    if (!(wireVecHandle->size() > 0))
    {
        // Reset our quad object
        fQuadSet->Reset(TEveQuadSet::kQT_RectangleXY, kFALSE, wireVecHandle->size());
        
        for(const auto& wire : *wireVecHandle)
        {
            // Recover the channel id and use to get plane/view and wire #
            raw::ChannelID_t channel = wire.Channel();
        
            // get the WireID for this hit
            std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);
        
            // for now, just take the first option returned from ChannelToWire
            geo::WireID wid  = wids[0];
        
            // We need to know the plane to look up parameters
            geo::PlaneID::PlaneID_t plane   = wid.Plane;
            geo::WireID::WireID_t   wireNum = wid.Wire;
        
            // Ok, for now collection plane only...
            if (plane != 2) continue;

            // Recover the ROI's
            const recob::Wire::RegionsOfInterest_t& signalROI = wire.SignalROI();

            // and loop over them
            for(const auto& range : signalROI.get_ranges())
            {
                // Recover the actual signal vector
                const std::vector<float>& signal = range.data();
            
                // ROI start time
                raw::TDCtick_t roiFirstBinTick = range.begin_index();
            
                // Now loop through the ticks and draw the box for them
                for(size_t tick = 0; tick < signal.size(); tick++)
                {
                    int adcVal = signal.at(tick);
                    
                    adcVal = std::min(80,std::max(-20,adcVal));
                
                    fQuadSet->AddQuad(wireNum, roiFirstBinTick + tick, 0., 1., 1.);
                    fQuadSet->QuadValue(adcVal);
                    
//                    std::string quadId = "w_" + std::to_string(wireNum) + "_t_" + std::to_string(roiFirstBinTick+tick);
                    
//                    fQuadSet->QuadId(new TNamed(quadId.c_str(),"assigned name"));
                }
            }
        }
    }
    
    // Try the RawDigits
    art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
    event.getByLabel("digitfilter",rawDigitHandle);
    
    //    fEvtXYScene->DestroyElements();
    //    fEventManager->DestroyElements();
    
    if (rawDigitHandle->size() > 0)
    {
        //get pedestal conditions
        const lariov::DetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider();
        
        // Reset our quad object
        fQuadSet->Reset(TEveQuadSet::kQT_RectangleXY, kFALSE, 6400*rawDigitHandle->size());
        
        for(size_t rdIter = 0; rdIter < rawDigitHandle->size(); ++rdIter)
        {
            // get the reference to the current raw::RawDigit
            art::Ptr<raw::RawDigit> digitVec(rawDigitHandle, rdIter);
            
            raw::ChannelID_t channel = digitVec->Channel();
            
            std::vector<geo::WireID> wids     = fGeometry->ChannelToWire(channel);
            geo::PlaneID::PlaneID_t  thePlane = wids[0].Plane;
            geo::WireID::WireID_t    wireNum  = wids[0].Wire;
            
            if (thePlane != 2) continue;

            size_t dataSize = digitVec->Samples();
                
            // vector holding uncompressed adc values
            std::vector<short> rawadc(dataSize);
            
            // uncompress the data
            raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
            
            // loop over all adc values and subtract the pedestal
            // When we have a pedestal database, can provide the digit timestamp as the third argument of GetPedestalMean
            float pedestal = pedestalRetrievalAlg.PedMean(channel);
            
            // Get the pedestal subtracted data, centered in the deconvolution vector
            std::vector<float> rawAdcLessPedVec(dataSize);
            
            std::transform(rawadc.begin(),rawadc.end(),rawAdcLessPedVec.begin(),[pedestal](const short& adc){return std::round(float(adc) - pedestal);});
            
            for(size_t tick = 0; tick < rawAdcLessPedVec.size(); tick++)
            {
                int adcVal = rawAdcLessPedVec.at(tick);
                
                adcVal = std::min(200,std::max(-50,adcVal));
                
                fQuadSet->AddQuad(wireNum, tick, 0., 1., 1.);
                fQuadSet->QuadValue(adcVal);
            }
        }
    }
    
    // What follows is initially taken from a root tutorial (qaudset.C)
////    TEveQuadSet* fQuadSet = new TEveQuadSet(TEveQuadSet::kQT_RectangleXY, kFALSE, 1000, "RectangleXY");
    
////    fQuadSet->SetOwnIds(kTRUE);
////    fQuadSet->SetPalette(fPalette);
////    fQuadSet->SetFrame(fFrameBox);
//    fQuadSet->Reset(TEveQuadSet::kQT_RectangleXY, kFALSE, 1000);
    
    // NOT FILLED HERE REALLY
//    TRandom r(0);
//    for(int idx = 0; idx < 1000; idx += 10)
//    {
//        fQuadSet->AddQuad(r.Uniform(-10, 9), r.Uniform(-10, 9), 0, r.Uniform(0.2, 1), r.Uniform(0.2, 1));
//        fQuadSet->QuadValue(r.Uniform(0, 130));
//
//        std::string quadId = "tickid_" + std::to_string(idx);
//
//        fQuadSet->QuadId(new TNamed(quadId.c_str(),"assigned name"));
//    }
    
    fQuadSet->RefitPlex();
    
    // what is this for?
//    TEveTrans& translation = fQuadSet->RefMainTrans();
    
//    translation.RotateLF(1, 3, 0.5*TMath::Pi());
//    translation.SetPos(0., 0., 0.);
    
//    gEve->GetDefaultGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
    
//    TGLCameraOverlay* co = gEve->GetDefaultGLViewer()->GetCameraOverlay();
    
//    co->SetShowOrthographic(kTRUE);
//    co->SetOrthographicMode(TGLCameraOverlay::kGridFront);
    
    // Uncomment these two lines to get internal highlight / selection.
    // q->SetPickable(1);
    // q->SetAlwaysSecSelect(1);
    
//    TEveRGBAPaletteOverlay *po = new TEveRGBAPaletteOverlay(palette, 0.55, 0.1, 0.4, 0.05);
    
//    gEve->GetDefaultGLViewer()->AddOverlayElement(po);
    
    gEve->AddElement(fQuadSet, fEventManager);
    
    gEve->SetCurrentEvent(fEventManager);
    
//    fXYMgr->ImportElements(fEventManager, fEvtXYScene);
//    fXYView->Redraw(kTRUE);
//    gEve->Redraw3D(kTRUE);
//    gEve->DoRedraw3D();
    
    return;
}

DEFINE_ART_CLASS_TOOL(DrawWire2D)
}

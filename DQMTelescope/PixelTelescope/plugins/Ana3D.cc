// -*- C++ -*-
//
// Package:    DQMTelescope/Ana3D
// Class:      Ana3D
//
// class Ana3D Ana3D.cc DQMTelescope/Ana3D/plugins/Ana3D.cc

// Description: [one line class summary]

// Implementation:
//     [Notes on implementation]

// Original Author:  Jeremy Andrea
//      Updated by:  Nikkie Deelen
//         Created:  03.08.2018

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"

#include <cstring>
#include <string> 
#include <TH2F.h>
#include <TTree.h>
#include <TString.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection ;

class Ana3D : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:

    explicit Ana3D ( const edm::ParameterSet& ) ;
    ~Ana3D ( ) ;

    static void fillDescriptions ( edm::ConfigurationDescriptions& descriptions ) ;
   
  private:

    virtual void beginJob ( ) override ;
    virtual void analyze ( const edm::Event&, const edm::EventSetup& ) override ;
    virtual void endJob ( ) override ;

    // ----------member data ---------------------------

    edm::EDGetTokenT<TrackCollection> tracksToken_ ;  //used to select what tracks to read from configuration file
      
    edm::Service<TFileService> fs ;
     
    std::vector<TH1F *> DQM_ClusterCharge;
    std::vector<TH1F *> DQM_ClusterSize_X   ;  
    std::vector<TH1F *> DQM_ClusterSize_Y   ; 
    std::vector<TH1F *> DQM_ClusterSize_XY ;
    std::vector<TH1F *> DQM_NumbOfCluster;
    std::vector<TH2F *> DQM_ClusterPosition ;
      
    std::vector<TH2F *> DQM_DigiPosition ;
    std::vector<TH1F *> DQM_NumbOfDigi ;

    TH1F * DQM_NumbOfCluster_Tot ;
    TH1F * DQM_ClusterCharge_Tot ;
      
    TH1F * DQM_NumbOfDigi_Tot ;
      
    edm::EDGetTokenT<edm::DetSetVector<PixelDigi> >         pixeldigiToken_ ;
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > pixelclusterToken_ ;
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelRecHit> >  pixelhitToken_ ;
      
    std::vector<uint32_t > list_of_modules ;
    std::map<int, int> modulesNbr_to_idx ;
    std::map<int , TString> detId_to_moduleName ;
    std::vector<TString> names_of_modules ;

    TTree* cluster3DTree ;

    // Declaration of leaves types
    Int_t      tree_runNumber ;
    Int_t      tree_lumiSection ;
    Int_t      tree_event ;
    Int_t      tree_detId ;
    Int_t      tree_cluster ;
    Double_t   tree_x ;
    Double_t   tree_y ;
    Double_t   tree_z ;
    TString    tree_modName;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Ana3D::Ana3D(const edm::ParameterSet& iConfig)
 :
  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))
{
  //now do what ever initialization is needed
  TFileDirectory subDir = fs->mkdir(  "clusters3D" );
  cluster3DTree = subDir.make<TTree>("cluster3DTree", "3D Cluster Tree");

  // Set branch addresses.
  cluster3DTree->Branch("runNumber",&tree_runNumber);
  cluster3DTree->Branch("lumiSection",&tree_lumiSection);
  cluster3DTree->Branch("event",&tree_event);
  cluster3DTree->Branch("detId",&tree_detId);
  cluster3DTree->Branch("modName",&tree_modName);
  cluster3DTree->Branch("cluster",&tree_cluster);
  cluster3DTree->Branch("x",&tree_x);
  cluster3DTree->Branch("y",&tree_y);
  cluster3DTree->Branch("z",&tree_z);

    pixeldigiToken_    = consumes<edm::DetSetVector<PixelDigi> >        (iConfig.getParameter<edm::InputTag>("PixelDigisLabel"))   ;
    pixelclusterToken_ = consumes<edmNew::DetSetVector<SiPixelCluster> >(iConfig.getParameter<edm::InputTag>("PixelClustersLabel"));
    pixelhitToken_     = consumes<edmNew::DetSetVector<SiPixelRecHit> > (iConfig.getParameter<edm::InputTag>("PixelHitsLabel"))    ;

    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344200196, "M3090") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344201220, "M3124") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344462340, "M3082") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344463364, "M3175") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344724484, "M3009") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344725508, "M3057") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344986628, "M3027") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344987652, "M3074") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(352588804, "M3192") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(352589828, "M3204") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(352850948, "M3226") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(352851972, "M3265") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(353113092, "M3023") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(353114116, "M3239") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(353375236, "M3164") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(353376260, "M3173") ) ;
}

Ana3D::~Ana3D()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
Ana3D::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   EventID myEvId = iEvent.id();
   
  // Get collection of digi
  edm::Handle<edm::DetSetVector<PixelDigi> > pixeldigis;
  iEvent.getByToken(pixeldigiToken_,pixeldigis  );
  
  // Get collection of cluster
  edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixelclusters;
  iEvent.getByToken(pixelclusterToken_,pixelclusters  );
 
  // Get collection or RecHits
  edm::Handle< edmNew::DetSetVector<SiPixelRecHit> > pixelhits;
  iEvent.getByToken(pixelhitToken_,pixelhits  );
  
  // Get the geometry of the tracker for converting the LocalPoint to a GlobalPoint
  edm::ESHandle<TrackerGeometry> tracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
  const TrackerGeometry *tkgeom = &(*tracker); 
  
  // // Choose the CPE Estimator that will be used to estimate the LocalPoint of the cluster
  edm::ESHandle<PixelClusterParameterEstimator> cpEstimator;
  iSetup.get<TkPixelCPERecord>().get("PixelCPEGeneric", cpEstimator);
  const PixelClusterParameterEstimator &cpe(*cpEstimator); 
  
  //---------------------------------
  // loop on clusters
  //---------------------------------

  // First loop on modules with a cluster
  for( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter=pixelclusters->begin(); DSViter!=pixelclusters->end();DSViter++   ) {

    edmNew::DetSet<SiPixelCluster>::const_iterator begin=DSViter->begin();
    edmNew::DetSet<SiPixelCluster>::const_iterator end  =DSViter->end();

    // Surface of (detId) and use the surface to convert 2D to 3D      
    auto detId = DetId(DSViter->detId());
    int iCluster=0;
    const PixelGeomDetUnit *pixdet = (const PixelGeomDetUnit*) tkgeom->idToDetUnit(detId);
    LocalPoint lp(-9999., -9999., -9999.);

    // Then loop on the clusters of the module
    for(edmNew::DetSet<SiPixelCluster>::const_iterator itCluster=begin;itCluster!=end;++itCluster) {
      PixelClusterParameterEstimator::ReturnType params = cpe.getParameters(*itCluster,*pixdet);
      lp = std::get<0>(params);

      const Surface& surface = tracker->idToDet(detId)->surface();
      GlobalPoint gp = surface.toGlobal(lp);

      double x=0, y=0, z=0;
        x = gp.x();
        y = gp.y();
	z = gp.z();
	
        // Let's fill in the tree
        tree_runNumber = myEvId.run();
        tree_lumiSection = myEvId.luminosityBlock();
        tree_event = myEvId.event();
        tree_detId = detId;
        tree_cluster = iCluster++;	
        tree_x = x;
        tree_y = y;
	tree_z = z;
        tree_modName = detId_to_moduleName[detId];
        cluster3DTree->Fill();
      } //end for clusters of the first detector
    } //end for first detectors
}


// ------------ method called once each job just before starting event loop  ------------
void
Ana3D::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
Ana3D::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Ana3D::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ana3D);

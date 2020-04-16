#include "RecoEcal/EgammaClusterProducers/interface/PFECALSuperClusterProducer.h"

#include <memory>

#include "RecoEcal/EgammaClusterAlgos/interface/PFECALSuperClusterAlgo.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"

#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "TVector2.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;
using namespace edm;

namespace {
  const std::string ClusterType__BOX("Box");
  const std::string ClusterType__Mustache("Mustache");
  const std::string ClusterType__DeepSC("DeepSC");

  const std::string EnergyWeight__Raw("Raw");
  const std::string EnergyWeight__CalibratedNoPS("CalibratedNoPS");
  const std::string EnergyWeight__CalibratedTotal("CalibratedTotal");
}

PFECALSuperClusterProducer::PFECALSuperClusterProducer(const edm::ParameterSet& iConfig)
{
    
  verbose_ = 
    iConfig.getUntrackedParameter<bool>("verbose",false);

  superClusterAlgo_.setUseRegression(iConfig.getParameter<bool>("useRegression")); 

  isOOTCollection_ = iConfig.getParameter<bool>("isOOTCollection");
  superClusterAlgo_.setIsOOTCollection(isOOTCollection_);

  superClusterAlgo_.setTokens(iConfig,consumesCollector());

  std::string _typename = iConfig.getParameter<std::string>("ClusteringType");
  if( _typename == ClusterType__BOX ) {
    _theclusteringtype = PFECALSuperClusterAlgo::kBOX;
  } else if ( _typename == ClusterType__Mustache ) {
    _theclusteringtype = PFECALSuperClusterAlgo::kMustache;
  } else if ( _typename == ClusterType__DeepSC ) {
    _theclusteringtype = PFECALSuperClusterAlgo::kDeepSC;
  } else {
    throw cms::Exception("InvalidClusteringType") 
      << "You have not chosen a valid clustering type," 
      << " please choose from \"Box\" or \"Mustache\"!";
  }

  std::string _weightname = iConfig.getParameter<std::string>("EnergyWeight");
  if( _weightname == EnergyWeight__Raw ) {
    _theenergyweight = PFECALSuperClusterAlgo::kRaw;
  } else if ( _weightname == EnergyWeight__CalibratedNoPS ) {
    _theenergyweight = PFECALSuperClusterAlgo::kCalibratedNoPS;
  } else if ( _weightname == EnergyWeight__CalibratedTotal) {
    _theenergyweight = PFECALSuperClusterAlgo::kCalibratedTotal;
  } else {
    throw cms::Exception("InvalidClusteringType") 
      << "You have not chosen a valid energy weighting scheme," 
      << " please choose from \"Raw\", \"CalibratedNoPS\", or"
      << " \"CalibratedTotal\"!";
  }
  

  // parameters for clustering
  bool seedThresholdIsET = iConfig.getParameter<bool>("seedThresholdIsET");

  bool useDynamicDPhi = iConfig.getParameter<bool>("useDynamicDPhiWindow");

  double threshPFClusterSeedBarrel = iConfig.getParameter<double>("thresh_PFClusterSeedBarrel");
  double threshPFClusterBarrel = iConfig.getParameter<double>("thresh_PFClusterBarrel");

  double threshPFClusterSeedEndcap = iConfig.getParameter<double>("thresh_PFClusterSeedEndcap");
  double threshPFClusterEndcap = iConfig.getParameter<double>("thresh_PFClusterEndcap");
  
  double phiwidthSuperClusterBarrel = iConfig.getParameter<double>("phiwidth_SuperClusterBarrel");
  double etawidthSuperClusterBarrel = iConfig.getParameter<double>("etawidth_SuperClusterBarrel");

  double phiwidthSuperClusterEndcap = iConfig.getParameter<double>("phiwidth_SuperClusterEndcap");
  double etawidthSuperClusterEndcap = iConfig.getParameter<double>("etawidth_SuperClusterEndcap");

  //double threshPFClusterMustacheOutBarrel = iConfig.getParameter<double>("thresh_PFClusterMustacheOutBarrel");
  //double threshPFClusterMustacheOutEndcap = iConfig.getParameter<double>("thresh_PFClusterMustacheOutEndcap");

  double doSatelliteClusterMerge = 
    iConfig.getParameter<bool>("doSatelliteClusterMerge");
  double satelliteClusterSeedThreshold = 
    iConfig.getParameter<double>("satelliteClusterSeedThreshold");
  double satelliteMajorityFraction = 
    iConfig.getParameter<double>("satelliteMajorityFraction");
  bool dropUnseedable = iConfig.getParameter<bool>("dropUnseedable");

  if(_typename == ClusterType__DeepSC){
     NNfileName_ = iConfig.getUntrackedParameter<std::string>("NNfileName","deepSC_model.pb");
     NNevaluationParams_ = iConfig.getUntrackedParameter<std::vector<std::string>>("NNevaluationParams",std::vector<std::string>({std::string("dense_21_input:0"),std::string("dense_24/Sigmoid:0")}));
     NNnormalizationParams_mean_ = iConfig.getUntrackedParameter<std::vector<double>>("NNnormalizationParams_mean",std::vector<double>({0.,0.,0.,0.,0.,0.,0.}));
     NNnormalizationParams_std_ = iConfig.getUntrackedParameter<std::vector<double>>("NNnormalizationParams_std",std::vector<double>({1.,1.,1.,1.,1.,1.,1.})); 
     NNscoreWP_ = iConfig.getUntrackedParameter<double>("NNscoreWP",-1.);
  }

  // Mustache parameters (dEta parabolas deparation)
  double p00 = iConfig.getParameter<double>("p00");
  double p01 = iConfig.getParameter<double>("p01");
  double p02 = iConfig.getParameter<double>("p02");
  double p10 = iConfig.getParameter<double>("p10");
  double p11 = iConfig.getParameter<double>("p11");
  double p12 = iConfig.getParameter<double>("p12");
  double w00 = iConfig.getParameter<double>("w00");
  double w01 = iConfig.getParameter<double>("w01");
  double w10 = iConfig.getParameter<double>("w10");
  double w11 = iConfig.getParameter<double>("w11");
  // dynamic dPhi window 
  double yoffsetEB   = iConfig.getParameter<double>("yoffsetEB");
  double scaleEB     = iConfig.getParameter<double>("scaleEB");
  double xoffsetEB   = iConfig.getParameter<double>("xoffsetEB");
  double widthEB     = iConfig.getParameter<double>("widthEB");
  double yoffsetEE_0 = iConfig.getParameter<double>("yoffsetEE_0");
  double scaleEE_0   = iConfig.getParameter<double>("scaleEE_0");
  double xoffsetEE_0 = iConfig.getParameter<double>("xoffsetEE_0");
  double widthEE_0   = iConfig.getParameter<double>("widthEE_0");
  double yoffsetEE_1 = iConfig.getParameter<double>("yoffsetEE_1");
  double scaleEE_1   = iConfig.getParameter<double>("scaleEE_1");
  double xoffsetEE_1 = iConfig.getParameter<double>("xoffsetEE_1");
  double widthEE_1   = iConfig.getParameter<double>("widthEE_1");
  double yoffsetEE_2 = iConfig.getParameter<double>("yoffsetEE_2");
  double scaleEE_2   = iConfig.getParameter<double>("scaleEE_2");
  double xoffsetEE_2 = iConfig.getParameter<double>("xoffsetEE_2");
  double widthEE_2   = iConfig.getParameter<double>("widthEE_2");



  mustacheParams_.resize(10);
  dynPhiWindowParams_.resize(16);

  mustacheParams_[0]= p00;
  mustacheParams_[1]= p01;
  mustacheParams_[2]= p02;
  mustacheParams_[3]= p10;
  mustacheParams_[4]= p11;
  mustacheParams_[5]= p12;
  mustacheParams_[6]= w00;
  mustacheParams_[7]= w01;
  mustacheParams_[8]= w10;
  mustacheParams_[9]= w11;


  dynPhiWindowParams_[0] = yoffsetEB;
  dynPhiWindowParams_[1] = scaleEB;
  dynPhiWindowParams_[2] = xoffsetEB;
  dynPhiWindowParams_[3] = widthEB;
  dynPhiWindowParams_[4] = yoffsetEE_0;
  dynPhiWindowParams_[5] = scaleEE_0;
  dynPhiWindowParams_[6] = xoffsetEE_0;
  dynPhiWindowParams_[7] = widthEE_0;
  dynPhiWindowParams_[8] = yoffsetEE_1;
  dynPhiWindowParams_[9] =  scaleEE_1;
  dynPhiWindowParams_[10] = xoffsetEE_1;
  dynPhiWindowParams_[11] = widthEE_1;
  dynPhiWindowParams_[12] = yoffsetEE_2;
  dynPhiWindowParams_[13] = scaleEE_2;
  dynPhiWindowParams_[14] = xoffsetEE_2;
  dynPhiWindowParams_[15] =  widthEE_2;

  //std::cout << " In the producer dynPhiWindowParams_ size " << dynPhiWindowParams_.size() << std::endl;

  superClusterAlgo_.setVerbosityLevel(verbose_);
  superClusterAlgo_.setClusteringType(_theclusteringtype);
  superClusterAlgo_.setEnergyWeighting(_theenergyweight);
  superClusterAlgo_.setUseETForSeeding(seedThresholdIsET);
  superClusterAlgo_.setUseDynamicDPhi(useDynamicDPhi);

  superClusterAlgo_.setThreshSuperClusterEt( iConfig.getParameter<double>("thresh_SCEt") );
  
  superClusterAlgo_.setThreshPFClusterSeedBarrel( threshPFClusterSeedBarrel );
  superClusterAlgo_.setThreshPFClusterBarrel( threshPFClusterBarrel );

  superClusterAlgo_.setThreshPFClusterSeedEndcap( threshPFClusterSeedEndcap );
  superClusterAlgo_.setThreshPFClusterEndcap( threshPFClusterEndcap );

  superClusterAlgo_.setPhiwidthSuperClusterBarrel( phiwidthSuperClusterBarrel );
  superClusterAlgo_.setEtawidthSuperClusterBarrel( etawidthSuperClusterBarrel );

  superClusterAlgo_.setPhiwidthSuperClusterEndcap( phiwidthSuperClusterEndcap );
  superClusterAlgo_.setEtawidthSuperClusterEndcap( etawidthSuperClusterEndcap );

  superClusterAlgo_.setSatelliteMerging( doSatelliteClusterMerge );
  superClusterAlgo_.setSatelliteThreshold( satelliteClusterSeedThreshold );
  superClusterAlgo_.setMajorityFraction( satelliteMajorityFraction );
  superClusterAlgo_.setDropUnseedable( dropUnseedable );
  //superClusterAlgo_.setThreshPFClusterMustacheOutBarrel( threshPFClusterMustacheOutBarrel );
  //superClusterAlgo_.setThreshPFClusterMustacheOutEndcap( threshPFClusterMustacheOutEndcap );
  superClusterAlgo_.setMustacheParameters(mustacheParams_);
  superClusterAlgo_.setDynPhiWindowParameters(dynPhiWindowParams_);

  if(_typename == ClusterType__DeepSC){
     superClusterAlgo_.setNNSession( NNfileName_ );
     superClusterAlgo_.setNNEvaluationParams(NNevaluationParams_);
     superClusterAlgo_.setNNnormalizationParams( NNnormalizationParams_mean_, NNnormalizationParams_std_ );
     superClusterAlgo_.setNNscoreWP( NNscoreWP_ );
  }
 
  //Load the ECAL energy calibration
  thePFEnergyCalibration_ = 
    std::make_shared<PFEnergyCalibration>();
  superClusterAlgo_.setPFClusterCalibration(thePFEnergyCalibration_);

  bool applyCrackCorrections_ = iConfig.getParameter<bool>("applyCrackCorrections");
  superClusterAlgo_.setCrackCorrections(applyCrackCorrections_);
  
  PFBasicClusterCollectionBarrel_ = iConfig.getParameter<string>("PFBasicClusterCollectionBarrel");
  PFSuperClusterCollectionBarrel_ = iConfig.getParameter<string>("PFSuperClusterCollectionBarrel");

  PFBasicClusterCollectionEndcap_ = iConfig.getParameter<string>("PFBasicClusterCollectionEndcap");
  PFSuperClusterCollectionEndcap_ = iConfig.getParameter<string>("PFSuperClusterCollectionEndcap");

  PFBasicClusterCollectionPreshower_ = iConfig.getParameter<string>("PFBasicClusterCollectionPreshower");
  PFSuperClusterCollectionEndcapWithPreshower_ = iConfig.getParameter<string>("PFSuperClusterCollectionEndcapWithPreshower");

  PFClusterAssociationEBEE_ = "PFClusterAssociationEBEE";
  PFClusterAssociationES_ = "PFClusterAssociationES";
  
  produces<reco::SuperClusterCollection>(PFSuperClusterCollectionBarrel_);  
  produces<reco::SuperClusterCollection>(PFSuperClusterCollectionEndcapWithPreshower_);
  produces<reco::CaloClusterCollection>(PFBasicClusterCollectionBarrel_);
  produces<reco::CaloClusterCollection>(PFBasicClusterCollectionEndcap_);
  produces<reco::CaloClusterCollection>(PFBasicClusterCollectionPreshower_);
  produces<edm::ValueMap<reco::CaloClusterPtr> >(PFClusterAssociationEBEE_);
  produces<edm::ValueMap<reco::CaloClusterPtr> >(PFClusterAssociationES_);
}



PFECALSuperClusterProducer::~PFECALSuperClusterProducer() {}

void PFECALSuperClusterProducer::
beginLuminosityBlock(LuminosityBlock const& iL, EventSetup const& iE) {
  superClusterAlgo_.update(iE); 
}


void PFECALSuperClusterProducer::produce(edm::Event& iEvent, 
				const edm::EventSetup& iSetup) {
  
  // do clustering
  superClusterAlgo_.loadAndSortPFClusters(iEvent, iSetup);
  superClusterAlgo_.run();
  
  //build collections of output CaloClusters from the used PFClusters
  auto caloClustersEB = std::make_unique<reco::CaloClusterCollection>();
  auto caloClustersEE = std::make_unique<reco::CaloClusterCollection>();
  auto caloClustersES = std::make_unique<reco::CaloClusterCollection>();
  
  std::map<reco::CaloClusterPtr, unsigned int> pfClusterMapEB; //maps of pfclusters to caloclusters 
  std::map<reco::CaloClusterPtr, unsigned int> pfClusterMapEE;
  std::map<reco::CaloClusterPtr, unsigned int> pfClusterMapES;
  
  //fill calocluster collections and maps
  for( const auto& ebsc : *(superClusterAlgo_.getEBOutputSCCollection()) ) {
    for (reco::CaloCluster_iterator pfclus = ebsc.clustersBegin(); pfclus!=ebsc.clustersEnd(); ++pfclus) {
      if (!pfClusterMapEB.count(*pfclus)) {
        reco::CaloCluster caloclus(**pfclus);
        caloClustersEB->push_back(caloclus);
        pfClusterMapEB[*pfclus] = caloClustersEB->size() - 1;
      }
      else {
        throw cms::Exception("PFECALSuperClusterProducer::produce")
            << "Found an EB pfcluster matched to more than one EB supercluster!" 
            << std::dec << std::endl;
      }
    }
  }
  for( const auto& eesc : *(superClusterAlgo_.getEEOutputSCCollection()) ) {
    for (reco::CaloCluster_iterator pfclus = eesc.clustersBegin(); pfclus!=eesc.clustersEnd(); ++pfclus) {
      if (!pfClusterMapEE.count(*pfclus)) {
        reco::CaloCluster caloclus(**pfclus);
        caloClustersEE->push_back(caloclus);
        pfClusterMapEE[*pfclus] = caloClustersEE->size() - 1;
      }
      else {
        throw cms::Exception("PFECALSuperClusterProducer::produce")
            << "Found an EE pfcluster matched to more than one EE supercluster!" 
            << std::dec << std::endl;
      }
    }
    for (reco::CaloCluster_iterator pfclus = eesc.preshowerClustersBegin(); pfclus!=eesc.preshowerClustersEnd(); ++pfclus) {
      if (!pfClusterMapES.count(*pfclus)) {
        reco::CaloCluster caloclus(**pfclus);
        caloClustersES->push_back(caloclus);
        pfClusterMapES[*pfclus] = caloClustersES->size() - 1;
      }
      else {
        throw cms::Exception("PFECALSuperClusterProducer::produce")
            << "Found an ES pfcluster matched to more than one EE supercluster!" 
            << std::dec << std::endl;
      }
    }
  }
  
  //create ValueMaps from output CaloClusters back to original PFClusters
  auto pfClusterAssociationEBEE = std::make_unique<edm::ValueMap<reco::CaloClusterPtr>>();
  auto pfClusterAssociationES = std::make_unique<edm::ValueMap<reco::CaloClusterPtr>>();    
  
  //vectors to fill ValueMaps
  std::vector<reco::CaloClusterPtr> clusptrsEB(caloClustersEB->size());
  std::vector<reco::CaloClusterPtr> clusptrsEE(caloClustersEE->size());
  std::vector<reco::CaloClusterPtr> clusptrsES(caloClustersES->size());  
  
  //put calocluster output collections in event and get orphan handles to create ptrs
  const edm::OrphanHandle<reco::CaloClusterCollection> &caloClusHandleEB = iEvent.put(std::move(caloClustersEB),PFBasicClusterCollectionBarrel_);
  const edm::OrphanHandle<reco::CaloClusterCollection> &caloClusHandleEE = iEvent.put(std::move(caloClustersEE),PFBasicClusterCollectionEndcap_);
  const edm::OrphanHandle<reco::CaloClusterCollection> &caloClusHandleES = iEvent.put(std::move(caloClustersES),PFBasicClusterCollectionPreshower_);
    
  //relink superclusters to output caloclusters and fill vectors for ValueMaps
  for( auto& ebsc : *(superClusterAlgo_.getEBOutputSCCollection()) ) {
    reco::CaloClusterPtr seedptr(caloClusHandleEB,pfClusterMapEB[ebsc.seed()]);
    ebsc.setSeed(seedptr);
    
    reco::CaloClusterPtrVector clusters;
    for (reco::CaloCluster_iterator pfclus = ebsc.clustersBegin(); pfclus!=ebsc.clustersEnd(); ++pfclus) {
      int caloclusidx = pfClusterMapEB[*pfclus];
      reco::CaloClusterPtr clusptr(caloClusHandleEB,caloclusidx);
      clusters.push_back(clusptr);
      clusptrsEB[caloclusidx] = *pfclus;
    }
    ebsc.setClusters(clusters);
  }
  for( auto& eesc : *(superClusterAlgo_.getEEOutputSCCollection()) ) {
    reco::CaloClusterPtr seedptr(caloClusHandleEE,pfClusterMapEE[eesc.seed()]);
    eesc.setSeed(seedptr);
    
    reco::CaloClusterPtrVector clusters;
    for (reco::CaloCluster_iterator pfclus = eesc.clustersBegin(); pfclus!=eesc.clustersEnd(); ++pfclus) {
      int caloclusidx = pfClusterMapEE[*pfclus];
      reco::CaloClusterPtr clusptr(caloClusHandleEE,caloclusidx);
      clusters.push_back(clusptr);
      clusptrsEE[caloclusidx] = *pfclus;
    }
    eesc.setClusters(clusters);
    
    reco::CaloClusterPtrVector psclusters;
    for (reco::CaloCluster_iterator pfclus = eesc.preshowerClustersBegin(); pfclus!=eesc.preshowerClustersEnd(); ++pfclus) {
      int caloclusidx = pfClusterMapES[*pfclus];
      reco::CaloClusterPtr clusptr(caloClusHandleES,caloclusidx);
      psclusters.push_back(clusptr);
      clusptrsES[caloclusidx] = *pfclus;
    }
    eesc.setPreshowerClusters(psclusters);  
  }
    
  //fill association maps from output CaloClusters back to original PFClusters
  edm::ValueMap<reco::CaloClusterPtr>::Filler fillerEBEE(*pfClusterAssociationEBEE);
  fillerEBEE.insert(caloClusHandleEB, clusptrsEB.begin(), clusptrsEB.end());
  fillerEBEE.insert(caloClusHandleEE, clusptrsEE.begin(), clusptrsEE.end());
  fillerEBEE.fill();
  
  edm::ValueMap<reco::CaloClusterPtr>::Filler fillerES(*pfClusterAssociationES);
  fillerES.insert(caloClusHandleES, clusptrsES.begin(), clusptrsES.end());  
  fillerES.fill();  

  //store in the event
  iEvent.put(std::move(pfClusterAssociationEBEE),PFClusterAssociationEBEE_);
  iEvent.put(std::move(pfClusterAssociationES),PFClusterAssociationES_);
  iEvent.put(std::move(superClusterAlgo_.getEBOutputSCCollection()),
	     PFSuperClusterCollectionBarrel_);
  iEvent.put(std::move(superClusterAlgo_.getEEOutputSCCollection()), 
	     PFSuperClusterCollectionEndcapWithPreshower_);
}

void PFECALSuperClusterProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("PFSuperClusterCollectionEndcap","particleFlowSuperClusterECALEndcap");
  desc.add<bool>("doSatelliteClusterMerge",false);
  desc.add<double>("thresh_PFClusterBarrel",0.0);
  desc.add<std::string>("PFBasicClusterCollectionBarrel","particleFlowBasicClusterECALBarrel");
  desc.add<bool>("useRegression",true);
  desc.add<double>("satelliteMajorityFraction",0.5);
  desc.add<double>("thresh_PFClusterEndcap",0.0);
  desc.add<edm::InputTag>("ESAssociation",edm::InputTag("particleFlowClusterECAL"));
  desc.add<std::string>("PFBasicClusterCollectionPreshower","particleFlowBasicClusterECALPreshower");
  desc.add<bool>("use_preshower",true);
  desc.addUntracked<bool>("verbose",false);
  desc.add<double>("thresh_SCEt",4.0);
  desc.add<double>("etawidth_SuperClusterEndcap",0.04);
  desc.add<double>("phiwidth_SuperClusterEndcap",0.6);
  desc.add<bool>("useDynamicDPhiWindow",true);
  desc.add<std::string>("PFSuperClusterCollectionBarrel","particleFlowSuperClusterECALBarrel");
  {
    edm::ParameterSetDescription psd0;
    psd0.add<bool>("isHLT", false);
    psd0.add<bool>("applySigmaIetaIphiBug", false);
    psd0.add<edm::InputTag>("ecalRecHitsEE",edm::InputTag("ecalRecHit","EcalRecHitsEE"));
    psd0.add<edm::InputTag>("ecalRecHitsEB",edm::InputTag("ecalRecHit","EcalRecHitsEB"));
    psd0.add<std::string>("regressionKeyEB","pfscecal_EBCorrection_offline_v1");
    psd0.add<std::string>("regressionKeyEE","pfscecal_EECorrection_offline_v1");
    psd0.add<std::string>("uncertaintyKeyEB","pfscecal_EBUncertainty_offline_v1");
    psd0.add<std::string>("uncertaintyKeyEE","pfscecal_EEUncertainty_offline_v1");
    psd0.add<edm::InputTag>("vertexCollection",edm::InputTag("offlinePrimaryVertices")); 
    psd0.add<double>("eRecHitThreshold", 1.);
    desc.add<edm::ParameterSetDescription>("regressionConfig",psd0);
  }
  desc.add<bool>("applyCrackCorrections",false);
  desc.add<double>("satelliteClusterSeedThreshold",50.0);
  desc.add<double>("etawidth_SuperClusterBarrel",0.04);
  desc.add<std::string>("PFBasicClusterCollectionEndcap","particleFlowBasicClusterECALEndcap");
  desc.add<edm::InputTag>("PFClusters",edm::InputTag("particleFlowClusterECAL"));
  desc.add<double>("thresh_PFClusterSeedBarrel",1.0);
  desc.add<std::string>("ClusteringType","Mustache");
  desc.add<std::string>("EnergyWeight","Raw");
  desc.add<edm::InputTag>("BeamSpot",edm::InputTag("offlineBeamSpot"));
  desc.add<double>("thresh_PFClusterSeedEndcap",1.0);
  desc.add<double>("phiwidth_SuperClusterBarrel",0.6);
  desc.add<double>("thresh_PFClusterES",0.0);
  desc.add<bool>("seedThresholdIsET",true);
  desc.add<bool>("isOOTCollection",false);
  desc.add<edm::InputTag>("barrelRecHits",edm::InputTag("ecalRecHit","EcalRecHitsEB"));
  desc.add<edm::InputTag>("endcapRecHits",edm::InputTag("ecalRecHit","EcalRecHitsEE"));
  desc.add<std::string>("PFSuperClusterCollectionEndcapWithPreshower","particleFlowSuperClusterECALEndcapWithPreshower");
  desc.add<bool>("dropUnseedable",false);
  // old mustache parameters, before optimization 2019-2020
  desc.add<double>("p00",-0.107537);
  desc.add<double>("p01",0.590969);
  desc.add<double>("p02",-0.076494);
  desc.add<double>("p10",-0.0268843);
  desc.add<double>("p11",0.147742);
  desc.add<double>("p12",-0.0191235);
  desc.add<double>("w00",-0.00571429);
  desc.add<double>("w01",-0.002);
  desc.add<double>("w10", 0.0135714);
  desc.add<double>("w11",0.001);
  desc.add<double>("yoffsetEB",7.151e-02);
  desc.add<double>("scaleEB",5.656e-01);
  desc.add<double>("xoffsetEB",2.931e-01);
  desc.add<double>("widthEB",2.976e-01);
  desc.add<double>("yoffsetEE_0",5.058e-02);
  desc.add<double>("scaleEE_0",7.131e-01);
  desc.add<double>("xoffsetEE_0",1.668e-02);
  desc.add<double>("widthEE_0",4.114e-01);
  desc.add<double>("yoffsetEE_1",-9.913e-02);
  desc.add<double>("scaleEE_1",4.404e+01);
  desc.add<double>("xoffsetEE_1",-5.326e+00);
  desc.add<double>("widthEE_1",1.184e+00);
  desc.add<double>("yoffsetEE_2",-6.346e-01);
  desc.add<double>("scaleEE_2",1.317e+01);
  desc.add<double>("xoffsetEE_2",-7.037e+00);
  desc.add<double>("widthEE_2",2.836e+00);
  desc.addUntracked<std::string>("NNfileName","deepSC_model.pb");
  desc.addUntracked<std::vector<std::string>>("NNevaluationParams",std::vector<std::string>({std::string("dense_21_input:0"),std::string("dense_24/Sigmoid:0")}));
  desc.addUntracked<std::vector<double>>("NNnormalizationParams_mean",std::vector<double>({0.,0.,0.,0.,0.,0.,0.}));
  desc.addUntracked<std::vector<double>>("NNnormalizationParams_std",std::vector<double>({1.,1.,1.,1.,1.,1.,1.}));
  desc.addUntracked<double>("NNscoreWP",-1.); 
  descriptions.add("particleFlowSuperClusterECALMustache",desc);
}

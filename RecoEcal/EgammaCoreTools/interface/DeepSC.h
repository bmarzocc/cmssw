#ifndef RecoEcal_EgammaCoreTools_DeepSC_h
#define RecoEcal_EgammaCoreTools_DeepSC_h

#include <vector>
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

namespace reco {
  class DeepSC {
    
  public:

    explicit DeepSC();
    ~DeepSC();
   
    double DeltaPhi(double seed_phi, double cluster_phi);
    double DeltaEta(double seed_eta, double cluster_eta); 
    std::pair<double,double> ComputeCovariances(const CaloCluster cluster, const EcalRecHitCollection* recHits, const CaloSubdetectorGeometry* geometry);
    void ComputeVariables(const CaloCluster* seed, const CaloCluster* cluster, const CaloTopology *topology, const CaloSubdetectorGeometry* ebGeom, const CaloSubdetectorGeometry* eeGeom, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE);
    void SetNNVarVal(std::vector<double> vars);
    void NormalizeNNVars(const CaloCluster* seed);
    float EvaluateNN(const CaloCluster* seed);
    bool InSuperCluster(const CaloCluster* seed, const CaloCluster* cluster, const CaloTopology *topology, const CaloSubdetectorGeometry* ebGeom, const CaloSubdetectorGeometry* eeGeom, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE);
    bool InSuperCluster(std::vector<double> clusterVariables_);
    std::vector<double> NNclusterVars(){ return NNclusterVars_; }; 
    void SetXtalsInWindow(DetId seedDetId, double etawidthSuperCluster, double phiwidthSuperCluster, const CaloTopology *topology, const CaloGeometry *geometry);
    void ClearXtalsInWindow();
    std::vector<DetId> XtalsInWindow(){ return xtals_inWindow_; };
    std::pair<double,double> GetMaximumDetaDphi(DetId* seedDetId, std::vector<DetId>* idMatrix, const CaloGeometry *geometry);  
    void DeepSCID(const CaloClusterPtrVector& clusters, 
		  int & nclusters, float & EoutsideDeepSC,
                  const CaloTopology *topology, 
                  const CaloSubdetectorGeometry* ebGeom, 
                  const CaloSubdetectorGeometry* eeGeom,
                  const EcalRecHitCollection *recHitsEB, 
                  const EcalRecHitCollection *recHitsEE);
    void DeepSCID(const std::vector<const CaloCluster*>&, 
		  int & nclusers,
		  float & EoutsideDeepSC,
                  const CaloTopology *topology, 
                  const CaloSubdetectorGeometry* ebGeom, 
                  const CaloSubdetectorGeometry* eeGeom,
                  const EcalRecHitCollection *recHitsEB, 
                  const EcalRecHitCollection *recHitsEE); 
    void DeepSCID(const reco::SuperCluster& sc, 
		  int & nclusters, 
		  float & EoutsideDeepSC,
                  const CaloTopology *topology, 
                  const CaloSubdetectorGeometry* ebGeom, 
                  const CaloSubdetectorGeometry* eeGeom,
                  const EcalRecHitCollection *recHitsEB, 
                  const EcalRecHitCollection *recHitsEE);


    void DeepSCClust(const std::vector<CaloCluster>& clusters,
                     const CaloTopology *topology, 
                     const CaloSubdetectorGeometry* ebGeom, 
                     const CaloSubdetectorGeometry* eeGeom,
                     const EcalRecHitCollection *recHitsEB, 
                     const EcalRecHitCollection *recHitsEE,  
		     std::vector<unsigned int>& insideDeep, 
		     std::vector<unsigned int>& outsideDeep);
    
    void FillDeepSCVar(const std::vector<CaloCluster>& clusters,
                       const CaloTopology *topology, 
                       const CaloSubdetectorGeometry* ebGeom, 
                       const CaloSubdetectorGeometry* eeGeom,
                       const EcalRecHitCollection *recHitsEB, 
                       const EcalRecHitCollection *recHitsEE);
    //return Functions for DeepSC Variables:
    float DeepSCE(){return Energy_In_DeepSC_;}
    float DeepSCEOut(){return Energy_Outside_DeepSC_;}
    float DeepSCEtOut(){return Et_Outside_DeepSC_;}
    float LowestDeepClust(){return LowestClusterEInDeepSC_;}
    int InsideDeep(){return included_;}
    int OutsideDeep(){return excluded_;}

  private:

    GlobalPoint cell_;
    std::string NNinput_string_Barrel_;
    std::string NNinput_string_Endcap_;
    std::string NNoutput_string_Barrel_;
    std::string NNoutput_string_Endcap_;
    std::vector<double> x_mean_Barrel_;
    std::vector<double> x_mean_Endcap_;
    std::vector<double> x_std_Barrel_;
    std::vector<double> x_std_Endcap_;
    tensorflow::Session* session_NN_Barrel_;
    tensorflow::Session* session_NN_Endcap_;
    std::vector<double> NNWPs_Barrel_;
    std::vector<double> NNWPs_Endcap_;
    EcalClusterTools egmTools_;
    std::vector<DetId> xtals_inWindow_;
    std::vector<float> full5x5_locCov_;
    std::vector<float> showerShapes_; 
    std::pair<double,double> widths_;
    std::vector<double> NNclusterVars_;  
    std::vector<double> vectorVar_;
    std::vector<double> etaWindows_Barrel_; 
    std::vector<double> etaWindows_Endcap_;  
    std::vector<double> phiWindows_Barrel_; 
    std::vector<double> phiWindows_Endcap_; 
    float Energy_In_DeepSC_;
    float Energy_Outside_DeepSC_;
    float Et_Outside_DeepSC_;
    float LowestClusterEInDeepSC_;
    int excluded_;
    int included_; 
    template<class RandomAccessPtrIterator>
    void DeepSCID(const RandomAccessPtrIterator&,
		  const RandomAccessPtrIterator&,
		  int& nclusters,
		  float& EoutsideDeepSC,
                  const CaloTopology *topology, 
                  const CaloSubdetectorGeometry* ebGeom, 
                  const CaloSubdetectorGeometry* eeGeom, 
                  const EcalRecHitCollection *recHitsEB, 
                  const EcalRecHitCollection *recHitsEE);
    
  };

  
}

#endif

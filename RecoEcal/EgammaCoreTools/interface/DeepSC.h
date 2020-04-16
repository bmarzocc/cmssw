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

    void SetXtalsInWindow(DetId seedDetId, double etawidthSuperCluster, double phiwidthSuperCluster, const CaloTopology *topology, const CaloGeometry *geometry);
    void ClearXtalsInWindow();
    std::vector<DetId> XtalsInWindow(){ return xtals_inWindow_; };
    std::pair<double,double> GetMaximumDetaDphi(DetId* seedDetId, std::vector<DetId>* idMatrix, const CaloGeometry *geometry);
    std::vector<double> NNclusterVars(){ return NNclusterVars_; };
    
    void ClearNNVars();
    void ClearNNnormalizeParams();
    void SetNNGraphSession(std::string NNfileName_);
    void SetNNEvaluationParams(std::string NNinput, std::string NNoutput);
    void SetNNVarVal(std::vector<double> vars);
    void SetNormalizeNNParams(std::vector<double> var_mean_, std::vector<double> var_std_);
    void SetClusterVariables(const edm::Ptr<reco::PFCluster> cluster, const CaloTopology *topology, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE);
    void NormalizeNNVars();
    float EvaluateNN();

  private:

    GlobalPoint cell_;
    std::string NNinput_string_;
    std::string NNoutput_string_;
    std::vector<DetId> xtals_inWindow_;
    std::vector<double> vectorVar_;
    std::vector<double> x_mean_, x_std_;
    tensorflow::GraphDef* graphDef_NN_;
    tensorflow::Session* session_NN_;
    std::vector<float> full5x5_locCov_;
    std::vector<double> NNclusterVars_;
    EcalClusterTools egmTools_;
    
  };

  
}

#endif

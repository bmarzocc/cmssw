/**
   \file
   Tools for manipulating ECAL Clusters as graphs
   \author Davide Valsecchi, Badder Marzocchi
   \date 05 October 2020
*/

#ifndef RecoEcal_EgammaCoreTools_EcalClustersGraph_h
#define RecoEcal_EgammaCoreTools_EcalClustersGraph_h

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

#include "RecoEcal/EgammaClusterAlgos/interface/PFECALSuperClusterAlgo.h"
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace reco;

namespace ublas = boost::numeric::ublas;

class EcalClustersGraph {

  typedef std::shared_ptr<PFECALSuperClusterAlgo::CalibratedPFCluster> CalibratedClusterPtr;
  typedef std::vector<PFECALSuperClusterAlgo::CalibratedClusterPtr> CalibratedClusterPtrVector;

private:

   CalibratedClusterPtrVector clusters_;
   int nSeeds_;
   int nCls_;

   // Adjacency matrix defining which clusters are inside the seeds windows.
   // row: seeds (Et ordered), column: clusters (Et ordered)
   ublas::matrix<int> inWindows_;
   // Adjacency matrix defining how much each cluster is linked to the seed
   // row: seeds (Et ordered), column: clusters (Et ordered)
   ublas::matrix<double> scoreMatrix_;

   //To compute the input variables
   const CaloTopology* topology_;
   const CaloSubdetectorGeometry* ebGeom_;
   const CaloSubdetectorGeometry* eeGeom_;
   const EcalRecHitCollection* recHitsEB_;
   const EcalRecHitCollection* recHitsEE_;
   const std::vector<double>* meanVals_;
   const std::vector<double>* stdVals_;
   std::vector<float> full5x5_locCov_;
   std::pair<double,double> widths_;
   std::vector<double> showerShapes_;
   std::vector<double> NNclusterVars_;
   std::vector<std::vector<std::vector<double>>> NNwindowVars_;

public:

   EcalClustersGraph(CalibratedClusterPtrVector clusters, int nSeeds, const CaloTopology *topology, const CaloSubdetectorGeometry* ebGeom, const CaloSubdetectorGeometry* eeGeom, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE, const std::vector<double>* meanVals, const std::vector<double>* stdVals);

   double deltaPhi(double seed_phi, double cluster_phi);
   double deltaEta(double seed_eta, double cluster_eta); 
   std::pair<double,double> computeCovariances(const CaloCluster* cluster);
   std::vector<double> computeShowerShapes(const CaloCluster* cluster);
   void computeVariables(const CaloCluster* seed, const CaloCluster* cluster);
   void initWindows(double etaWidth, double phiWidth);
   void clearWindows();
   void fillVariables();

};

#endif

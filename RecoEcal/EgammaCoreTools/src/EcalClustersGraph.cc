#include "RecoEcal/EgammaCoreTools/interface/EcalClustersGraph.h"
#include <algorithm>
#include <cmath>
#include "TVector2.h"
#include "TMath.h"


using namespace std;
using namespace reco;

EcalClustersGraph::EcalClustersGraph(CalibratedClusterPtrVector clusters, int nSeeds, const CaloTopology *topology, const CaloSubdetectorGeometry* ebGeom, const CaloSubdetectorGeometry* eeGeom, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE, const std::vector<double>* meanVals, const std::vector<double>* stdVals):
    clusters_(clusters), nSeeds_(nSeeds), topology_(topology), ebGeom_(ebGeom), eeGeom_(eeGeom), recHitsEB_(recHitsEB), recHitsEE_(recHitsEE), meanVals_(meanVals), stdVals_(stdVals){
       nCls_ = clusters_.size();
       inWindows_ = ublas::matrix<int> (nSeeds_, nCls_);
       scoreMatrix_ = ublas::matrix<int> (nSeeds_, nCls_);

       //test
       std::cout << "ClustersGraph created. Nseeds " << nSeeds_ << ", nClusters " << nCls_ << endl;
    }

double EcalClustersGraph::deltaPhi(double seed_phi, double cluster_phi)
{
     double dphi = seed_phi - cluster_phi;
     if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
     if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
     return dphi;
} 

double EcalClustersGraph::deltaEta(double seed_eta, double cluster_eta)
{
     double deta = 0.;
     if(seed_eta > 0.) deta = cluster_eta - seed_eta;
     if(seed_eta <= 0.) deta = seed_eta - cluster_eta;
     return deta;
}

void EcalClustersGraph::initWindows(double etaWidth, double phiWidth){

    for (int is=0; is < nSeeds_; is++){
        double seed_eta = clusters_.at(is)->eta();
        double seed_phi = clusters_.at(is)->phi();
        inWindows_(is,is) = 1;

        for (int icl=is+1 ; icl < nCls_; icl++){
            double cl_eta = clusters_.at(icl)->eta();
            double cl_phi = clusters_.at(icl)->phi();
            double dphi = deltaPhi(seed_phi,cl_phi); 
            double deta = deltaEta(seed_eta,cl_eta);

            if( deta <= etaWidth && dphi <= phiWidth){
                inWindows_(is, icl) = 1;
                //Save also symmetric part of the adj matrix
                if (icl < nSeeds_)  inWindows_(icl,is) = 1;
            }else{
                 inWindows_(is, icl) = 0;
                 if (icl < nSeeds_)  inWindows_(icl,is) = 0;
            }
        }
    }

    //test
    /*for (int is=0; is < nSeeds_; is++){
        for (int ic=0; ic < nCls_; ic++){
            cout << inWindows_(is,ic) << " ";
        }
        cout << endl;
    }*/

}

void EcalClustersGraph::clearWindows()
{
   inWindows_.clear();
}

std::pair<double,double> EcalClustersGraph::computeCovariances(const CaloCluster* cluster) 
{

     double etaWidth = 0.;
     double phiWidth = 0.;
     double numeratorEtaWidth = 0;
     double numeratorPhiWidth = 0;

     double clEnergy = cluster->energy();
     double denominator = clEnergy;

     double clEta = cluster->position().eta();
     double clPhi = cluster->position().phi();

     std::shared_ptr<const CaloCellGeometry> this_cell;
     EcalRecHitCollection::const_iterator rHit;

     const std::vector<std::pair<DetId, float> >& detId = cluster->hitsAndFractions();
     // Loop over recHits associated with the given SuperCluster
     for (std::vector<std::pair<DetId, float> >::const_iterator hit = detId.begin(); hit != detId.end(); ++hit) {     
       if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_BARREL){
          rHit = recHitsEB_->find((*hit).first);
          //FIXME: THIS IS JUST A WORKAROUND A FIX SHOULD BE APPLIED
          if (rHit == recHitsEB_->end()) {
              continue;
          }  
       }else if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_ENDCAP){
          rHit = recHitsEE_->find((*hit).first);
          //FIXME: THIS IS JUST A WORKAROUND A FIX SHOULD BE APPLIED
          if (rHit == recHitsEE_->end()) {
              continue;
          }  
       }  

       if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_BARREL){
          this_cell = ebGeom_->getGeometry(rHit->id());
       }else if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_ENDCAP){
          this_cell = eeGeom_->getGeometry(rHit->id());
       } 
       if (this_cell == nullptr) {
         //edm::LogInfo("SuperClusterShapeAlgo") << "pointer to the cell in Calculate_Covariances is NULL!";
         continue;
       }

       GlobalPoint position = this_cell->getPosition();
       //take into account energy fractions
       double energyHit = rHit->energy() * hit->second; 

       //form differences
       double dPhi = position.phi() - clPhi;
       if (dPhi > +Geom::pi()) {
         dPhi = Geom::twoPi() - dPhi;
       }
       if (dPhi < -Geom::pi()) {
         dPhi = Geom::twoPi() + dPhi;
       }

       double dEta = position.eta() - clEta;

       if (energyHit > 0) {
         numeratorEtaWidth += energyHit * dEta * dEta;
         numeratorPhiWidth += energyHit * dPhi * dPhi;
       }

       etaWidth = sqrt(numeratorEtaWidth / denominator);
       phiWidth = sqrt(numeratorPhiWidth / denominator);
     }

     return std::make_pair(etaWidth,phiWidth);
} 

void EcalClustersGraph::computeVariables(const CaloCluster* seed, const CaloCluster* cluster)
{
     NNclusterVars_.clear();
     NNclusterVars_.resize(meanVals_->size());

     double zSide = 0;
     if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_ENDCAP && seed->eta()<0.) zSide = -1.;   
     if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_ENDCAP && seed->eta()>0.) zSide = +1.; 

     NNclusterVars_[0] = (((cluster==seed) ? 1 : 0)-meanVals_->at(0))/stdVals_->at(0); //isSeed
     NNclusterVars_[1] = (deltaEta(seed->eta(), cluster->eta())-meanVals_->at(1))/stdVals_->at(1); //cl_dEta  
     NNclusterVars_[2] = (deltaPhi(seed->phi(), cluster->phi())-meanVals_->at(2))/stdVals_->at(2); //cl_dPhi
     NNclusterVars_[3] = (cluster->energy()-meanVals_->at(3))/stdVals_->at(3); //cl_energy
     NNclusterVars_[4] = (cluster->energy()/TMath::CosH(cluster->eta())-meanVals_->at(4))/stdVals_->at(4); //cl_et
     if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_BARREL){ 
        showerShapes_ = computeShowerShapes(cluster);
        NNclusterVars_[5] = (showerShapes_[0]-meanVals_->at(5))/stdVals_->at(5); //cl_r9
        NNclusterVars_[6] = (showerShapes_[1]-meanVals_->at(6))/stdVals_->at(6); //cl_sigmaietaieta
        NNclusterVars_[7] = (showerShapes_[2]-meanVals_->at(7))/stdVals_->at(7); //cl_sigmaietaiphi
        NNclusterVars_[8] = (showerShapes_[3]-meanVals_->at(8))/stdVals_->at(8); //cl_sigmaiphiiphi     
        NNclusterVars_[9] = (showerShapes_[4]-meanVals_->at(9))/stdVals_->at(9); //cl_swiss_cross      
        NNclusterVars_[10] = (showerShapes_[5]-meanVals_->at(10))/stdVals_->at(10); //cl_nXtals 
        NNclusterVars_[11] = (showerShapes_[6]-meanVals_->at(11))/stdVals_->at(11); //cl_etaWidth 
        NNclusterVars_[12] = (showerShapes_[7]-meanVals_->at(12))/stdVals_->at(12); //cl_phiWidth 
     }else if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_ENDCAP){
        showerShapes_ = computeShowerShapes(cluster);
        NNclusterVars_[5] = (showerShapes_[0]-meanVals_->at(5))/stdVals_->at(5); //cl_r9
        NNclusterVars_[6] = (showerShapes_[1]-meanVals_->at(6))/stdVals_->at(6); //cl_sigmaietaieta
        NNclusterVars_[7] = (showerShapes_[2]-meanVals_->at(7))/stdVals_->at(7); //cl_sigmaietaiphi
        NNclusterVars_[8] = (showerShapes_[3]-meanVals_->at(8))/stdVals_->at(8); //cl_sigmaiphiiphi     
        NNclusterVars_[9] = (showerShapes_[4]-meanVals_->at(9))/stdVals_->at(9); //cl_swiss_cross      
        NNclusterVars_[10] = (showerShapes_[5]-meanVals_->at(10))/stdVals_->at(10); //cl_nXtals 
        NNclusterVars_[11] = (showerShapes_[6]-meanVals_->at(11))/stdVals_->at(11); //cl_etaWidth 
        NNclusterVars_[12] = (showerShapes_[7]-meanVals_->at(12))/stdVals_->at(12); //cl_phiWidth 
     }
     NNclusterVars_[13] = (seed->eta()-meanVals_->at(13))/stdVals_->at(13); //seed_eta
     NNclusterVars_[14] = (zSide-meanVals_->at(14))/stdVals_->at(14); //seed_iz
     NNclusterVars_[15] = (seed->energy()-meanVals_->at(15))/stdVals_->at(15); //seed_energy
     NNclusterVars_[16] = (seed->energy()/TMath::CosH(seed->eta())-meanVals_->at(16))/stdVals_->at(16); //seed_et
     if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_BARREL){ 
        showerShapes_ = computeShowerShapes(seed);
        NNclusterVars_[17] = (showerShapes_[0]-meanVals_->at(17))/stdVals_->at(17); //seed_r9
        NNclusterVars_[18] = (showerShapes_[1]-meanVals_->at(18))/stdVals_->at(18); //seed_sigmaietaieta
        NNclusterVars_[19] = (showerShapes_[2]-meanVals_->at(19))/stdVals_->at(19); //seed_sigmaietaiphi
        NNclusterVars_[20] = (showerShapes_[3]-meanVals_->at(20))/stdVals_->at(20); //seed_sigmaiphiiphi     
        NNclusterVars_[21] = (showerShapes_[4]-meanVals_->at(21))/stdVals_->at(21); //seed_swiss_cross      
        NNclusterVars_[22] = (showerShapes_[5]-meanVals_->at(22))/stdVals_->at(22); //seed_nXtals 
        NNclusterVars_[23] = (showerShapes_[6]-meanVals_->at(23))/stdVals_->at(23); //seed_etaWidth 
        NNclusterVars_[24] = (showerShapes_[7]-meanVals_->at(24))/stdVals_->at(24); //seed_phiWidth        
     }else if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_ENDCAP){    
        showerShapes_ = computeShowerShapes(seed);
        NNclusterVars_[17] = (showerShapes_[0]-meanVals_->at(17))/stdVals_->at(17); //seed_r9
        NNclusterVars_[18] = (showerShapes_[1]-meanVals_->at(18))/stdVals_->at(18); //seed_sigmaietaieta
        NNclusterVars_[19] = (showerShapes_[2]-meanVals_->at(19))/stdVals_->at(19); //seed_sigmaietaiphi
        NNclusterVars_[20] = (showerShapes_[3]-meanVals_->at(20))/stdVals_->at(20); //seed_sigmaiphiiphi     
        NNclusterVars_[21] = (showerShapes_[4]-meanVals_->at(21))/stdVals_->at(21); //seed_swiss_cross      
        NNclusterVars_[22] = (showerShapes_[5]-meanVals_->at(22))/stdVals_->at(22); //seed_nXtals 
        NNclusterVars_[23] = (showerShapes_[6]-meanVals_->at(23))/stdVals_->at(23); //seed_etaWidth 
        NNclusterVars_[24] = (showerShapes_[7]-meanVals_->at(24))/stdVals_->at(24); //seed_phiWidth    
     }
}

std::vector<double> EcalClustersGraph::computeShowerShapes(const CaloCluster* cluster)
{
     std::vector<double> showerVars_;
     showerVars_.resize(8);
     float e1=1.; 
     float e4=0.;

     full5x5_locCov_.clear();      
     if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_BARREL){ 
        full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(*cluster, recHitsEB_, topology_);
        widths_ = computeCovariances(cluster);  
        e1 = noZS::EcalClusterTools::eMax(*cluster, recHitsEB_);
        e4 = noZS::EcalClusterTools::eTop(*cluster, recHitsEB_, topology_) +
             noZS::EcalClusterTools::eRight(*cluster, recHitsEB_, topology_) +
             noZS::EcalClusterTools::eBottom(*cluster, recHitsEB_, topology_) +
             noZS::EcalClusterTools::eLeft(*cluster, recHitsEB_, topology_); 
        showerVars_[0] = noZS::EcalClusterTools::e3x3(*cluster, recHitsEB_, topology_)/cluster->energy(); //r9
        showerVars_[1] = sqrt(full5x5_locCov_[0]); //sigmaietaieta
        showerVars_[2] = full5x5_locCov_[1]; //sigmaietaiphi
        showerVars_[3] = (!edm::isFinite(full5x5_locCov_[2])) ? 0. : sqrt(full5x5_locCov_[2]); //sigmaiphiiphi     
        showerVars_[4] = (e1!=0.) ? 1.-e4/e1 : -999.; //swiss_cross      
        showerVars_[5] = cluster->hitsAndFractions().size(); //nXtals 
        showerVars_[6] = widths_.first; //etaWidth 
        showerVars_[7] = widths_.second; //phiWidth   
     }else if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_ENDCAP){
        full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(*cluster, recHitsEE_, topology_);
        widths_ = computeCovariances(cluster);  
        e1 = noZS::EcalClusterTools::eMax(*cluster, recHitsEE_);
        e4 = noZS::EcalClusterTools::eTop(*cluster, recHitsEE_, topology_) +
             noZS::EcalClusterTools::eRight(*cluster, recHitsEE_, topology_) +
             noZS::EcalClusterTools::eBottom(*cluster, recHitsEE_, topology_) +
             noZS::EcalClusterTools::eLeft(*cluster, recHitsEE_, topology_); 
        showerVars_[0] = noZS::EcalClusterTools::e3x3(*cluster, recHitsEE_, topology_)/cluster->energy(); //r9
        showerVars_[1] = sqrt(full5x5_locCov_[0]); //sigmaietaieta
        showerVars_[2] = full5x5_locCov_[1]; //sigmaietaiphi
        showerVars_[3] = (!edm::isFinite(full5x5_locCov_[2])) ? 0. : sqrt(full5x5_locCov_[2]); //sigmaiphiiphi    
        showerVars_[4] = (e1!=0.) ? 1.-e4/e1 : -999.; //swiss_cross      
        showerVars_[5] = cluster->hitsAndFractions().size(); //nXtals 
        showerVars_[6] = widths_.first; //etaWidth 
        showerVars_[7] = widths_.second; //phiWidth  
     }

     return showerVars_; 
}



void EcalClustersGraph::fillVariables()
{
     NNwindowVars_.clear();
     NNwindowVars_.resize(nSeeds_);
     for (int is=0; is < nSeeds_; is++){
        for (int ic=0; ic < nCls_; ic++){
             if(inWindows_(is,ic)==1){ 
                computeVariables((reco::CaloCluster*)(*clusters_.at(is)).the_ptr().get(),(reco::CaloCluster*)(*clusters_.at(ic)).the_ptr().get());     
                NNwindowVars_[is].push_back(NNclusterVars_);
             }     
        }
    }
}



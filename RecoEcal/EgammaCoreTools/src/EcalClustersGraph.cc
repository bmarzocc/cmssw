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

std::vector<int> EcalClustersGraph::clusterPosition(const CaloCluster* cluster)
{
     std::vector<int> coordinates;
     coordinates.resize(3);
     int ieta = -999;
     int iphi = -999;
     int iz = -99;

     math::XYZPoint caloPos = cluster->position();
     if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_BARREL){  
        EBDetId eb_id(ebGeom_->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
        ieta = eb_id.ieta();
        iphi = eb_id.iphi();
        iz = 0;
     }else if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_ENDCAP){ 
        EEDetId ee_id(eeGeom_->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));    
        ieta = ee_id.ix();
        iphi = ee_id.iy();
        if(ee_id.zside()<0) iz=-1;
        if(ee_id.zside()>0) iz=1;   
     }  

     coordinates[0] = ieta;
     coordinates[1] = iphi;
     coordinates[2] = iz;  
     return coordinates;     
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

std::vector<double> EcalClustersGraph::dynamicWindow(double seedEta)
{
     std::vector<double> window;
     window.resize(3);

     double eta = fabs(seedEta);
     double deta_down = 0.;
     double deta_up = 0.;
     double dphi = 0.;
     
     //deta_down
     if(eta < 2.1) 
        deta_down = -0.075;
     else if(eta >= 2.1 && eta < 2.5)
        deta_down = -0.1875 *eta + 0.31875;
     else if(eta >=2.5)
        deta_down = -0.15;

     //deta_up
     if(eta >= 0 && eta < 0.1)
        deta_up = 0.075;
     else if(eta >= 0.1 && eta < 1.3)
        deta_up = 0.0758929 -0.0178571* eta + 0.0892857*(eta*eta); 
     else if(eta >= 1.3 && eta < 1.7)
        deta_up = 0.2;
     else if(eta >=1.7 && eta < 1.9)
        deta_up = 0.625 -0.25*eta;
     else if(eta >= 1.9)
        deta_up = 0.15;
            
     //dphi
     if(eta < 1.9)
        dphi = 0.6;
     else if(eta >= 1.9 && eta < 2.7)
        dphi = 1.075 -0.25 * eta;
     else if(eta >= 2.7)
        dphi = 0.4;
     
     window[0] = deta_down;
     window[1] = deta_up;    
     window[2] = dphi;  
     
     return window;

}

void EcalClustersGraph::initWindows(){

    for (int is=0; is < nSeeds_; is++){
        std::vector<int> seedLocal = clusterPosition((*clusters_.at(is)).the_ptr().get());
        double seed_eta = clusters_.at(is)->eta();
        double seed_phi = clusters_.at(is)->phi();
        inWindows_(is,is) = 1;
        std::vector<double> width = dynamicWindow(seed_eta);

        for (int icl=is+1 ; icl < nCls_; icl++){
            std::vector<int> clusterLocal = clusterPosition((*clusters_.at(icl)).the_ptr().get());
            double cl_eta = clusters_.at(icl)->eta();
            double cl_phi = clusters_.at(icl)->phi();
            double dphi = deltaPhi(seed_phi,cl_phi); 
            double deta = deltaEta(seed_eta,cl_eta);
           
            int isIn=0;
            if(seedLocal[2]==clusterLocal[2] && deta>=width[0] && deta<=width[1] && fabs(dphi)<=width[2]) isIn = 1; 
            inWindows_(is, icl) = isIn;
            //Save also symmetric part of the adj matrix
            if (icl < nSeeds_) inWindows_(icl,is) = isIn;
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

void EcalClustersGraph::computeVariables(const CaloCluster* seed, const CaloCluster* cluster)
{
     NNclusterVars_.clear();
     NNclusterVars_.resize(meanVals_->size());
     
     showerShapes_ = computeShowerShapes(cluster,false);
     std::vector<int> clusterLocal = clusterPosition(cluster);

     NNclusterVars_[0] = (((cluster==seed) ? 1 : 0)-meanVals_->at(0))/stdVals_->at(0); //isSeed
     NNclusterVars_[1] = (cluster->energy()-meanVals_->at(1))/stdVals_->at(1); //cl_energy
     NNclusterVars_[2] = (cluster->energy()/TMath::CosH(cluster->eta())-meanVals_->at(2))/stdVals_->at(2); //cl_et
     NNclusterVars_[3] = (cluster->eta()-meanVals_->at(3))/stdVals_->at(3); //cl_eta  
     NNclusterVars_[4] = (cluster->phi()-meanVals_->at(4))/stdVals_->at(4); //cl_phi  
     NNclusterVars_[5] = (clusterLocal[0]-meanVals_->at(5))/stdVals_->at(5); //cl_ieta/ix  
     NNclusterVars_[6] = (clusterLocal[1]-meanVals_->at(6))/stdVals_->at(6); //cl_iphi/iy  
     NNclusterVars_[7] = (clusterLocal[2]-meanVals_->at(7))/stdVals_->at(7); //cl_iz  
     NNclusterVars_[8] = (deltaEta(seed->eta(), cluster->eta())-meanVals_->at(8))/stdVals_->at(8); //cl_dEta  
     NNclusterVars_[9] = (deltaPhi(seed->phi(), cluster->phi())-meanVals_->at(9))/stdVals_->at(9); //cl_dPhi
     NNclusterVars_[10] = ((seed->energy()-cluster->energy())-meanVals_->at(10))/stdVals_->at(10); //cl_dEnergy
     NNclusterVars_[11] = ((seed->energy()/TMath::CosH(seed->eta())-cluster->energy()/TMath::CosH(cluster->eta()))-meanVals_->at(11))/stdVals_->at(11); //cl_dEt
     NNclusterVars_[12] = (showerShapes_[0]-meanVals_->at(12))/stdVals_->at(12); //cl_r9
     NNclusterVars_[13] = (showerShapes_[1]-meanVals_->at(13))/stdVals_->at(13); //cl_sigmaietaieta
     NNclusterVars_[14] = (showerShapes_[2]-meanVals_->at(14))/stdVals_->at(14); //cl_sigmaietaiphi
     NNclusterVars_[15] = (showerShapes_[3]-meanVals_->at(15))/stdVals_->at(15); //cl_sigmaiphiiphi     
     NNclusterVars_[16] = (showerShapes_[4]-meanVals_->at(16))/stdVals_->at(16); //cl_swiss_cross      
     NNclusterVars_[17] = (showerShapes_[5]-meanVals_->at(17))/stdVals_->at(17); //cl_nXtals 
     NNclusterVars_[18] = (showerShapes_[6]-meanVals_->at(18))/stdVals_->at(18); //cl_etaWidth 
     NNclusterVars_[19] = (showerShapes_[7]-meanVals_->at(19))/stdVals_->at(19); //cl_phiWidth 
     
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

std::vector<double> EcalClustersGraph::computeShowerShapes(const CaloCluster* cluster, bool full5x5=false)
{
     std::vector<double> showerVars_;
     showerVars_.resize(8);
     float e1=1.; 
     float e4=0.;

     locCov_.clear();     
     if(full5x5){ 
        if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_BARREL){ 
           locCov_ = noZS::EcalClusterTools::localCovariances(*cluster, recHitsEB_, topology_);
           widths_ = computeCovariances(cluster);  
           e1 = noZS::EcalClusterTools::eMax(*cluster, recHitsEB_);
           e4 = noZS::EcalClusterTools::eTop(*cluster, recHitsEB_, topology_) +
                noZS::EcalClusterTools::eRight(*cluster, recHitsEB_, topology_) +
                noZS::EcalClusterTools::eBottom(*cluster, recHitsEB_, topology_) +
                noZS::EcalClusterTools::eLeft(*cluster, recHitsEB_, topology_); 
           showerVars_[0] = noZS::EcalClusterTools::e3x3(*cluster, recHitsEB_, topology_)/cluster->energy(); //r9
           showerVars_[1] = sqrt(locCov_[0]); //sigmaietaieta
           showerVars_[2] = locCov_[1]; //sigmaietaiphi
           showerVars_[3] = (!edm::isFinite(locCov_[2])) ? 0. : sqrt(locCov_[2]); //sigmaiphiiphi     
           showerVars_[4] = (e1!=0.) ? 1.-e4/e1 : -999.; //swiss_cross      
           showerVars_[5] = cluster->hitsAndFractions().size(); //nXtals 
           showerVars_[6] = widths_.first; //etaWidth 
           showerVars_[7] = widths_.second; //phiWidth   
        }else if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_ENDCAP){
           locCov_ = noZS::EcalClusterTools::localCovariances(*cluster, recHitsEE_, topology_);
           widths_ = computeCovariances(cluster);  
           e1 = noZS::EcalClusterTools::eMax(*cluster, recHitsEE_);
           e4 = noZS::EcalClusterTools::eTop(*cluster, recHitsEE_, topology_) +
                noZS::EcalClusterTools::eRight(*cluster, recHitsEE_, topology_) +
                noZS::EcalClusterTools::eBottom(*cluster, recHitsEE_, topology_) +
                noZS::EcalClusterTools::eLeft(*cluster, recHitsEE_, topology_); 
           showerVars_[0] = noZS::EcalClusterTools::e3x3(*cluster, recHitsEE_, topology_)/cluster->energy(); //r9
           showerVars_[1] = sqrt(locCov_[0]); //sigmaietaieta
           showerVars_[2] = locCov_[1]; //sigmaietaiphi
           showerVars_[3] = (!edm::isFinite(locCov_[2])) ? 0. : sqrt(locCov_[2]); //sigmaiphiiphi    
           showerVars_[4] = (e1!=0.) ? 1.-e4/e1 : -999.; //swiss_cross      
           showerVars_[5] = cluster->hitsAndFractions().size(); //nXtals 
           showerVars_[6] = widths_.first; //etaWidth 
           showerVars_[7] = widths_.second; //phiWidth  
        }
     }else{
        if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_BARREL){ 
           locCov_ = EcalClusterTools::localCovariances(*cluster, recHitsEB_, topology_);
           widths_ = computeCovariances(cluster);  
           e1 = EcalClusterTools::eMax(*cluster, recHitsEB_);
           e4 = EcalClusterTools::eTop(*cluster, recHitsEB_, topology_) +
                EcalClusterTools::eRight(*cluster, recHitsEB_, topology_) +
                EcalClusterTools::eBottom(*cluster, recHitsEB_, topology_) +
                EcalClusterTools::eLeft(*cluster, recHitsEB_, topology_); 
           showerVars_[0] = EcalClusterTools::e3x3(*cluster, recHitsEB_, topology_)/cluster->energy(); //r9
           showerVars_[1] = sqrt(locCov_[0]); //sigmaietaieta
           showerVars_[2] = locCov_[1]; //sigmaietaiphi
           showerVars_[3] = (!edm::isFinite(locCov_[2])) ? 0. : sqrt(locCov_[2]); //sigmaiphiiphi     
           showerVars_[4] = (e1!=0.) ? 1.-e4/e1 : -999.; //swiss_cross      
           showerVars_[5] = cluster->hitsAndFractions().size(); //nXtals 
           showerVars_[6] = widths_.first; //etaWidth 
           showerVars_[7] = widths_.second; //phiWidth   
        }else if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_ENDCAP){
           locCov_ = EcalClusterTools::localCovariances(*cluster, recHitsEE_, topology_);
           widths_ = computeCovariances(cluster);  
           e1 = EcalClusterTools::eMax(*cluster, recHitsEE_);
           e4 = EcalClusterTools::eTop(*cluster, recHitsEE_, topology_) +
                EcalClusterTools::eRight(*cluster, recHitsEE_, topology_) +
                EcalClusterTools::eBottom(*cluster, recHitsEE_, topology_) +
                EcalClusterTools::eLeft(*cluster, recHitsEE_, topology_); 
           showerVars_[0] = EcalClusterTools::e3x3(*cluster, recHitsEE_, topology_)/cluster->energy(); //r9
           showerVars_[1] = sqrt(locCov_[0]); //sigmaietaieta
           showerVars_[2] = locCov_[1]; //sigmaietaiphi
           showerVars_[3] = (!edm::isFinite(locCov_[2])) ? 0. : sqrt(locCov_[2]); //sigmaiphiiphi    
           showerVars_[4] = (e1!=0.) ? 1.-e4/e1 : -999.; //swiss_cross      
           showerVars_[5] = cluster->hitsAndFractions().size(); //nXtals 
           showerVars_[6] = widths_.first; //etaWidth 
           showerVars_[7] = widths_.second; //phiWidth  
        }
     }

     return showerVars_; 
}

void EcalClustersGraph::fillHits(const CaloCluster* cluster)
{
     NNclusterHits_.clear();
     std::vector<double> info_;
     
     const std::vector<std::pair<DetId,float> > &hitsAndFractions = cluster->hitsAndFractions(); 
     for(unsigned int i=0; i<hitsAndFractions.size(); i++){
         info_.clear();
         info_.resize(6);         
         if(hitsAndFractions[i].first.subdetId()==EcalBarrel){ 
            double energy = (*recHitsEB_->find(hitsAndFractions[i].first)).energy();     
            EBDetId eb_id(hitsAndFractions[i].first); 
            info_[0] = eb_id.ieta(); //ieta
            info_[1] = eb_id.iphi(); //iphi
            info_[2] = 0.; //iz
            info_[3] = energy; //energy
            info_[4] = energy*hitsAndFractions[i].second; //energy * fraction
            info_[5] = hitsAndFractions[i].second; //fraction 
         }else if(hitsAndFractions[i].first.subdetId()==EcalEndcap){  
            double energy = (*recHitsEE_->find(hitsAndFractions[i].first)).energy();     
            EEDetId ee_id(hitsAndFractions[i].first); 
            info_[0] = ee_id.ix(); //ix
            info_[1] = ee_id.iy(); //iy
            info_[3] = energy; //energy
            info_[4] = energy*hitsAndFractions[i].second;//energy * fraction
            info_[5] = hitsAndFractions[i].second; //fraction  
            if(ee_id.zside()<0) info_[2] = -1.; //iz
            if(ee_id.zside()>0) info_[2] = +1.; //iz 
         }
         NNclusterHits_.push_back(info_);          
     } 
}

void EcalClustersGraph::fillVariables()
{
     NNwindowVars_.clear();
     NNwindowVars_.resize(nSeeds_);
     for (int is=0; is < nSeeds_; is++){
        for (int ic=0; ic < nCls_; ic++){
             if(inWindows_(is,ic)==1){ 
                computeVariables((*clusters_.at(is)).the_ptr().get(),(*clusters_.at(ic)).the_ptr().get()); 
                fillHits((*clusters_.at(ic)).the_ptr().get());     
                NNwindowVars_[is].push_back(std::make_pair(NNclusterVars_,NNclusterHits_));
             }     
        }
    }
}



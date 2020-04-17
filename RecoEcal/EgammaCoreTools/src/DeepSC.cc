#include "RecoEcal/EgammaCoreTools/interface/DeepSC.h"
#include "TMath.h"
#include "TVector2.h"
#include <cmath>
using namespace std;

namespace reco {  
  
  void DeepSC::SetXtalsInWindow(DetId seedDetId, double etawidthSuperCluster, double phiwidthSuperCluster, const CaloTopology *topology, const CaloGeometry *geometry)
  {   
    std::vector<DetId> v_id;
    double maxDeta = 0.;
    double maxDphi = 0.;
    int window_index = 0; 
    while(maxDeta<etawidthSuperCluster || maxDphi<phiwidthSuperCluster){
      window_index++;
      v_id = EcalClusterTools::matrixDetId(topology, seedDetId, window_index);
      std::pair<double,double> DetaDphi_pair = GetMaximumDetaDphi(&seedDetId,&v_id,geometry);  
      maxDeta = DetaDphi_pair.first;
      maxDphi = DetaDphi_pair.second; 
    }

    for(unsigned int iXtal=0; iXtal<v_id.size(); iXtal++){
      double eta_seed = (geometry->getPosition(seedDetId)).eta();  
      double phi_seed = (geometry->getPosition(seedDetId)).phi();
      double eta_xtal = (geometry->getPosition(v_id.at(iXtal))).eta();
      double phi_xtal = (geometry->getPosition(v_id.at(iXtal))).phi(); 
      if(fabs(eta_xtal-eta_seed)<etawidthSuperCluster && fabs(TVector2::Phi_mpi_pi(phi_xtal-phi_seed))<phiwidthSuperCluster && seedDetId.subdetId()==v_id.at(iXtal).subdetId()) xtals_inWindow_.push_back(v_id.at(iXtal));
    }
  }

  std::pair<double,double> DeepSC::GetMaximumDetaDphi(DetId* seedDetId, std::vector<DetId>* idMatrix, const CaloGeometry *geometry)
  {
    double maxDeta=0.;
    double maxDphi=0.;
    for(unsigned int iXtal=0; iXtal<idMatrix->size(); iXtal++){
        double eta_seed = (geometry->getPosition(*seedDetId)).eta();  
        double phi_seed = (geometry->getPosition(*seedDetId)).phi();
        double eta_xtal = (geometry->getPosition(idMatrix->at(iXtal))).eta();
        double phi_xtal = (geometry->getPosition(idMatrix->at(iXtal))).phi(); 
        if(fabs(eta_xtal-eta_seed)>maxDeta) maxDeta=fabs(eta_xtal-eta_seed);
        if(fabs(TVector2::Phi_mpi_pi(phi_xtal-phi_seed))>maxDphi) maxDphi=fabs(TVector2::Phi_mpi_pi(phi_xtal-phi_seed));
    }
    return std::make_pair(maxDeta,maxDphi);
  }
 
  void DeepSC::ClearXtalsInWindow(){ xtals_inWindow_.clear(); }  
  
  void DeepSC::SetNNGraphSession(std::string NNfileName_)
  {
    graphDef_NN_ = tensorflow::loadGraphDef(NNfileName_.c_str());
    session_NN_ = tensorflow::createSession(graphDef_NN_); 
  }

  void DeepSC::SetNNEvaluationParams(std::string NNinput, std::string NNoutput)
  {
    NNinput_string_ = NNinput; 
    NNoutput_string_ = NNoutput;  
    //std::cout << "SetNNEvaluationParams: " << NNinput_string_ << " - " << NNoutput_string_ << std::endl;
  }

  void DeepSC::ClearNNnormalizeParams()
  {
    x_mean_.clear();
    x_std_.clear();
  }

  void DeepSC::ClearNNVars()
  {
    vectorVar_.clear();
  }

  void DeepSC::SetNormalizeNNParams(std::vector<double> var_mean_, std::vector<double> var_std_)
  {
    if(var_mean_.size()==var_std_.size()){
       x_mean_ = var_mean_;
       x_std_ = var_std_;

      //for(unsigned int iVar=0; iVar<x_mean_.size(); iVar++)
      //    std::cout << "DeepSC::SetNormalizeNNParams: " << iVar << " - var_mean=" << x_mean_.at(iVar) << " - var_std=" << x_std_.at(iVar) << std::endl;
    }else{
       std::cout << "DeepSC::SetNormalizeNNParams WARNING: sizes don't match!" << std::endl;
    }
  }

  void DeepSC::SetNNVarVal(std::vector<double> vars)
  {
    vectorVar_.resize(x_mean_.size());
    if(vars.size()==x_mean_.size()){
       for(unsigned int iVar=0; iVar<vectorVar_.size(); iVar++){   
           vectorVar_[iVar] = vars.at(iVar); 
           //std::cout << "DeepSC::SetNNVarVal: "  << iVar << " - vars=" << vectorVar_[iVar] << std::endl; 
       }
    }else{
       std::cout << "DeepSC::SetNNVarVal WARNING: sizes don't match!" << std::endl;
    }  
  }  

  void DeepSC::NormalizeNNVars()
  {
    for(unsigned int iVar=0; iVar<vectorVar_.size(); iVar++){
        vectorVar_[iVar] = (vectorVar_[iVar] - x_mean_[iVar])/x_std_[iVar];  
        //std::cout << "DeepSC::NormalizeNNVars: "  << iVar << " - vars=" << vectorVar_[iVar] << std::endl;
    }    
  }  

 void DeepSC::SetClusterVariables(const edm::Ptr<reco::PFCluster> cluster, const CaloTopology *topology, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE)
  {
    NNclusterVars_.clear();
    NNclusterVars_.resize(7);
    full5x5_locCov_.clear();
    const reco::BasicCluster caloBC_(*cluster);  
    float e1=1.; 
    float e4=0.;

    if(cluster->layer() == PFLayer::ECAL_BARREL){ 
       full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(caloBC_, recHitsEB, topology);
       e1 = EcalClusterTools::eMax(caloBC_, recHitsEB);
       e4 = EcalClusterTools::eTop(caloBC_, recHitsEB, topology) +
            EcalClusterTools::eRight(caloBC_, recHitsEB, topology) +
            EcalClusterTools::eBottom(caloBC_, recHitsEB, topology) +
            EcalClusterTools::eLeft(caloBC_, recHitsEB, topology);  
       NNclusterVars_[0] = cluster->energy(); // energy    
       NNclusterVars_[1] = noZS::EcalClusterTools::e3x3(caloBC_, recHitsEB, topology)/cluster->energy(); //r9
       NNclusterVars_[2] = sqrt(full5x5_locCov_[0]); //sigmaietaieta
       NNclusterVars_[3] = full5x5_locCov_[1]; //sigmaietaiphi
       NNclusterVars_[4] = (!edm::isFinite(full5x5_locCov_[2])) ? 0. : sqrt(full5x5_locCov_[2]); //sigmaiphiiphi     
       NNclusterVars_[5] = 1.-e4/e1; //swiss_cross      
       NNclusterVars_[6] = cluster->hitsAndFractions().size(); //nXtals 
    }else if(cluster->layer() == PFLayer::ECAL_ENDCAP){
       full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(caloBC_, recHitsEE, topology);
       e1 = EcalClusterTools::eMax(caloBC_, recHitsEE);
       e4 = EcalClusterTools::eTop(caloBC_, recHitsEE, topology) +
            EcalClusterTools::eRight(caloBC_, recHitsEE, topology) +
            EcalClusterTools::eBottom(caloBC_, recHitsEE, topology) +
            EcalClusterTools::eLeft(caloBC_, recHitsEE, topology);
       NNclusterVars_[0] = cluster->energy(); // energy    
       NNclusterVars_[1] = noZS::EcalClusterTools::e3x3(caloBC_, recHitsEE, topology)/cluster->energy(); //r9
       NNclusterVars_[2] = sqrt(full5x5_locCov_[0]); //sigmaietaieta
       NNclusterVars_[3] = full5x5_locCov_[1]; //sigmaietaiphi
       NNclusterVars_[4] = (!edm::isFinite(full5x5_locCov_[2])) ? 0. : sqrt(full5x5_locCov_[2]); //sigmaiphiiphi     
       NNclusterVars_[5] = 1.-e4/e1; //swiss_cross      
       NNclusterVars_[6] = cluster->hitsAndFractions().size(); //nXtals
    }

    //std::cout << "ClusterVariables: " << NNclusterVars_[0] << " - " << NNclusterVars_[1] << " - " << NNclusterVars_[2] << " - " << NNclusterVars_[3] << " - " << NNclusterVars_[4] << " - " << NNclusterVars_[5] << " - " << NNclusterVars_[6] << std::endl; 
  }

  float DeepSC::EvaluateNN()
  {
    unsigned int shape = vectorVar_.size();
    tensorflow::Tensor NNinput(tensorflow::DT_FLOAT, {1,shape});
    for(unsigned int i = 0; i < shape; i++){
        NNinput.matrix<float>()(0,i) =  float(vectorVar_[i]);
    }

    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run(session_NN_, { {NNinput_string_.c_str(), NNinput} } , { NNoutput_string_.c_str() }, &outputs);
    float NNscore = outputs[0].matrix<float>()(0, 0);
    return NNscore;
  }
}

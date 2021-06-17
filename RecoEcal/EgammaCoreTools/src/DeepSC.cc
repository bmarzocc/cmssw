#include "RecoEcal/EgammaCoreTools/interface/DeepSC.h"
#include "TMath.h"
#include "TVector2.h"
#include <cmath>
using namespace std;

namespace reco {  
  
  DeepSC::DeepSC()
  {
     session_NN_Barrel_ = tensorflow::createSession(tensorflow::loadGraphDef(edm::FileInPath("RecoEcal/EgammaCoreTools/data/deepSC_model_v17_joindet_elegamma_EBEE.pb").fullPath().c_str()));
     session_NN_Endcap_ = tensorflow::createSession(tensorflow::loadGraphDef(edm::FileInPath("RecoEcal/EgammaCoreTools/data/deepSC_model_v17_joindet_elegamma_EBEE.pb").fullPath().c_str()));
     x_mean_Barrel_ = std::vector<double>({6.84241156e-03, 1.62242679e-03, 5.81495577e+01, 2.57215845e+01,
-7.09402501e-04, -1.27142875e-04, 1.30375508e+00, 5.67249500e-01,
1.00772582e+00, 1.35803461e-02, -4.29317013e-06, 1.71072024e-02,
4.90466869e-01, 5.10511982e+00, 8.82101138e-03, 1.04095965e-02,
1.92096066e+00, 1.31476120e-02, 1.62948213e-05, 1.42948806e-02,
5.92920497e-01, 1.49597644e+00, 3.36213188e-03, 3.06446267e-03});
     x_std_Barrel_ = std::vector<double>({1.31333380e+00, 5.06988411e-01, 9.21157365e+01, 2.98580765e+01,
1.10279784e-01, 3.30488055e-01, 2.62605247e+00, 1.16284769e+00,
1.17047757e-01, 1.11969442e-02, 1.86572967e-04, 1.31036359e-02,
4.01511744e-01, 5.67007350e+00, 6.14304203e-03, 7.24808860e-03,
7.81094814e+00, 1.70392176e-02, 3.05995567e-04, 1.80176053e-02,
1.99316624e+00, 1.88845046e+00, 4.12315715e-03, 4.79639033e-03});
     x_mean_Endcap_ = std::vector<double>({6.84241156e-03, 1.62242679e-03, 5.81495577e+01, 2.57215845e+01,
-7.09402501e-04, -1.27142875e-04, 1.30375508e+00, 5.67249500e-01,
1.00772582e+00, 1.35803461e-02, -4.29317013e-06, 1.71072024e-02,
4.90466869e-01, 5.10511982e+00, 8.82101138e-03, 1.04095965e-02,
1.92096066e+00, 1.31476120e-02, 1.62948213e-05, 1.42948806e-02,
5.92920497e-01, 1.49597644e+00, 3.36213188e-03, 3.06446267e-03});
     x_std_Endcap_ = std::vector<double>({1.31333380e+00, 5.06988411e-01, 9.21157365e+01, 2.98580765e+01,
1.10279784e-01, 3.30488055e-01, 2.62605247e+00, 1.16284769e+00,
1.17047757e-01, 1.11969442e-02, 1.86572967e-04, 1.31036359e-02,
4.01511744e-01, 5.67007350e+00, 6.14304203e-03, 7.24808860e-03,
7.81094814e+00, 1.70392176e-02, 3.05995567e-04, 1.80176053e-02,
1.99316624e+00, 1.88845046e+00, 4.12315715e-03, 4.79639033e-03}); 
     NNinput_string_Barrel_ = std::string("dense_4_input:0");
     NNoutput_string_Barrel_ = std::string("dense_6/Sigmoid:0");
     NNinput_string_Endcap_ = std::string("dense_4_input:0");
     NNoutput_string_Endcap_ = std::string("dense_6/Sigmoid:0");
     NNWPs_Barrel_ = std::vector<double>({0.75,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85});
     NNWPs_Endcap_ = std::vector<double>({0.75,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85});  
     etaWindows_Barrel_ = std::vector<double>({0.2,0.,0.}); // p0 + p1*seedEta + p2*seedEta*seedEta
     phiWindows_Barrel_ = std::vector<double>({0.6,0.,0.}); // p0 + p1*seedEta + p2*seedEta*seedEta   
     etaWindows_Endcap_ = std::vector<double>({0.2,0.,0.}); // p0 + p1*seedEta + p2*seedEta*seedEta
     phiWindows_Endcap_ = std::vector<double>({0.6,0.,0.}); // p0 + p1*seedEta + p2*seedEta*seedEta
  }
  
  DeepSC::~DeepSC()
  {
   
  }

  double DeepSC::DeltaPhi(double seed_phi, double cluster_phi)
  {
     double dphi = seed_phi - cluster_phi;
     if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
     if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
     return dphi;
  } 

  double DeepSC::DeltaEta(double seed_eta, double cluster_eta)
  {
     double deta = 0.;
     if(seed_eta > 0.) deta = cluster_eta - seed_eta;
     if(seed_eta <= 0.) deta = seed_eta - cluster_eta;
     return deta;
  }

  std::pair<double,double> DeepSC::ComputeCovariances(const CaloCluster cluster, const EcalRecHitCollection* recHits, const CaloSubdetectorGeometry* geometry) {

     double etaWidth = 0.;
     double phiWidth = 0.;
     double numeratorEtaWidth = 0;
     double numeratorPhiWidth = 0;

     double clEnergy = cluster.energy();
     double denominator = clEnergy;

     double clEta = cluster.position().eta();
     double clPhi = cluster.position().phi();

     const std::vector<std::pair<DetId, float> >& detId = cluster.hitsAndFractions();
     // Loop over recHits associated with the given SuperCluster
     for (std::vector<std::pair<DetId, float> >::const_iterator hit = detId.begin(); hit != detId.end(); ++hit) {
       EcalRecHitCollection::const_iterator rHit = recHits->find((*hit).first);
       //FIXME: THIS IS JUST A WORKAROUND A FIX SHOULD BE APPLIED
       if (rHit == recHits->end()) {
         continue;
       }
       auto this_cell = geometry->getGeometry(rHit->id());
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

  void DeepSC::ComputeVariables(const CaloCluster* seed, const CaloCluster* cluster, const CaloTopology *topology, const CaloSubdetectorGeometry* ebGeom, const CaloSubdetectorGeometry* eeGeom, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE)
  {
    NNclusterVars_.clear();
    NNclusterVars_.resize(24);
    float e1=1.; 
    float e4=0.;

    double zSide = 0;
    if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_ENDCAP && seed->eta()<0.) zSide = -1.;   
    if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_ENDCAP && seed->eta()>0.) zSide = +1.; 

    NNclusterVars_[0] = seed->eta();
    //NNclusterVars_[1] = seed->phi();
    NNclusterVars_[1] = zSide;
    NNclusterVars_[2] = seed->energy();
    NNclusterVars_[3] = seed->energy()/TMath::CosH(seed->eta());
    NNclusterVars_[4] = DeltaEta(seed->eta(), cluster->eta());  
    NNclusterVars_[5] = DeltaPhi(seed->phi(), cluster->phi()); 
    NNclusterVars_[6] = cluster->energy();
    NNclusterVars_[7] = cluster->energy()/TMath::CosH(cluster->eta());
    
    full5x5_locCov_.clear();
    const reco::BasicCluster seedBC_(*seed);  
    if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_BARREL){ 
       full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(seedBC_, recHitsEB, topology);
       widths_ = ComputeCovariances(seedBC_, recHitsEB, ebGeom); 
       e1 = noZS::EcalClusterTools::eMax(seedBC_, recHitsEB);
       e4 = noZS::EcalClusterTools::eTop(seedBC_, recHitsEB, topology) +
            noZS::EcalClusterTools::eRight(seedBC_, recHitsEB, topology) +
            noZS::EcalClusterTools::eBottom(seedBC_, recHitsEB, topology) +
            noZS::EcalClusterTools::eLeft(seedBC_, recHitsEB, topology);  
       NNclusterVars_[8] = noZS::EcalClusterTools::e3x3(seedBC_, recHitsEB, topology)/seed->energy(); //r9
       NNclusterVars_[9] = sqrt(full5x5_locCov_[0]); //sigmaietaieta
       NNclusterVars_[10] = full5x5_locCov_[1]; //sigmaietaiphi
       NNclusterVars_[11] = (!edm::isFinite(full5x5_locCov_[2])) ? 0. : sqrt(full5x5_locCov_[2]); //sigmaiphiiphi     
       NNclusterVars_[12] = 1.-e4/e1; //swiss_cross      
       NNclusterVars_[13] = seed->hitsAndFractions().size(); //nXtals 
       NNclusterVars_[14] = widths_.first; //etaWidth 
       NNclusterVars_[15] = widths_.second; //phiWidth 
    }else if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_ENDCAP){
       full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(seedBC_, recHitsEE, topology);
       widths_ = ComputeCovariances(seedBC_, recHitsEE, eeGeom); 
       e1 = noZS::EcalClusterTools::eMax(seedBC_, recHitsEE);
       e4 = noZS::EcalClusterTools::eTop(seedBC_, recHitsEE, topology) +
            noZS::EcalClusterTools::eRight(seedBC_, recHitsEE, topology) +
            noZS::EcalClusterTools::eBottom(seedBC_, recHitsEE, topology) +
            noZS::EcalClusterTools::eLeft(seedBC_, recHitsEE, topology);     
       NNclusterVars_[8] = noZS::EcalClusterTools::e3x3(seedBC_, recHitsEE, topology)/seed->energy(); //r9
       NNclusterVars_[9] = sqrt(full5x5_locCov_[0]); //sigmaietaieta
       NNclusterVars_[10] = full5x5_locCov_[1]; //sigmaietaiphi
       NNclusterVars_[11] = (!edm::isFinite(full5x5_locCov_[2])) ? 0. : sqrt(full5x5_locCov_[2]); //sigmaiphiiphi     
       NNclusterVars_[12] = 1.-e4/e1; //swiss_cross      
       NNclusterVars_[13] = seed->hitsAndFractions().size(); //nXtals
       NNclusterVars_[14] = widths_.first; //etaWidth 
       NNclusterVars_[15] = widths_.second; //phiWidth    
    }
     
    full5x5_locCov_.clear(); 
    const reco::BasicCluster caloBC_(*cluster); 
    if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_BARREL){ 
       full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(caloBC_, recHitsEB, topology);
       widths_ = ComputeCovariances(caloBC_, recHitsEB, ebGeom);  
       e1 = noZS::EcalClusterTools::eMax(caloBC_, recHitsEB);
       e4 = noZS::EcalClusterTools::eTop(caloBC_, recHitsEB, topology) +
            noZS::EcalClusterTools::eRight(caloBC_, recHitsEB, topology) +
            noZS::EcalClusterTools::eBottom(caloBC_, recHitsEB, topology) +
            noZS::EcalClusterTools::eLeft(caloBC_, recHitsEB, topology);   
       NNclusterVars_[16] = noZS::EcalClusterTools::e3x3(caloBC_, recHitsEB, topology)/cluster->energy(); //r9
       NNclusterVars_[17] = sqrt(full5x5_locCov_[0]); //sigmaietaieta
       NNclusterVars_[18] = full5x5_locCov_[1]; //sigmaietaiphi
       NNclusterVars_[19] = (!edm::isFinite(full5x5_locCov_[2])) ? 0. : sqrt(full5x5_locCov_[2]); //sigmaiphiiphi     
       NNclusterVars_[20] = 1.-e4/e1; //swiss_cross      
       NNclusterVars_[21] = cluster->hitsAndFractions().size(); //nXtals 
       NNclusterVars_[22] = widths_.first; //etaWidth 
       NNclusterVars_[23] = widths_.second; //phiWidth    
    }else if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_ENDCAP){
       full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(caloBC_, recHitsEE, topology);
       widths_ = ComputeCovariances(caloBC_, recHitsEE, eeGeom);  
       e1 = noZS::EcalClusterTools::eMax(caloBC_, recHitsEE);
       e4 = noZS::EcalClusterTools::eTop(caloBC_, recHitsEE, topology) +
            noZS::EcalClusterTools::eRight(caloBC_, recHitsEE, topology) +
            noZS::EcalClusterTools::eBottom(caloBC_, recHitsEE, topology) +
            noZS::EcalClusterTools::eLeft(caloBC_, recHitsEE, topology);
       NNclusterVars_[16] = noZS::EcalClusterTools::e3x3(caloBC_, recHitsEE, topology)/cluster->energy(); //r9
       NNclusterVars_[17] = sqrt(full5x5_locCov_[0]); //sigmaietaieta
       NNclusterVars_[18] = full5x5_locCov_[1]; //sigmaietaiphi
       NNclusterVars_[19] = (!edm::isFinite(full5x5_locCov_[2])) ? 0. : sqrt(full5x5_locCov_[2]); //sigmaiphiiphi     
       NNclusterVars_[20] = 1.-e4/e1; //swiss_cross      
       NNclusterVars_[21] = cluster->hitsAndFractions().size(); //nXtals
       NNclusterVars_[22] = widths_.first; //etaWidth 
       NNclusterVars_[23] = widths_.second; //phiWidth    
    } 
  }

  void DeepSC::SetNNVarVal(std::vector<double> vars)
  {
    vectorVar_.clear();
    vectorVar_.resize(vars.size());
    for(unsigned int iVar=0; iVar<vectorVar_.size(); iVar++){   
        vectorVar_[iVar] = vars.at(iVar); 
    } 
  }  

  void DeepSC::NormalizeNNVars(const CaloCluster* seed)
  {
    if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_BARREL){ 
       for(unsigned int iVar=0; iVar<vectorVar_.size(); iVar++)
           vectorVar_[iVar] = (vectorVar_[iVar] - x_mean_Barrel_[iVar])/x_std_Barrel_[iVar];  
    }else if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_ENDCAP){
       for(unsigned int iVar=0; iVar<vectorVar_.size(); iVar++)
           vectorVar_[iVar] = (vectorVar_[iVar] - x_mean_Endcap_[iVar])/x_std_Endcap_[iVar];  
    }     
  }  

  float DeepSC::EvaluateNN(const CaloCluster* seed)
  {
    unsigned int shape = vectorVar_.size();
    tensorflow::Tensor NNinput(tensorflow::DT_FLOAT, {1,shape});
    for(unsigned int i = 0; i < shape; i++){
        NNinput.matrix<float>()(0,i) =  float(vectorVar_[i]);
    }

    std::vector<tensorflow::Tensor> outputs;
    if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_BARREL){
       tensorflow::run(session_NN_Barrel_, { {NNinput_string_Barrel_.c_str(), NNinput} } , { NNoutput_string_Barrel_.c_str() }, &outputs);
    }else if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_ENDCAP){
       tensorflow::run(session_NN_Endcap_, { {NNinput_string_Endcap_.c_str(), NNinput} } , { NNoutput_string_Endcap_.c_str() }, &outputs);
    }
    float NNscore = outputs[0].matrix<float>()(0, 0);
    return NNscore;
  }

  bool DeepSC::InSuperCluster(const CaloCluster* seed, const CaloCluster* cluster, const CaloTopology *topology, const CaloSubdetectorGeometry* ebGeom, const CaloSubdetectorGeometry* eeGeom, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE)
  {
     ComputeVariables(seed, cluster, topology, ebGeom, eeGeom, recHitsEB, recHitsEE); 
     SetNNVarVal(NNclusterVars_);
     NormalizeNNVars(seed);

     bool isIn = false; 
     bool pass = false;
     double etaWindow = 0.;
     double phiWindow = 0.;
     double seedEt = seed->energy()/TMath::CosH(seed->eta());
     double seedEta = seed->eta(); 
     double deta = DeltaEta(seed->eta(),cluster->eta());  
     double dphi = DeltaPhi(seed->phi(),cluster->phi());
     float score = EvaluateNN(seed);
     
     if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_BARREL){ 

        etaWindow = etaWindows_Barrel_[0] + etaWindows_Barrel_[1]*fabs(seedEta) + etaWindows_Barrel_[2]*fabs(seedEta)*fabs(seedEta);  
        phiWindow = phiWindows_Barrel_[0] + phiWindows_Barrel_[1]*fabs(seedEta) + phiWindows_Barrel_[2]*fabs(seedEta)*fabs(seedEta);
        if(deta < etaWindow && dphi < phiWindow) isIn = true;
        
        if(seedEt>0. && seedEt<=10. && score > NNWPs_Barrel_[0]) pass = true;  
        else if(seedEt>10. && seedEt<=20. && score > NNWPs_Barrel_[1]) pass = true;  
        else if(seedEt>20. && seedEt<=30. && score > NNWPs_Barrel_[2]) pass = true;  
        else if(seedEt>30. && seedEt<=40. && score > NNWPs_Barrel_[3]) pass = true;  
        else if(seedEt>40. && seedEt<=50. && score > NNWPs_Barrel_[4]) pass = true;  
        else if(seedEt>50. && seedEt<=60. && score > NNWPs_Barrel_[5]) pass = true;   
        else if(seedEt>60. && seedEt<=70. && score > NNWPs_Barrel_[6]) pass = true;  
        else if(seedEt>70. && seedEt<=80. && score > NNWPs_Barrel_[7]) pass = true;  
        else if(seedEt>80. && seedEt<=90. && score > NNWPs_Barrel_[8]) pass = true;  
        else if(seedEt>90. && score > NNWPs_Barrel_[9]) pass = true; 

     }else if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_ENDCAP){
         
        etaWindow = etaWindows_Endcap_[0] + etaWindows_Endcap_[1]*fabs(seedEta) + etaWindows_Endcap_[2]*fabs(seedEta)*fabs(seedEta);  
        phiWindow = phiWindows_Endcap_[0] + phiWindows_Endcap_[1]*fabs(seedEta) + phiWindows_Endcap_[2]*fabs(seedEta)*fabs(seedEta);
        if(deta < etaWindow && dphi < phiWindow) isIn = true;
        
        if(seedEt>0. && seedEt<=10. && score > NNWPs_Endcap_[0]) pass = true;  
        else if(seedEt>10. && seedEt<=20. && score > NNWPs_Endcap_[1]) pass = true;  
        else if(seedEt>20. && seedEt<=30. && score > NNWPs_Endcap_[2]) pass = true;  
        else if(seedEt>30. && seedEt<=40. && score > NNWPs_Endcap_[3]) pass = true;  
        else if(seedEt>40. && seedEt<=50. && score > NNWPs_Endcap_[4]) pass = true;  
        else if(seedEt>50. && seedEt<=60. && score > NNWPs_Endcap_[5]) pass = true;   
        else if(seedEt>60. && seedEt<=70. && score > NNWPs_Endcap_[6]) pass = true;  
        else if(seedEt>70. && seedEt<=80. && score > NNWPs_Endcap_[7]) pass = true;  
        else if(seedEt>80. && seedEt<=90. && score > NNWPs_Endcap_[8]) pass = true;  
        else if(seedEt>90. && score > NNWPs_Endcap_[9]) pass = true; 
     }
     return (isIn && pass);       
  } 

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

  void DeepSC::DeepSCID(const reco::SuperCluster& sc, 
			    int & nclusters, 
			    float & EoutsideDeepSC,
                            const CaloTopology *topology, 
                            const CaloSubdetectorGeometry* ebGeom, 
                            const CaloSubdetectorGeometry* eeGeom,
                            const EcalRecHitCollection *recHitsEB, 
                            const EcalRecHitCollection *recHitsEE) {
    DeepSCID(sc.clustersBegin(),sc.clustersEnd(),nclusters,EoutsideDeepSC,topology,ebGeom,eeGeom,recHitsEB,recHitsEE);
  }
  
  void DeepSC::DeepSCID(const CaloClusterPtrVector& clusters, 
			    int & nclusters, 
			    float & EoutsideDeepSC,
                            const CaloTopology *topology, 
                            const CaloSubdetectorGeometry* ebGeom, 
                            const CaloSubdetectorGeometry* eeGeom, 
                            const EcalRecHitCollection *recHitsEB, 
                            const EcalRecHitCollection *recHitsEE) {    
    DeepSCID(clusters.begin(),clusters.end(),nclusters,EoutsideDeepSC,topology,ebGeom,eeGeom,recHitsEB,recHitsEE);
  }
  
  void DeepSC::DeepSCID(const std::vector<const CaloCluster*>& clusters, 
			    int & nclusters, 
			    float & EoutsideDeepSC,   
                            const CaloTopology *topology, 
                            const CaloSubdetectorGeometry* ebGeom, 
                            const CaloSubdetectorGeometry* eeGeom, 
                            const EcalRecHitCollection *recHitsEB, 
                            const EcalRecHitCollection *recHitsEE) {
    DeepSCID(clusters.cbegin(),clusters.cend(),nclusters,EoutsideDeepSC,topology,ebGeom,eeGeom,recHitsEB,recHitsEE);
  }

  template<class RandomAccessPtrIterator>
  void DeepSC::DeepSCID(const RandomAccessPtrIterator& begin, 
			    const RandomAccessPtrIterator& end,
			    int & nclusters, 
			    float & EoutsideDeepSC,
                            const CaloTopology *topology, 
                            const CaloSubdetectorGeometry* ebGeom, 
                            const CaloSubdetectorGeometry* eeGeom,
                            const EcalRecHitCollection *recHitsEB, 
                            const EcalRecHitCollection *recHitsEE) {    
    nclusters = 0;
    EoutsideDeepSC = 0;
    
    unsigned int ncl = end-begin;
    if(!ncl) return;
    
    //loop over all clusters to find the one with highest energy
    RandomAccessPtrIterator icl = begin;
    RandomAccessPtrIterator clmax = end;
    float emax = 0;
    for( ; icl != end; ++icl){
      const float e = (*icl)->energy();
      if(e > emax){
	emax = e;
	clmax = icl;
      }
    }
    
    if(end == clmax) return;
    
    bool inDeep = false;
    icl = begin;
    for( ; icl != end; ++icl ){
      //auto seed = *(*clmax); 
      //auto cluster = *(*icl);  
      reco::CaloCluster* seed = new reco::CaloCluster(*(*clmax));
      reco::CaloCluster* cluster = new reco::CaloCluster(*(*icl));   
      inDeep=InSuperCluster(seed, cluster, topology, ebGeom, eeGeom, recHitsEB, recHitsEE);
      
      nclusters += (int)!inDeep;
      EoutsideDeepSC += (!inDeep)*((*icl)->energy()); 
    }
  }
  
  void DeepSC::DeepSCClust(const std::vector<CaloCluster>& clusters,
                           const CaloTopology *topology, 
                           const CaloSubdetectorGeometry* ebGeom, 
                           const CaloSubdetectorGeometry* eeGeom,
                           const EcalRecHitCollection *recHitsEB, 
                           const EcalRecHitCollection *recHitsEE,  
			   std::vector<unsigned int>& insideDeep, 
			   std::vector<unsigned int>& outsideDeep){  
    unsigned int ncl = clusters.size();
    if(!ncl) return;
    
    //loop over all clusters to find the one with highest energy
    float emax = 0;
    int imax = -1;
    for(unsigned int i=0; i<ncl; ++i){
      float e = (clusters[i]).energy();
      if(e > emax){
	emax = e;
	imax = i;
      }
    }
    
    if(imax<0) return;
 
    for(unsigned int k=0; k<ncl; k++){
      
      bool inDeep=InSuperCluster(&clusters[imax], &clusters[k], topology, ebGeom, eeGeom, recHitsEB, recHitsEE);

      //return indices of Clusters outside the DeepSC
      if (!(inDeep)){
	outsideDeep.push_back(k);
      }
      else{//return indices of Clusters inside the DeepSC
	insideDeep.push_back(k);
      }
    }
  }
  
  void DeepSC::FillDeepSCVar(const std::vector<CaloCluster>& clusters, const CaloTopology *topology, const CaloSubdetectorGeometry* ebGeom, const CaloSubdetectorGeometry* eeGeom, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE){
    Energy_In_DeepSC_=0;
    Energy_Outside_DeepSC_=0;
    LowestClusterEInDeepSC_=0;
    excluded_=0;
    included_=0;
    std::multimap<float, unsigned int>OrderedClust;
    std::vector<unsigned int> insideDeep;
    std::vector<unsigned int> outsideDeep;
    DeepSCClust(clusters, topology, ebGeom, eeGeom, recHitsEB, recHitsEE, insideDeep, outsideDeep);
    included_=insideDeep.size(); excluded_=outsideDeep.size();
    for(unsigned int i=0; i<insideDeep.size(); ++i){
      unsigned int index=insideDeep[i];
      Energy_In_DeepSC_=clusters[index].energy()+Energy_In_DeepSC_;
      OrderedClust.insert(make_pair(clusters[index].energy(), index));
    }
    for(unsigned int i=0; i<outsideDeep.size(); ++i){
      unsigned int index=outsideDeep[i];
      Energy_Outside_DeepSC_=clusters[index].energy()+Energy_Outside_DeepSC_;
      Et_Outside_DeepSC_=clusters[index].energy()*sin(clusters[index].position().theta())
	+Et_Outside_DeepSC_;
    }
    std::multimap<float, unsigned int>::iterator it;
    it=OrderedClust.begin();
    unsigned int lowEindex=(*it).second; 
    LowestClusterEInDeepSC_=clusters[lowEindex].energy();
    
  }
}

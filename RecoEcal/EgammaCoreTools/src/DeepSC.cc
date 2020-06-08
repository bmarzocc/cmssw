#include "RecoEcal/EgammaCoreTools/interface/DeepSC.h"
#include "TMath.h"
#include "TVector2.h"
#include <cmath>
using namespace std;

namespace reco {  
  
  DeepSC::DeepSC()
  {
     session_NN_Barrel_ = tensorflow::createSession(tensorflow::loadGraphDef(edm::FileInPath("RecoEcal/EgammaCoreTools/data/deepSC_model_v14_finalscore_v2_EB.pb").fullPath().c_str()));
     session_NN_Endcap_ = tensorflow::createSession(tensorflow::loadGraphDef(edm::FileInPath("RecoEcal/EgammaCoreTools/data/deepSC_model_v14_finalscore_v6_EE.pb").fullPath().c_str()));
     x_mean_Barrel_ = std::vector<double>({.60707421e-03, -1.73794814e-02,  0.00000000e+00,  3.56214761e+01,
2.49357538e+01,  3.11332870e-04,  1.98443393e-06,  7.73690247e-01,
5.41722921e-01,  1.00503862e+00,  8.74048387e-03,  2.70028418e-07,
1.15299516e-02,  4.90702664e-01,  4.67579586e+00,  6.67593469e-03,
8.21538199e-03,  1.62043356e+00,  7.03778100e-03,  1.82998814e-07,
7.76987369e-03,  7.07775815e-01,  1.41931161e+00,  2.88880743e-03,
2.77690441e-03});
     x_std_Barrel_ = std::vector<double>({8.70911527e-01, 1.80781524e+00, 1.00000000e+00, 4.43159827e+01,
2.97445738e+01, 1.09634038e-01, 3.31674900e-01, 1.93701364e+00,
1.23390516e+00, 8.53776838e-02, 4.78313451e-03, 6.11027544e-05,
5.66862122e-03, 4.01735208e-01, 5.49507876e+00, 2.94984302e-03,
4.01382913e-03, 6.81211508e+00, 8.17118200e-03, 1.08280997e-04,
8.54969802e-03, 1.42626514e+00, 1.70514757e+00, 2.06380546e-03,
2.74318441e-03});
     x_mean_Endcap_ = std::vector<double>({3.56700000e-03, -5.61881480e-04,  5.42432152e-03,  1.23144162e+02,
2.80286361e+01, -3.77895168e-03, -2.52080067e-04,  2.83572379e+00,
6.40789285e-01,  1.01572938e+00,  2.75843990e-02, -1.78424510e-05,
3.32103293e-02,  4.89525266e-01,  6.47583705e+00,  1.50177114e-02,
1.67321469e-02,  2.79373242e+00,  3.08011135e-02,  6.28599546e-05,
3.31705605e-02,  2.61118938e-01,  1.72031559e+00,  4.72343147e-03,
3.88800860e-03});
     x_std_Endcap_ = std::vector<double>({2.12599275e+00, 1.81089171e+00, 9.99985288e-01, 1.46843692e+02,
3.01028994e+01, 1.12075884e-01, 3.27203361e-01, 3.59426937e+00,
9.23646139e-01, 1.90638817e-01, 1.25780060e-02, 3.53072883e-04,
1.50233813e-02, 4.06099128e-01, 6.06635914e+00, 8.35589745e-03,
1.01933224e-02, 1.01608434e+01, 2.27313472e-02, 5.72635095e-04,
2.39377164e-02, 3.07407340e+00, 2.33897667e+00, 7.15246184e-03,
8.16679482e-03}); 
     NNinput_string_Barrel_ = std::string("dense_12_input:0");
     NNoutput_string_Barrel_ = std::string("dense_15/Sigmoid:0");
     NNinput_string_Endcap_ = std::string("dense_9_input:0");
     NNoutput_string_Endcap_ = std::string("dense_12/Sigmoid:0");
     NNWPs_Barrel_ = std::vector<double>({0.75,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85});
     NNWPs_Endcap_ = std::vector<double>({0.75,0.90,0.85,0.85,0.85,0.80,0.85,0.85,0.95,0.90});  
     etaWindows_Barrel_ = std::vector<double>({0.2,0.,0.}); // p0 + p1*seedEta + p2*seedEta*dEta
     phiWindows_Barrel_ = std::vector<double>({0.6,0.,0.}); // p0 + p1*seedEta + p2*seedEta*dEta   
     etaWindows_Endcap_ = std::vector<double>({0.2,0.,0.}); // p0 + p1*seedEta + p2*seedEta*dEta
     phiWindows_Endcap_ = std::vector<double>({0.6,0.,0.}); // p0 + p1*seedEta + p2*seedEta*dEta
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
    NNclusterVars_.resize(25);
    float e1=1.; 
    float e4=0.;

    double zSide = 0;
    if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_ENDCAP && seed->eta()<0.) zSide = -1.;   
    if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_ENDCAP && seed->eta()>0.) zSide = +1.; 

    NNclusterVars_[0] = seed->eta();
    NNclusterVars_[1] = seed->phi();
    NNclusterVars_[2] = zSide;
    NNclusterVars_[3] = seed->energy();
    NNclusterVars_[4] = seed->energy()/TMath::CosH(seed->eta());
    NNclusterVars_[5] = DeltaEta(seed->eta(), cluster->eta());  
    NNclusterVars_[6] = DeltaPhi(seed->phi(), cluster->phi()); 
    NNclusterVars_[7] = cluster->energy();
    NNclusterVars_[8] = cluster->energy()/TMath::CosH(cluster->eta());
    
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
       NNclusterVars_[9] = noZS::EcalClusterTools::e3x3(seedBC_, recHitsEB, topology)/seed->energy(); //r9
       NNclusterVars_[10] = sqrt(full5x5_locCov_[0]); //sigmaietaieta
       NNclusterVars_[11] = full5x5_locCov_[1]; //sigmaietaiphi
       NNclusterVars_[12] = (!edm::isFinite(full5x5_locCov_[2])) ? 0. : sqrt(full5x5_locCov_[2]); //sigmaiphiiphi     
       NNclusterVars_[13] = 1.-e4/e1; //swiss_cross      
       NNclusterVars_[14] = seed->hitsAndFractions().size(); //nXtals 
       NNclusterVars_[15] = widths_.first; //etaWidth 
       NNclusterVars_[16] = widths_.second; //phiWidth 
    }else if(PFLayer::fromCaloID(seed->caloID()) == PFLayer::ECAL_ENDCAP){
       full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(seedBC_, recHitsEE, topology);
       widths_ = ComputeCovariances(seedBC_, recHitsEE, eeGeom); 
       e1 = noZS::EcalClusterTools::eMax(seedBC_, recHitsEE);
       e4 = noZS::EcalClusterTools::eTop(seedBC_, recHitsEE, topology) +
            noZS::EcalClusterTools::eRight(seedBC_, recHitsEE, topology) +
            noZS::EcalClusterTools::eBottom(seedBC_, recHitsEE, topology) +
            noZS::EcalClusterTools::eLeft(seedBC_, recHitsEE, topology);     
       NNclusterVars_[9] = noZS::EcalClusterTools::e3x3(seedBC_, recHitsEE, topology)/seed->energy(); //r9
       NNclusterVars_[10] = sqrt(full5x5_locCov_[0]); //sigmaietaieta
       NNclusterVars_[11] = full5x5_locCov_[1]; //sigmaietaiphi
       NNclusterVars_[12] = (!edm::isFinite(full5x5_locCov_[2])) ? 0. : sqrt(full5x5_locCov_[2]); //sigmaiphiiphi     
       NNclusterVars_[13] = 1.-e4/e1; //swiss_cross      
       NNclusterVars_[14] = seed->hitsAndFractions().size(); //nXtals
       NNclusterVars_[15] = widths_.first; //etaWidth 
       NNclusterVars_[16] = widths_.second; //phiWidth    
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
       NNclusterVars_[17] = noZS::EcalClusterTools::e3x3(caloBC_, recHitsEB, topology)/cluster->energy(); //r9
       NNclusterVars_[18] = sqrt(full5x5_locCov_[0]); //sigmaietaieta
       NNclusterVars_[19] = full5x5_locCov_[1]; //sigmaietaiphi
       NNclusterVars_[20] = (!edm::isFinite(full5x5_locCov_[2])) ? 0. : sqrt(full5x5_locCov_[2]); //sigmaiphiiphi     
       NNclusterVars_[21] = 1.-e4/e1; //swiss_cross      
       NNclusterVars_[22] = cluster->hitsAndFractions().size(); //nXtals 
       NNclusterVars_[23] = widths_.first; //etaWidth 
       NNclusterVars_[24] = widths_.second; //phiWidth    
    }else if(PFLayer::fromCaloID(cluster->caloID()) == PFLayer::ECAL_ENDCAP){
       full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(caloBC_, recHitsEE, topology);
       widths_ = ComputeCovariances(caloBC_, recHitsEE, eeGeom);  
       e1 = noZS::EcalClusterTools::eMax(caloBC_, recHitsEE);
       e4 = noZS::EcalClusterTools::eTop(caloBC_, recHitsEE, topology) +
            noZS::EcalClusterTools::eRight(caloBC_, recHitsEE, topology) +
            noZS::EcalClusterTools::eBottom(caloBC_, recHitsEE, topology) +
            noZS::EcalClusterTools::eLeft(caloBC_, recHitsEE, topology);
       NNclusterVars_[17] = noZS::EcalClusterTools::e3x3(caloBC_, recHitsEE, topology)/cluster->energy(); //r9
       NNclusterVars_[18] = sqrt(full5x5_locCov_[0]); //sigmaietaieta
       NNclusterVars_[19] = full5x5_locCov_[1]; //sigmaietaiphi
       NNclusterVars_[20] = (!edm::isFinite(full5x5_locCov_[2])) ? 0. : sqrt(full5x5_locCov_[2]); //sigmaiphiiphi     
       NNclusterVars_[21] = 1.-e4/e1; //swiss_cross      
       NNclusterVars_[22] = cluster->hitsAndFractions().size(); //nXtals
       NNclusterVars_[23] = widths_.first; //etaWidth 
       NNclusterVars_[24] = widths_.second; //phiWidth    
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

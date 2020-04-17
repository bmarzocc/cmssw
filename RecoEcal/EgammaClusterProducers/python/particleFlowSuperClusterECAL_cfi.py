import FWCore.ParameterSet.Config as cms
import os

particleFlowSuperClusterECALBox = cms.EDProducer(
    "PFECALSuperClusterProducer",
    # verbosity 
    verbose = cms.untracked.bool(False),
    #clustering type: "Box" or "Mustache"
    ClusteringType = cms.string("Box"),
    #energy weighting: "Raw", "CalibratedNoPS", "CalibratedTotal"
    EnergyWeight = cms.string("Raw"),

    #this overrides both dphi cuts below if true!
    useDynamicDPhiWindow = cms.bool(False),
    
    #PFClusters collection
    PFClusters = cms.InputTag("particleFlowClusterECAL"),
    ESAssociation = cms.InputTag("particleFlowClusterECAL"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),    
                                                  
    PFBasicClusterCollectionBarrel = cms.string("particleFlowBasicClusterECALBarrel"),                                       
    PFSuperClusterCollectionBarrel = cms.string("particleFlowSuperClusterECALBarrel"),
    PFBasicClusterCollectionEndcap = cms.string("particleFlowBasicClusterECALEndcap"),                                       
    PFSuperClusterCollectionEndcap = cms.string("particleFlowSuperClusterECALEndcap"),
    PFBasicClusterCollectionPreshower = cms.string("particleFlowBasicClusterECALPreshower"),
    PFSuperClusterCollectionEndcapWithPreshower = cms.string("particleFlowSuperClusterECALEndcapWithPreshower"),                                          

    #use preshower ?
    use_preshower = cms.bool(True),

    # are the seed thresholds Et or Energy?
    seedThresholdIsET = cms.bool(True),

    # regression setup
    useRegression = cms.bool(False), #regressions are mustache only
    regressionConfig = cms.PSet(
       regressionKeyEB = cms.string('pfscecal_EBCorrection_offline_v2'),
       uncertaintyKeyEB = cms.string('pfscecal_EBUncertainty_offline_v2'),
       regressionKeyEE = cms.string('pfscecal_EECorrection_offline_v2'),
       uncertaintyKeyEE = cms.string('pfscecal_EEUncertainty_offline_v2'),
       vertexCollection = cms.InputTag("offlinePrimaryVertices"),
       ecalRecHitsEB = cms.InputTag('ecalRecHit','EcalRecHitsEB'),
       ecalRecHitsEE = cms.InputTag('ecalRecHit','EcalRecHitsEE'),
       applySigmaIetaIphiBug = cms.bool(False)
       ),

    #threshold for final SuperCluster Et
    thresh_SCEt = cms.double(4.0),    
    
    # threshold in ECAL
    thresh_PFClusterSeedBarrel = cms.double(3.0),
    thresh_PFClusterBarrel = cms.double(0.5),

    thresh_PFClusterSeedEndcap = cms.double(5.0),
    thresh_PFClusterEndcap = cms.double(0.5),

    # window width in ECAL
    phiwidth_SuperClusterBarrel = cms.double(0.28),
    etawidth_SuperClusterBarrel = cms.double(0.04),

    phiwidth_SuperClusterEndcap = cms.double(0.28),
    etawidth_SuperClusterEndcap = cms.double(0.04),

    # turn on merging of the seed cluster to its nearest neighbors
    # that share a rechit
    doSatelliteClusterMerge = cms.bool(False),
    satelliteClusterSeedThreshold = cms.double(50.0),
    satelliteMajorityFraction = cms.double(0.5),
    dropUnseedable = cms.bool(False),
    #thresh_PFClusterMustacheOutBarrel = cms.double(0.),
    #thresh_PFClusterMustacheOutEndcap = cms.double(0.), 
   
    #corrections  
    applyCrackCorrections = cms.bool(False),  
                                     
)

particleFlowSuperClusterECALMustache = cms.EDProducer(
    "PFECALSuperClusterProducer",
    # verbosity 
    verbose = cms.untracked.bool(False),
    #clustering type: "Box" or "Mustache"
    ClusteringType = cms.string("Mustache"),
    #energy weighting: "Raw", "CalibratedNoPS", "CalibratedTotal"
    EnergyWeight = cms.string("Raw"),

    #this overrides both dphi cuts below if true!
    useDynamicDPhiWindow = cms.bool(True), 
                                              
    #PFClusters collection
    PFClusters = cms.InputTag("particleFlowClusterECAL"),
    ESAssociation = cms.InputTag("particleFlowClusterECAL"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
                                              
    PFBasicClusterCollectionBarrel = cms.string("particleFlowBasicClusterECALBarrel"),                                       
    PFSuperClusterCollectionBarrel = cms.string("particleFlowSuperClusterECALBarrel"),
    PFBasicClusterCollectionEndcap = cms.string("particleFlowBasicClusterECALEndcap"),                                       
    PFSuperClusterCollectionEndcap = cms.string("particleFlowSuperClusterECALEndcap"),
    PFBasicClusterCollectionPreshower = cms.string("particleFlowBasicClusterECALPreshower"),
    PFSuperClusterCollectionEndcapWithPreshower = cms.string("particleFlowSuperClusterECALEndcapWithPreshower"),                                          

    # are the seed thresholds Et or Energy?
    seedThresholdIsET = cms.bool(True),
    # regression setup
    useRegression = cms.bool(True),
    regressionConfig = cms.PSet(
       regressionKeyEB = cms.string('pfscecal_EBCorrection_offline_v2'),
       uncertaintyKeyEB = cms.string('pfscecal_EBUncertainty_offline_v2'),
       regressionKeyEE = cms.string('pfscecal_EECorrection_offline_v2'),
       uncertaintyKeyEE = cms.string('pfscecal_EEUncertainty_offline_v2'),
       vertexCollection = cms.InputTag("offlinePrimaryVertices"),
       ecalRecHitsEB = cms.InputTag('ecalRecHit','EcalRecHitsEB'),
       ecalRecHitsEE = cms.InputTag('ecalRecHit','EcalRecHitsEE'),
       applySigmaIetaIphiBug = cms.bool(False)
       ),
       
    #threshold for final SuperCluster Et
    thresh_SCEt = cms.double(4.0),
    
    # threshold in ECAL
    thresh_PFClusterSeedBarrel = cms.double(1.0),
    thresh_PFClusterBarrel = cms.double(0.0),

    thresh_PFClusterSeedEndcap = cms.double(1.0),
    thresh_PFClusterEndcap = cms.double(0.0),

    # window width in ECAL ( these don't mean anything for Mustache )
    phiwidth_SuperClusterBarrel = cms.double(0.6),
    etawidth_SuperClusterBarrel = cms.double(0.04),

    phiwidth_SuperClusterEndcap = cms.double(0.6),
    etawidth_SuperClusterEndcap = cms.double(0.04),

    # threshold in preshower
    thresh_PFClusterES = cms.double(0.),           

    # turn on merging of the seed cluster to its nearest neighbors
    # that share a rechit
    doSatelliteClusterMerge = cms.bool(False),
    satelliteClusterSeedThreshold = cms.double(50.0),
    satelliteMajorityFraction = cms.double(0.5),
    dropUnseedable = cms.bool(False),
    #thresh_PFClusterMustacheOutBarrel = cms.double(0.),
    #thresh_PFClusterMustacheOutEndcap = cms.double(0.), 
    # old mustache parameters, before optimization 2019-2020
    p00 = cms.double(-0.107537),
    p01 = cms.double(0.590969),
    p02 = cms.double(-0.076494),
    p10 = cms.double(-0.0268843),
    p11 = cms.double(0.147742),
    p12 = cms.double(-0.0191235),
    w00 = cms.double(-0.00571429),
    w01 = cms.double(-0.002),
    w10 = cms.double(0.0135714),
    w11 = cms.double(0.001),
    yoffsetEB = cms.double(7.151e-02),
    scaleEB = cms.double(5.656e-01),
    xoffsetEB = cms.double(2.931e-01),
    widthEB = cms.double(2.976e-01),
    yoffsetEE_0 = cms.double(5.058e-02),
    scaleEE_0 = cms.double(7.131e-01),
    xoffsetEE_0 = cms.double(1.668e-02),
    widthEE_0 = cms.double(4.114e-01),
    yoffsetEE_1 = cms.double(-9.913e-02),
    scaleEE_1 = cms.double(4.404e+01),
    xoffsetEE_1 = cms.double(-5.326e+00),
    widthEE_1 = cms.double(1.184e+00),
    yoffsetEE_2 = cms.double(-6.346e-01),
    scaleEE_2 = cms.double(1.317e+01),
    xoffsetEE_2 = cms.double(-7.037e+00),
    widthEE_2 = cms.double(2.836e+00),

    #corrections
    applyCrackCorrections = cms.bool(False),
               
)

particleFlowSuperClusterECALMustacheNewParams = cms.EDProducer(
    "PFECALSuperClusterProducer",
    # verbosity 
    verbose = cms.untracked.bool(False),
    #clustering type: "Box" or "Mustache"
    ClusteringType = cms.string("Mustache"),
    #energy weighting: "Raw", "CalibratedNoPS", "CalibratedTotal"
    EnergyWeight = cms.string("Raw"),

    #this overrides both dphi cuts below if true!
    useDynamicDPhiWindow = cms.bool(True), 
                                              
    #PFClusters collection
    PFClusters = cms.InputTag("particleFlowClusterECAL"),
    ESAssociation = cms.InputTag("particleFlowClusterECAL"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
                                              
    PFBasicClusterCollectionBarrel = cms.string("particleFlowBasicClusterECALBarrelMustacheNewParams"),                                       
    PFSuperClusterCollectionBarrel = cms.string("particleFlowSuperClusterECALBarrelMustacheNewParams"),
    PFBasicClusterCollectionEndcap = cms.string("particleFlowBasicClusterECALEndcapMustacheNewParams"),                                       
    PFSuperClusterCollectionEndcap = cms.string("particleFlowSuperClusterECALEndcapMustacheNewParams"),
    PFBasicClusterCollectionPreshower = cms.string("particleFlowBasicClusterECALPreshowerMustacheNewParams"),
    PFSuperClusterCollectionEndcapWithPreshower = cms.string("particleFlowSuperClusterECALEndcapWithPreshowerMustacheNewParams"),                                          

    # are the seed thresholds Et or Energy?
    seedThresholdIsET = cms.bool(True),
    # regression setup
    useRegression = cms.bool(True),
    regressionConfig = cms.PSet(
       regressionKeyEB = cms.string('pfscecal_EBCorrection_offline_v2'),
       uncertaintyKeyEB = cms.string('pfscecal_EBUncertainty_offline_v2'),
       regressionKeyEE = cms.string('pfscecal_EECorrection_offline_v2'),
       uncertaintyKeyEE = cms.string('pfscecal_EEUncertainty_offline_v2'),
       vertexCollection = cms.InputTag("offlinePrimaryVertices"),
       ecalRecHitsEB = cms.InputTag('ecalRecHit','EcalRecHitsEB'),
       ecalRecHitsEE = cms.InputTag('ecalRecHit','EcalRecHitsEE'),
       applySigmaIetaIphiBug = cms.bool(False)
       ),
       
    #threshold for final SuperCluster Et
    thresh_SCEt = cms.double(4.0),
    
    # threshold in ECAL
    thresh_PFClusterSeedBarrel = cms.double(1.0),
    thresh_PFClusterBarrel = cms.double(0.0),

    thresh_PFClusterSeedEndcap = cms.double(1.0),
    thresh_PFClusterEndcap = cms.double(0.0),

    # window width in ECAL ( these don't mean anything for Mustache )
    phiwidth_SuperClusterBarrel = cms.double(0.6),
    etawidth_SuperClusterBarrel = cms.double(0.04),

    phiwidth_SuperClusterEndcap = cms.double(0.6),
    etawidth_SuperClusterEndcap = cms.double(0.04),

    # threshold in preshower
    thresh_PFClusterES = cms.double(0.),           

    # turn on merging of the seed cluster to its nearest neighbors
    # that share a rechit
    doSatelliteClusterMerge = cms.bool(False),
    satelliteClusterSeedThreshold = cms.double(50.0),
    satelliteMajorityFraction = cms.double(0.5),
    dropUnseedable = cms.bool(False),
    #thresh_PFClusterMustacheOutBarrel = cms.double(0.),
    #thresh_PFClusterMustacheOutEndcap = cms.double(0.), 
    p00 = cms.double(-0.107537),
    p01 = cms.double(0.590969),
    p02 = cms.double(-0.076494),
    p10 = cms.double(-0.0268843),
    p11 = cms.double(0.147742),
    p12 = cms.double(-0.0191235),
    w00 = cms.double(-0.00016325),
    w01 = cms.double(0.00330476),
    w10 = cms.double(0.00663276),
    w11 = cms.double(-0.00461881),
    yoffsetEB = cms.double(-0.0374196),
    scaleEB = cms.double(2.53626),
    xoffsetEB = cms.double(-1.18163),
    widthEB = cms.double(0.782166),
    yoffsetEE_0 = cms.double(0.083572),
    scaleEE_0 = cms.double(0.656847),
    xoffsetEE_0 = cms.double(0.0207699),
    widthEE_0 = cms.double(0.298995),
    yoffsetEE_1 = cms.double(-0.101296),
    scaleEE_1 = cms.double(32.2549),
    xoffsetEE_1 = cms.double(-5.5264),
    widthEE_1 = cms.double(1.29536),
    yoffsetEE_2 = cms.double(-0.42682),
    scaleEE_2 = cms.double(10.9146),
    xoffsetEE_2 = cms.double(-7.44059),
    widthEE_2 = cms.double(2.85026),                          

    #corrections
    applyCrackCorrections = cms.bool(False),
               
)

particleFlowSuperClusterECALDeepSC = cms.EDProducer(
    "PFECALSuperClusterProducer",
    # verbosity 
    verbose = cms.untracked.bool(False),
    #clustering type: "Box" or "Mustache"
    ClusteringType = cms.string("DeepSC"),
    #energy weighting: "Raw", "CalibratedNoPS", "CalibratedTotal"
    EnergyWeight = cms.string("Raw"),

    #this overrides both dphi cuts below if true!
    useDynamicDPhiWindow = cms.bool(False), 
                                              
    #PFClusters collection
    PFClusters = cms.InputTag("particleFlowClusterECAL"),
    ESAssociation = cms.InputTag("particleFlowClusterECAL"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
                                              
    PFBasicClusterCollectionBarrel = cms.string("particleFlowDeepBasicClusterECALBarrel"),                                       
    PFSuperClusterCollectionBarrel = cms.string("particleFlowDeepSuperClusterECALBarrel"),
    PFBasicClusterCollectionEndcap = cms.string("particleFlowDeepBasicClusterECALEndcap"),                                       
    PFSuperClusterCollectionEndcap = cms.string("particleFlowDeepSuperClusterECALEndcap"),
    PFBasicClusterCollectionPreshower = cms.string("particleFlowBDeepBasicClusterECALPreshower"),
    PFSuperClusterCollectionEndcapWithPreshower = cms.string("particleFlowDeepSuperClusterECALEndcapWithPreshower"),                                          

    # are the seed thresholds Et or Energy?
    seedThresholdIsET = cms.bool(True),
    # regression setup
    useRegression = cms.bool(True),
    regressionConfig = cms.PSet(
       regressionKeyEB = cms.string('pfscecal_EBCorrection_offline_v2'),
       uncertaintyKeyEB = cms.string('pfscecal_EBUncertainty_offline_v2'),
       regressionKeyEE = cms.string('pfscecal_EECorrection_offline_v2'),
       uncertaintyKeyEE = cms.string('pfscecal_EEUncertainty_offline_v2'),
       vertexCollection = cms.InputTag("offlinePrimaryVertices"),
       ecalRecHitsEB = cms.InputTag('ecalRecHit','EcalRecHitsEB'),
       ecalRecHitsEE = cms.InputTag('ecalRecHit','EcalRecHitsEE'),
       applySigmaIetaIphiBug = cms.bool(False)
       ),
       
    #threshold for final SuperCluster Et
    thresh_SCEt = cms.double(4.0),
    
    # threshold in ECAL
    thresh_PFClusterSeedBarrel = cms.double(1.0),
    thresh_PFClusterBarrel = cms.double(0.0),

    thresh_PFClusterSeedEndcap = cms.double(1.0),
    thresh_PFClusterEndcap = cms.double(0.0),

    # window width in ECAL ( these don't mean anything for Mustache )
    phiwidth_SuperClusterBarrel = cms.double(0.7),
    etawidth_SuperClusterBarrel = cms.double(0.3),

    phiwidth_SuperClusterEndcap = cms.double(0.7),
    etawidth_SuperClusterEndcap = cms.double(0.3),

    # threshold in preshower
    thresh_PFClusterES = cms.double(0.),           

    # turn on merging of the seed cluster to its nearest neighbors
    # that share a rechit
    doSatelliteClusterMerge = cms.bool(False),
    satelliteClusterSeedThreshold = cms.double(50.0),
    satelliteMajorityFraction = cms.double(0.5),
    dropUnseedable = cms.bool(False),
    #thresh_PFClusterMustacheOutBarrel = cms.double(0.),
    #thresh_PFClusterMustacheOutEndcap = cms.double(0.), 

    #corrections
    applyCrackCorrections = cms.bool(False),

    #DeepSC NN showerShapes v2 inputs
    NNfileName = cms.untracked.string("/afs/cern.ch/work/b/bmarzocc/DeepClustering/deepcluster/models/v2_showervars/model1_showervars_100_80_50.pb"),
    NNevaluationParams = cms.untracked.vstring("dense_5_input:0","dense_8/Sigmoid:0"),
    NNnormalizationParams_mean = cms.untracked.vdouble(4.90602126e-03, -5.35902050e-03,  1.32942831e-03, -3.21680264e-03,
1.34150907e-04,  9.92061190e+01,  1.92094031e+00,  3.38077598e+00,
1.48469587e-02,  1.50669136e-05,  1.68650608e-02,  8.44421380e-01,
1.84096328e+00),
    NNnormalizationParams_std = cms.untracked.vdouble(1.36047253e+00, 1.80968595e+00, 5.36392643e-01, 1.55475234e-01,
3.67102008e-01, 1.06488331e+02, 3.98013560e+00, 1.35784224e+01,
1.75594624e-02, 3.13839208e-04, 1.90925831e-02, 3.19459787e-01,
2.70549407e+00), #seedEta, seedPhi, seedZside, clusterDeta, clusterDphi, seedEnergy, clusterEnergy, clusterR9, clusterSigmaIetaIeta, clusterSigmaIetaIphi, clusterSigmaIphiIphi, clusterSwissCross, clusterNXtals 
    NNscoreWP = cms.untracked.double(0.886)

    #DeepSC NN showerShapes v3 inputs
    #NNfileName = cms.untracked.string("/afs/cern.ch/work/b/bmarzocc/DeepClustering/deepcluster/models/v3_showervars/v3_newscore/model1_newscore_100_80_50.pb"),
    #NNevaluationParams = cms.untracked.vstring("dense_5_input:0","dense_8/Sigmoid:0"),
    #NNnormalizationParams_mean = cms.untracked.vdouble(7.55142242e-03, -1.72542329e-02,  2.20115368e-03, -2.16739912e-03,
    #2.41509171e-05,  9.15054607e+01,  2.29035380e+00,  2.35499840e+00,
    #1.36106254e-02,  1.28558875e-05,  1.58615069e-02,  8.13872733e-01,
    #1.91021365e+00),
    #NNnormalizationParams_std = cms.untracked.vdouble(1.27738454e+00, 1.81072314e+00, 4.75175842e-01, 1.50758418e-01,
    #3.57268056e-01, 9.80173185e+01, 4.50022007e+00, 8.11772595e+00,
    #1.63615644e-02, 2.86207324e-04, 1.74513492e-02, 3.39357575e-01,
    #2.80740256e+00), #seedEta, seedPhi, seedZside, clusterDeta, clusterDphi, seedEnergy, clusterEnergy, clusterR9, clusterSigmaIetaIeta, clusterSigmaIetaIphi, clusterSigmaIphiIphi, clusterSwissCross, clusterNXtals 
    #NNscoreWP = cms.untracked.double(0.898)
                                              
)

particleFlowSuperClusterECALDeepSCLWP = cms.EDProducer(
    "PFECALSuperClusterProducer",
    # verbosity 
    verbose = cms.untracked.bool(False),
    #clustering type: "Box" or "Mustache"
    ClusteringType = cms.string("DeepSC"),
    #energy weighting: "Raw", "CalibratedNoPS", "CalibratedTotal"
    EnergyWeight = cms.string("Raw"),

    #this overrides both dphi cuts below if true!
    useDynamicDPhiWindow = cms.bool(False), 
                                              
    #PFClusters collection
    PFClusters = cms.InputTag("particleFlowClusterECAL"),
    ESAssociation = cms.InputTag("particleFlowClusterECAL"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
                                              
    PFBasicClusterCollectionBarrel = cms.string("particleFlowDeepBasicClusterLWPECALBarrel"),                                       
    PFSuperClusterCollectionBarrel = cms.string("particleFlowDeepSuperClusterLWPECALBarrel"),
    PFBasicClusterCollectionEndcap = cms.string("particleFlowDeepBasicClusterLWPECALEndcap"),                                       
    PFSuperClusterCollectionEndcap = cms.string("particleFlowDeepSuperClusterLWPECALEndcap"),
    PFBasicClusterCollectionPreshower = cms.string("particleFlowBDeepBasicClusterLWPECALPreshower"),
    PFSuperClusterCollectionEndcapWithPreshower = cms.string("particleFlowDeepSuperClusterLWPECALEndcapWithPreshower"),                                          

    # are the seed thresholds Et or Energy?
    seedThresholdIsET = cms.bool(True),
    # regression setup
    useRegression = cms.bool(True),
    regressionConfig = cms.PSet(
       regressionKeyEB = cms.string('pfscecal_EBCorrection_offline_v2'),
       uncertaintyKeyEB = cms.string('pfscecal_EBUncertainty_offline_v2'),
       regressionKeyEE = cms.string('pfscecal_EECorrection_offline_v2'),
       uncertaintyKeyEE = cms.string('pfscecal_EEUncertainty_offline_v2'),
       vertexCollection = cms.InputTag("offlinePrimaryVertices"),
       ecalRecHitsEB = cms.InputTag('ecalRecHit','EcalRecHitsEB'),
       ecalRecHitsEE = cms.InputTag('ecalRecHit','EcalRecHitsEE'),
       applySigmaIetaIphiBug = cms.bool(False)
       ),
       
    #threshold for final SuperCluster Et
    thresh_SCEt = cms.double(4.0),
    
    # threshold in ECAL
    thresh_PFClusterSeedBarrel = cms.double(1.0),
    thresh_PFClusterBarrel = cms.double(0.0),

    thresh_PFClusterSeedEndcap = cms.double(1.0),
    thresh_PFClusterEndcap = cms.double(0.0),

    # window width in ECAL ( these don't mean anything for Mustache )
    phiwidth_SuperClusterBarrel = cms.double(0.7),
    etawidth_SuperClusterBarrel = cms.double(0.3),

    phiwidth_SuperClusterEndcap = cms.double(0.7),
    etawidth_SuperClusterEndcap = cms.double(0.3),

    # threshold in preshower
    thresh_PFClusterES = cms.double(0.),           

    # turn on merging of the seed cluster to its nearest neighbors
    # that share a rechit
    doSatelliteClusterMerge = cms.bool(False),
    satelliteClusterSeedThreshold = cms.double(50.0),
    satelliteMajorityFraction = cms.double(0.5),
    dropUnseedable = cms.bool(False),
    #thresh_PFClusterMustacheOutBarrel = cms.double(0.),
    #thresh_PFClusterMustacheOutEndcap = cms.double(0.), 

    #corrections
    applyCrackCorrections = cms.bool(False),

    #DeepSC NN showerShapes v2 inputs
    NNfileName = cms.untracked.string("/afs/cern.ch/work/b/bmarzocc/DeepClustering/deepcluster/models/v2_showervars/model1_showervars_100_80_50.pb"),
    NNevaluationParams = cms.untracked.vstring("dense_5_input:0","dense_8/Sigmoid:0"),
    NNnormalizationParams_mean = cms.untracked.vdouble(4.90602126e-03, -5.35902050e-03,  1.32942831e-03, -3.21680264e-03,
1.34150907e-04,  9.92061190e+01,  1.92094031e+00,  3.38077598e+00,
1.48469587e-02,  1.50669136e-05,  1.68650608e-02,  8.44421380e-01,
1.84096328e+00),
    NNnormalizationParams_std = cms.untracked.vdouble(1.36047253e+00, 1.80968595e+00, 5.36392643e-01, 1.55475234e-01,
3.67102008e-01, 1.06488331e+02, 3.98013560e+00, 1.35784224e+01,
1.75594624e-02, 3.13839208e-04, 1.90925831e-02, 3.19459787e-01,
2.70549407e+00), #seedEta, seedPhi, seedZside, clusterDeta, clusterDphi, seedEnergy, clusterEnergy, clusterR9, clusterSigmaIetaIeta, clusterSigmaIetaIphi, clusterSigmaIphiIphi, clusterSwissCross, clusterNXtals 
    NNscoreWP = cms.untracked.double(0.726)

    #DeepSC NN showerShapes v3 inputs
    #NNfileName = cms.untracked.string("/afs/cern.ch/work/b/bmarzocc/DeepClustering/deepcluster/models/v3_showervars/v3_newscore/model1_newscore_100_80_50.pb"),
    #NNevaluationParams = cms.untracked.vstring("dense_5_input:0","dense_8/Sigmoid:0"),
    #NNnormalizationParams_mean = cms.untracked.vdouble(7.55142242e-03, -1.72542329e-02,  2.20115368e-03, -2.16739912e-03,
    #2.41509171e-05,  9.15054607e+01,  2.29035380e+00,  2.35499840e+00,
    #1.36106254e-02,  1.28558875e-05,  1.58615069e-02,  8.13872733e-01,
    #1.91021365e+00),
    #NNnormalizationParams_std = cms.untracked.vdouble(1.27738454e+00, 1.81072314e+00, 4.75175842e-01, 1.50758418e-01,
    #3.57268056e-01, 9.80173185e+01, 4.50022007e+00, 8.11772595e+00,
    #1.63615644e-02, 2.86207324e-04, 1.74513492e-02, 3.39357575e-01,
    #2.80740256e+00), #seedEta, seedPhi, seedZside, clusterDeta, clusterDphi, seedEnergy, clusterEnergy, clusterR9, clusterSigmaIetaIeta, clusterSigmaIetaIphi, clusterSigmaIphiIphi, clusterSwissCross, clusterNXtals 
    #NNscoreWP = cms.untracked.double(0.523)               
)

particleFlowSuperClusterECALDeepSCTWP = cms.EDProducer(
    "PFECALSuperClusterProducer",
    # verbosity 
    verbose = cms.untracked.bool(False),
    #clustering type: "Box" or "Mustache"
    ClusteringType = cms.string("DeepSC"),
    #energy weighting: "Raw", "CalibratedNoPS", "CalibratedTotal"
    EnergyWeight = cms.string("Raw"),

    #this overrides both dphi cuts below if true!
    useDynamicDPhiWindow = cms.bool(False), 
                                              
    #PFClusters collection
    PFClusters = cms.InputTag("particleFlowClusterECAL"),
    ESAssociation = cms.InputTag("particleFlowClusterECAL"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
                                              
    PFBasicClusterCollectionBarrel = cms.string("particleFlowDeepBasicClusterTWPECALBarrel"),                                       
    PFSuperClusterCollectionBarrel = cms.string("particleFlowDeepSuperClusterTWPECALBarrel"),
    PFBasicClusterCollectionEndcap = cms.string("particleFlowDeepBasicClusterTWPECALEndcap"),                                       
    PFSuperClusterCollectionEndcap = cms.string("particleFlowDeepSuperClusterTWPECALEndcap"),
    PFBasicClusterCollectionPreshower = cms.string("particleFlowBDeepBasicClusterTWPECALPreshower"),
    PFSuperClusterCollectionEndcapWithPreshower = cms.string("particleFlowDeepSuperClusterTWPECALEndcapWithPreshower"),                                          

    # are the seed thresholds Et or Energy?
    seedThresholdIsET = cms.bool(True),
    # regression setup
    useRegression = cms.bool(True),
    regressionConfig = cms.PSet(
       regressionKeyEB = cms.string('pfscecal_EBCorrection_offline_v2'),
       uncertaintyKeyEB = cms.string('pfscecal_EBUncertainty_offline_v2'),
       regressionKeyEE = cms.string('pfscecal_EECorrection_offline_v2'),
       uncertaintyKeyEE = cms.string('pfscecal_EEUncertainty_offline_v2'),
       vertexCollection = cms.InputTag("offlinePrimaryVertices"),
       ecalRecHitsEB = cms.InputTag('ecalRecHit','EcalRecHitsEB'),
       ecalRecHitsEE = cms.InputTag('ecalRecHit','EcalRecHitsEE'),
       applySigmaIetaIphiBug = cms.bool(False)
       ),
       
    #threshold for final SuperCluster Et
    thresh_SCEt = cms.double(4.0),
    
    # threshold in ECAL
    thresh_PFClusterSeedBarrel = cms.double(1.0),
    thresh_PFClusterBarrel = cms.double(0.0),

    thresh_PFClusterSeedEndcap = cms.double(1.0),
    thresh_PFClusterEndcap = cms.double(0.0),

    # window width in ECAL ( these don't mean anything for Mustache )
    phiwidth_SuperClusterBarrel = cms.double(0.7),
    etawidth_SuperClusterBarrel = cms.double(0.3),

    phiwidth_SuperClusterEndcap = cms.double(0.7),
    etawidth_SuperClusterEndcap = cms.double(0.3),

    # threshold in preshower
    thresh_PFClusterES = cms.double(0.),           

    # turn on merging of the seed cluster to its nearest neighbors
    # that share a rechit
    doSatelliteClusterMerge = cms.bool(False),
    satelliteClusterSeedThreshold = cms.double(50.0),
    satelliteMajorityFraction = cms.double(0.5),
    dropUnseedable = cms.bool(False),
    #thresh_PFClusterMustacheOutBarrel = cms.double(0.),
    #thresh_PFClusterMustacheOutEndcap = cms.double(0.), 

    #corrections
    applyCrackCorrections = cms.bool(False),

    #DeepSC NN showerShapes v2 inputs
    NNfileName = cms.untracked.string("/afs/cern.ch/work/b/bmarzocc/DeepClustering/deepcluster/models/v2_showervars/model1_showervars_100_80_50.pb"),
    NNevaluationParams = cms.untracked.vstring("dense_5_input:0","dense_8/Sigmoid:0"),
    NNnormalizationParams_mean = cms.untracked.vdouble(4.90602126e-03, -5.35902050e-03,  1.32942831e-03, -3.21680264e-03,
1.34150907e-04,  9.92061190e+01,  1.92094031e+00,  3.38077598e+00,
1.48469587e-02,  1.50669136e-05,  1.68650608e-02,  8.44421380e-01,
1.84096328e+00),
    NNnormalizationParams_std = cms.untracked.vdouble(1.36047253e+00, 1.80968595e+00, 5.36392643e-01, 1.55475234e-01,
3.67102008e-01, 1.06488331e+02, 3.98013560e+00, 1.35784224e+01,
1.75594624e-02, 3.13839208e-04, 1.90925831e-02, 3.19459787e-01,
2.70549407e+00), #seedEta, seedPhi, seedZside, clusterDeta, clusterDphi, seedEnergy, clusterEnergy, clusterR9, clusterSigmaIetaIeta, clusterSigmaIetaIphi, clusterSigmaIphiIphi, clusterSwissCross, clusterNXtals 
    NNscoreWP = cms.untracked.double(0.948)

    #DeepSC NN showerShapes v3 inputs
    #NNfileName = cms.untracked.string("/afs/cern.ch/work/b/bmarzocc/DeepClustering/deepcluster/models/v3_showervars/v3_newscore/model1_newscore_100_80_50.pb"),
    #NNevaluationParams = cms.untracked.vstring("dense_5_input:0","dense_8/Sigmoid:0"),
    #NNnormalizationParams_mean = cms.untracked.vdouble(7.55142242e-03, -1.72542329e-02,  2.20115368e-03, -2.16739912e-03,
    #2.41509171e-05,  9.15054607e+01,  2.29035380e+00,  2.35499840e+00,
    #1.36106254e-02,  1.28558875e-05,  1.58615069e-02,  8.13872733e-01,
    #1.91021365e+00),
    #NNnormalizationParams_std = cms.untracked.vdouble(1.27738454e+00, 1.81072314e+00, 4.75175842e-01, 1.50758418e-01,
    #3.57268056e-01, 9.80173185e+01, 4.50022007e+00, 8.11772595e+00,
    #1.63615644e-02, 2.86207324e-04, 1.74513492e-02, 3.39357575e-01,
    #2.80740256e+00), #seedEta, seedPhi, seedZside, clusterDeta, clusterDphi, seedEnergy, clusterEnergy, clusterR9, clusterSigmaIetaIeta, clusterSigmaIetaIphi, clusterSigmaIphiIphi, clusterSwissCross, clusterNXtals 
    #NNscoreWP = cms.untracked.double(0.948)                                                   
)

#define the default clustering type
particleFlowSuperClusterECAL = particleFlowSuperClusterECALMustache.clone()
particleFlowSuperClusterECALNewParams = particleFlowSuperClusterECALMustacheNewParams.clone() 
particleFlowDeepSuperClusterECAL = particleFlowSuperClusterECALDeepSC.clone()
particleFlowDeepSuperClusterLWPECAL = particleFlowSuperClusterECALDeepSCLWP.clone()
particleFlowDeepSuperClusterTWPECAL = particleFlowSuperClusterECALDeepSCTWP.clone()

from Configuration.Eras.Modifier_pp_on_AA_2018_cff import pp_on_AA_2018
pp_on_AA_2018.toModify(particleFlowSuperClusterECAL, useDynamicDPhiWindow = False)
pp_on_AA_2018.toModify(particleFlowSuperClusterECAL, phiwidth_SuperClusterBarrel = 0.20)
pp_on_AA_2018.toModify(particleFlowSuperClusterECAL, phiwidth_SuperClusterEndcap = 0.20)

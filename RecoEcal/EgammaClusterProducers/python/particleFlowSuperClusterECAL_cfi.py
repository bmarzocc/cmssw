import FWCore.ParameterSet.Config as cms

from RecoEcal.EgammaClusterProducers.particleFlowSuperClusterECALMustache_cfi import particleFlowSuperClusterECALMustache as _particleFlowSuperClusterECALMustache
from RecoEcal.EgammaClusterProducers.particleFlowSuperClusterDeepSC_cfi import particleFlowSuperClusterECALDeepSC as _particleFlowSuperClusterECALDeepSC

# define the default ECAL clustering (Mustache or Box)
#particleFlowSuperClusterECAL = _particleFlowSuperClusterECALMustache.clone()
particleFlowSuperClusterECAL = _particleFlowSuperClusterECALDeepSC.clone()

from Configuration.Eras.Modifier_pp_on_AA_2018_cff import pp_on_AA_2018
pp_on_AA_2018.toModify(particleFlowSuperClusterECAL, useDynamicDPhiWindow = False)
pp_on_AA_2018.toModify(particleFlowSuperClusterECAL, phiwidth_SuperClusterBarrel = 0.20)
pp_on_AA_2018.toModify(particleFlowSuperClusterECAL, phiwidth_SuperClusterEndcap = 0.20)

from Configuration.ProcessModifiers.egamma_lowPt_exclusive_cff import egamma_lowPt_exclusive
egamma_lowPt_exclusive.toModify(particleFlowSuperClusterECAL,
                                thresh_SCEt = 1.0,
                                thresh_PFClusterSeedBarrel = 0.5,
                                thresh_PFClusterSeedEndcap = 0.5)

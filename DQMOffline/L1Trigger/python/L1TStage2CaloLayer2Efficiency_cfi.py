import FWCore.ParameterSet.Config as cms
from DQMOffline.L1Trigger import L1TStage2CaloLayer2Offline_cfi as L1TStep1

variables = {
    'jet': L1TStep1.jetEfficiencyThresholds,
    'met': L1TStep1.metEfficiencyThresholds,
    'mht': L1TStep1.mhtEfficiencyThresholds,
    'ett': L1TStep1.ettEfficiencyThresholds,
    'htt': L1TStep1.httEfficiencyThresholds,
}

plots = {
    'jet': [
        "efficiencyJetEt_HB", "efficiencyJetEt_HE", "efficiencyJetEt_HF",
        "efficiencyJetEt_HB_HE"],
    'met': ['efficiencyMET', 'efficiencyETMHF'],
    'mht': ['efficiencyMHT'],
    'ett': ['efficiencyETT'],
    'htt': ['efficiencyHTT'],
}

allEfficiencyPlots = []
add_plot = allEfficiencyPlots.append
for variable, thresholds in variables.iteritems():
    for plot in plots[variable]:
        for threshold in thresholds:
            plotName = '{0}_threshold_{1}'.format(plot, threshold)
            add_plot(plotName)

from DQMOffline.L1Trigger.L1TEfficiencyHarvesting_cfi import l1tEfficiencyHarvesting
l1tStage2CaloLayer2Efficiency = l1tEfficiencyHarvesting.clone(
    plotCfgs=cms.untracked.VPSet(
        cms.untracked.PSet(
            numeratorDir=cms.untracked.string("L1T/L1TStage2CaloLayer2/efficiency_raw"),
            outputDir=cms.untracked.string("L1T/L1TStage2CaloLayer2"),
            numeratorSuffix=cms.untracked.string("_Num"),
            denominatorSuffix=cms.untracked.string("_Den"),
            plots=cms.untracked.vstring(allEfficiencyPlots)
        ),
        cms.untracked.PSet(
            numeratorDir=cms.untracked.string("L1TEMU/L1TStage2CaloLayer2/efficiency_raw"),
            outputDir=cms.untracked.string("L1TEMU/L1TStage2CaloLayer2"),
            numeratorSuffix=cms.untracked.string("_Num"),
            denominatorSuffix=cms.untracked.string("_Den"),
            plots=cms.untracked.vstring(allEfficiencyPlots)
        ),
    )
)

# modifications for the pp reference run
variables_HI = variables
variables_HI['jet'] = L1TStep1.jetEfficiencyThresholds_HI

allEfficiencyPlots_HI = []
add_plot = allEfficiencyPlots_HI.append
for variable, thresholds in variables_HI.iteritems():
    for plot in plots[variable]:
        for threshold in thresholds:
            plotName = '{0}_threshold_{1}'.format(plot, threshold)
            add_plot(plotName)

from Configuration.Eras.Modifier_ppRef_2017_cff import ppRef_2017
ppRef_2017.toModify(l1tStage2CaloLayer2Efficiency,
    plotCfgs = {
        0:dict(plots = allEfficiencyPlots_HI),
        1:dict(plots = allEfficiencyPlots_HI)
    }
)


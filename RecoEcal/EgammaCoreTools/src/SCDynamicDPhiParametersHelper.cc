// Implementation of the SC dynamic dPhi parameters interface

#include "RecoEcal/EgammaCoreTools/interface/SCDynamicDPhiParametersHelper.h"

#include <algorithm>
#include <utility>

using namespace reco;

SCDynamicDPhiParametersHelper::SCDynamicDPhiParametersHelper(EcalSCDynamicDPhiParameters &params)
    : parameters_(params) {}

SCDynamicDPhiParametersHelper::SCDynamicDPhiParametersHelper(EcalSCDynamicDPhiParameters &params,
                                                             const edm::ParameterSet &iConfig)
    : parameters_(params) {
  // dynamic dPhi parameters
  const auto dynamicDPhiPSets = iConfig.getParameter<std::vector<edm::ParameterSet>>("dynamicDPhiParameterSets");
  for (const auto &pSet : dynamicDPhiPSets) {
    EcalSCDynamicDPhiParameters::DynamicDPhiParameters params({pSet.getParameter<double>("eMin"),
                                                               pSet.getParameter<double>("etaMin"),
                                                               pSet.getParameter<double>("yoffset"),
                                                               pSet.getParameter<double>("scale"),
                                                               pSet.getParameter<double>("xoffset"),
                                                               pSet.getParameter<double>("width"),
                                                               pSet.getParameter<double>("saturation"),
                                                               pSet.getParameter<double>("cutoff")});
    addDynamicDPhiParameters(params);
    sortDynamicDPhiParametersCollection();
  }
}

void SCDynamicDPhiParametersHelper::addDynamicDPhiParameters(
    const EcalSCDynamicDPhiParameters::DynamicDPhiParameters &params) {
  parameters_.dynamicDPhiParametersCollection_.emplace_back(params);
}

void SCDynamicDPhiParametersHelper::sortDynamicDPhiParametersCollection() {
  std::sort(parameters_.dynamicDPhiParametersCollection_.begin(),
            parameters_.dynamicDPhiParametersCollection_.end(),
            [](const EcalSCDynamicDPhiParameters::DynamicDPhiParameters &p1,
               const EcalSCDynamicDPhiParameters::DynamicDPhiParameters &p2) {
              const auto p1Mins = std::make_pair(p1.eMin, p1.etaMin);
              const auto p2Mins = std::make_pair(p2.eMin, p2.etaMin);
              return p1Mins < p2Mins;
            });
}

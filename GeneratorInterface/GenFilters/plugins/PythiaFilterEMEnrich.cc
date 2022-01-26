#include "GeneratorInterface/GenFilters/plugins/PythiaFilterEMEnrich.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "TLorentzVector.h"

#include <iostream>

using namespace edm;
using namespace std;
using namespace HepMC;

PythiaFilterEMEnrich::PythiaFilterEMEnrich(const edm::ParameterSet& iConfig)
    : token_(consumes<edm::HepMCProduct>(
          edm::InputTag(iConfig.getUntrackedParameter("moduleLabel", std::string("generator")), "unsmeared"))),
      hepMCFilter_(new PythiaHepMCFilterEMEnrich(iConfig)) {}

PythiaFilterEMEnrich::~PythiaFilterEMEnrich() {}

bool PythiaFilterEMEnrich::filter(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  Handle<HepMCProduct> evt;
  iEvent.getByToken(token_, evt);

  const HepMC::GenEvent* myGenEvent = evt->GetEvent();

  return hepMCFilter_->filter(myGenEvent);
}

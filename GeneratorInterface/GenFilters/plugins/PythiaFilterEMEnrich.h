#ifndef PYTHIAFILTEREMENRICH_h
#define PYTHIAFILTEREMENRICH_h

//
// Package:    GeneratorInterface/GenFilters
// Class:      PythiaFilterEMEnrich
//
// Original Author:  Matteo Sani
//
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "GeneratorInterface/Core/interface/PythiaHepMCFilterEMEnrich.h"

namespace edm {
  class HepMCProduct;
}

class PythiaFilterEMEnrich : public edm::global::EDFilter<> {
public:
  explicit PythiaFilterEMEnrich(const edm::ParameterSet&);
  ~PythiaFilterEMEnrich() override;

  bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  const edm::EDGetTokenT<edm::HepMCProduct> token_;

  /** the actual implementation of the filter,
      adapted to be used with HepMCFilterDriver.
      
      We make this a pointer because EDFilter::filter() is const
      while BaseHepMCFilter::filter() is not.
 */
  std::unique_ptr<PythiaHepMCFilterEMEnrich> hepMCFilter_;
};
#endif

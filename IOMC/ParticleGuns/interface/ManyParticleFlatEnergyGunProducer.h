#ifndef IOMC_ParticleGun_ManyParticleFlatEnergyGunProducer_H
#define IOMC_ParticleGun_ManyParticleFlatEnergyGunProducer_H

#include "IOMC/ParticleGuns/interface/BaseFlatGunProducer.h"

#include "IOMC/ParticleGuns/interface/BaseFlatGunProducer.h"

namespace edm
{

  class ManyParticleFlatEnergyGunProducer : public BaseFlatGunProducer
  {

  public:
    ManyParticleFlatEnergyGunProducer(const ParameterSet &);
    ~ManyParticleFlatEnergyGunProducer() override;

  private:

    void produce(Event & e, const EventSetup& es) override;

  protected :
    // data members
    std::vector<double> vEMax;
    std::vector<double> vEMin;
    std::vector<double> vEtaMin;
    std::vector<double> vEtaMax;
    std::vector<double> vPhiMin;
    std::vector<double> vPhiMax;
    std::vector<int> vPartIDs;
  };
}

#endif

#ifndef IOMC_ParticleGun_CloseByParticleFlatEtGunProducer_H
#define IOMC_ParticleGun_CloseByParticleFlatEtGunProducer_H

#include "IOMC/ParticleGuns/interface/BaseFlatGunProducer.h"

namespace edm
{

  class CloseByParticleFlatEtGunProducer : public BaseFlatGunProducer
  {

  public:
    CloseByParticleFlatEtGunProducer(const ParameterSet &);
    ~CloseByParticleFlatEtGunProducer() override;

  private:

    void produce(Event & e, const EventSetup& es) override;

  protected :

    // data members
    double fPtMax,fPtMin,fRMin,fRMax,fZMin,fZMax,fDelta,fPhiMin,fPhiMax;
    int fNParticles;
    bool fPointing = false;
    bool fOverlapping = false;
    bool fRandomShoot = false;
    std::vector<int> fPartIDs;
  };
}

#endif

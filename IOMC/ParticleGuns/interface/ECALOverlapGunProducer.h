#ifndef IOMC_ParticleGun_ECALOverlapGunProducer_H
#define IOMC_ParticleGun_ECALOverlapGunProducer_H

#include "IOMC/ParticleGuns/interface/BaseFlatGunProducer.h"

namespace edm
{

  class ECALOverlapGunProducer : public BaseFlatGunProducer
  {

  public:
    ECALOverlapGunProducer(const ParameterSet &);
    ~ECALOverlapGunProducer() override;

  private:

    void produce(Event & e, const EventSetup& es) override;

  protected :

    // data members
    double fEnFix,fEnMin,fEnMax,fRMin,fRMax,fZMin,fZMax,fDeltaR,fDeltaPhi,fDeltaZ,fPhiMin,fPhiMax;
    int fNParticles;
    bool fPointing = false;
    bool fOverlapping = false;
    bool fRandomShoot = false;
    std::vector<int> fPartIDs;
  };
}

#endif

#include <ostream>

#include "IOMC/ParticleGuns/interface/ManyParticleFlatEnergyGunProducer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Math/interface/Vector3D.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CLHEP/Random/RandFlat.h"

using namespace edm;
using namespace std;

ManyParticleFlatEnergyGunProducer::ManyParticleFlatEnergyGunProducer(const ParameterSet& pset) :
   BaseFlatGunProducer(pset)
{

  ParameterSet defpset ;
  ParameterSet pgun_params =
  pset.getParameter<ParameterSet>("PGunParameters") ;

  vPartIDs = pgun_params.getParameter< vector<int> >("PartID"); 
  vEMin = pgun_params.getParameter< vector<double> >("EMin");
  vEMax = pgun_params.getParameter< vector<double> >("EMax");
  vEtaMin = pgun_params.getParameter< vector<double> >("EtaMin");
  vEtaMax = pgun_params.getParameter< vector<double> >("EtaMax");
  vPhiMin = pgun_params.getParameter< vector<double> >("PhiMin");
  vPhiMax = pgun_params.getParameter< vector<double> >("PhiMax");

  produces<HepMCProduct>("unsmeared");
  produces<GenEventInfoProduct>();

}

ManyParticleFlatEnergyGunProducer::~ManyParticleFlatEnergyGunProducer()
{
   // no need to cleanup GenEvent memory - done in HepMCProduct
}

void ManyParticleFlatEnergyGunProducer::produce(Event &e, const EventSetup& es)
{
   edm::Service<edm::RandomNumberGenerator> rng;
   CLHEP::HepRandomEngine* engine = &rng->getEngine(e.streamID());

   if ( fVerbosity > 0 )
     {
       LogDebug("CloseByParticleFlatEtGunProducer") << " CloseByParticleFlatEtGunProducer : Begin New Event Generation" << endl ;
     }
   fEvt = new HepMC::GenEvent() ;

   // loop over particles
   //
   int barcode = 1 ;

   HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(0.,0.,0.));

   for (unsigned int ip=0; ip<vEMin.size(); ++ip)
   {
     double rnd = CLHEP::RandFlat::shoot(engine,-0.5,(double)(vPartIDs.size()-0.5));
     int partID_index = (int)round(rnd);
     double phi = CLHEP::RandFlat::shoot(engine, vPhiMin.at(ip), vPhiMax.at(ip));
     double eta = CLHEP::RandFlat::shoot(engine, vEtaMin.at(ip), vEtaMax.at(ip));
     double energy = CLHEP::RandFlat::shoot(engine,vEMin.at(ip),vEMax.at(ip));
     const HepPDT::ParticleData *PData = fPDGTable->particle(HepPDT::ParticleID(vPartIDs.at(partID_index))) ;
     double mass   = PData->mass().value() ;
     double mom2 = energy * energy - mass * mass;  
     double mom = 0.;
     if (mom2 > 0.) {
       mom = sqrt(mom2);
     } else {
       mom = 0.;
     }
     double theta = 2. * atan(exp(-eta));
     double px = mom * sin(theta) * cos(phi);
     double py = mom * sin(theta) * sin(phi);
     double pz = mom * cos(theta);
     HepMC::FourVector p(px,py,pz,energy) ;
     HepMC::GenParticle* Part = new HepMC::GenParticle(p,vPartIDs.at(partID_index),1);
     Part->suggest_barcode( barcode ) ;
     barcode++ ;
     Vtx->add_particle_out(Part);

     if (fVerbosity > 0) {
       Vtx->print();
       Part->print();
     }
     fEvt->add_vertex(Vtx);
   }

   fEvt->set_event_number(e.id().event());
   fEvt->set_signal_process_id(20);

   if ( fVerbosity > 0 )
   {
      fEvt->print();
   }

   unique_ptr<HepMCProduct> BProduct(new HepMCProduct());
   BProduct->addHepMCData( fEvt );
   e.put(std::move(BProduct), "unsmeared");

   unique_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(fEvt));
   e.put(std::move(genEventInfo));

   if ( fVerbosity > 0 )
     {
       LogDebug("ManyParticleFlatEnergyGunProducer") << " ManyParticleFlatEnergyGunProducer : Event Generation Done " << endl;
     }
}

#include "FWCore/Framework/interface/HistoryAppender.h"
#include "FWCore/Framework/interface/ResetParameterHelper.h"
#include "DataFormats/Provenance/interface/ProcessConfiguration.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include <string>
#include <cassert>

namespace {
  edm::ProcessHistory initializeEmpty() {
    edm::ProcessHistory tmp{};
    tmp.setProcessHistoryID();
    return tmp;
  }
}  // namespace
static const edm::ProcessHistory s_emptyHistory = initializeEmpty();
static std::map<ParameterSetID, ParameterSetID> rewrittenIDs;
static bool initialized = false;

namespace edm {

  HistoryAppender::HistoryAppender() {}

  std::shared_ptr<ProcessHistory const> HistoryAppender::appendToProcessHistory(
      ProcessHistoryID const& inputPHID, ProcessHistory const* iInputProcessHistory, ProcessConfiguration const& pc) {
      
    // DEBUG
    //std::cout << "DEBUG-appendToProcessHistory: Registry size = " << pset::Registry::instance()->size() << std::endl; 
    // DEBUG 
      
    assert((iInputProcessHistory) == nullptr or (inputPHID == iInputProcessHistory->id()));
    if (m_cachedHistory.get() != nullptr and inputPHID == m_cachedInputPHID) {
      return m_cachedHistory;
    }

    ProcessHistory const* inputProcessHistory = iInputProcessHistory ? iInputProcessHistory : &s_emptyHistory;

    if (inputPHID.isValid()) {
      if (iInputProcessHistory == nullptr) {
        throw Exception(errors::LogicError) << "HistoryAppender::appendToProcessHistory\n"
                                            << "Input ProcessHistory has valid ID but is nullptr\n"
                                            << "Contact a Framework developer\n";
      }
    }

    if (!initialized) {
        auto const& registry = *pset::Registry::instance();
        for (auto const& p : registry) {
             ParameterSet copy(p.second);
             std::vector<std::pair<std::string,std::string>> changes;
             rewriteParameterSet(copy, changes);
             rewrittenIDs[p.first] = copy.registerIt().id();
        }
        initialized = true;
    }

    auto newProcessHistory = std::make_shared<ProcessHistory>();
    //*newProcessHistory = *inputProcessHistory;
    //HACK THE OLD-PROCESS-HISTORY:
    for(auto const& oldPC : *inputProcessHistory){
        auto newID = rewrittenIDs.at(oldPC.parameterSetID()); 
        if(oldPC.processName() == "GEN"){
           //std::cout << "HistoryAppender: " << oldPC.processName() << "  " << oldPC.parameterSetID() << std::endl;
           ProcessConfiguration gen2("inputSTEP1",newID,ReleaseVersion("CMSSW_14_0_19"),oldPC.passID());
           newProcessHistory->push_back(gen2);
        }else if(oldPC.processName() == "SIM"){
           //std::cout << "HistoryAppender: " << oldPC.processName() << "  " << oldPC.parameterSetID() << std::endl;
           ProcessConfiguration sim2("inputSTEP2",newID,ReleaseVersion("CMSSW_14_0_19"),oldPC.passID());
           newProcessHistory->push_back(sim2);
        }else if(oldPC.processName() == "HLT"){
           //std::cout << "HistoryAppender: " << oldPC.processName() << "  " << oldPC.parameterSetID() << std::endl;
           ProcessConfiguration hlt2("inputSTEP3",newID,ReleaseVersion("CMSSW_14_0_19"),oldPC.passID());
           newProcessHistory->push_back(hlt2);
        }else if(oldPC.processName() == "RECO"){
           //std::cout << "HistoryAppender: " << oldPC.processName() << "  " << oldPC.parameterSetID() << std::endl;
           ProcessConfiguration reco2("inputSTEP4",newID,ReleaseVersion("CMSSW_14_0_19"),oldPC.passID());
           newProcessHistory->push_back(reco2);   
        }else if(oldPC.processName() == "PAT"){
           //std::cout << "HistoryAppender: " << oldPC.processName() << "  " << oldPC.parameterSetID() << std::endl;
           ProcessConfiguration pat2("inputSTEP5",newID,ReleaseVersion("CMSSW_14_0_19"),oldPC.passID());
           newProcessHistory->push_back(pat2);              
        }else{
           newProcessHistory->push_back(oldPC);
        }
    }

    checkProcessHistory(*newProcessHistory, pc);
    //newProcessHistory->push_back(pc);
    //HACK THE NEW-PROCESS-HISTORY:
    ParameterSetID newID = rewrittenIDs.at(pc.parameterSetID());
    if(pc.processName() == "GEN"){
       ProcessConfiguration newPC(pc.processName(),newID,pc.releaseVersion(),pc.passID());
       newProcessHistory->push_back(newPC);
    }else if(pc.processName() == "SIM"){
       ProcessConfiguration newPC(pc.processName(),newID,pc.releaseVersion(),pc.passID());
       newProcessHistory->push_back(newPC);
    }else if(pc.processName() == "HLT"){
       ProcessConfiguration newPC(pc.processName(),newID,pc.releaseVersion(),pc.passID());
       newProcessHistory->push_back(newPC);
    }else if(pc.processName() == "RECO"){
       ProcessConfiguration newPC(pc.processName(),newID,pc.releaseVersion(),pc.passID());
       newProcessHistory->push_back(newPC);   
    }else if(pc.processName() == "PAT"){
       ProcessConfiguration newPC(pc.processName(),newID,pc.releaseVersion(),pc.passID());
       newProcessHistory->push_back(newPC);      
    }else{
       newProcessHistory->push_back(pc);
    }
    
    //force it to create the ID
    newProcessHistory->setProcessHistoryID();
    
    // DEBUG
    //std::cout << "\nDEBUG - HistoryAppender returns PHID = " << newProcessHistory->id() << std::endl;
    // DEBUG
    
    m_cachedInputPHID = inputPHID;
    m_cachedHistory = newProcessHistory;
    return m_cachedHistory;
  }

  void HistoryAppender::checkProcessHistory(ProcessHistory const& ph, ProcessConfiguration const& pc) const {
    std::string const& processName = pc.processName();
    for (auto const& process : ph) {
      if (processName == process.processName()) {
        throw edm::Exception(errors::Configuration, "Duplicate Process.")
            << "The process name " << processName << " was already in the ProcessHistory.\n"
            << "Please modify the configuration file to use a distinct process name.\n";
      }
    }
  }
}  // namespace edm

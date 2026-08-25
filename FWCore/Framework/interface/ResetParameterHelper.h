#include "DataFormats/Provenance/interface/BranchType.h"
#include "DataFormats/Provenance/interface/ParameterSetBlob.h"
#include "DataFormats/Provenance/interface/ParameterSetID.h"
#include "DataFormats/Provenance/interface/ProcessHistoryRegistry.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/Registry.h"

#include "TBranch.h"
#include "TTree.h"

#include <cassert>
#include <utility>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <map>

using namespace edm;

void rewriteParameterSet(
      edm::ParameterSet& ps,
      std::vector<std::pair<std::string, std::string>>& changes) {

    static const std::map<std::string, std::string> processMap = {
        {"GEN",  "inputSTEP1"},
        {"SIM",  "inputSTEP2"},
        {"HLT",  "inputSTEP3"},
        {"RECO", "inputSTEP4"},
        {"PAT",  "inputSTEP5"}
    };

    for (auto const& name : ps.getParameterNames()) {

      // tracked InputTag
      if (ps.existsAs<edm::InputTag>(name, true)) {

        auto tag = ps.getParameter<edm::InputTag>(name);
        auto it = processMap.find(tag.process());

        if (it != processMap.end()) {

          edm::InputTag newTag(tag.label(),
                               tag.instance(),
                               it->second);

          changes.emplace_back(tag.encode(), newTag.encode());

          std::cout << "REWRITE InputTag: "
          << tag.encode()
          << " -> "
          << newTag.encode()
          << std::endl;
          ps.addParameter(name, newTag);
        }
      }

      // untracked InputTag
      else if (ps.existsAs<edm::InputTag>(name, false)) {

        auto tag = ps.getUntrackedParameter<edm::InputTag>(name);
        auto it = processMap.find(tag.process());

        if (it != processMap.end()) {

          edm::InputTag newTag(tag.label(),
                               tag.instance(),
                               it->second);

          changes.emplace_back(tag.encode(), newTag.encode());

          ps.addUntrackedParameter(name, newTag);
        }
      }

      // tracked ParameterSet
      else if (ps.existsAs<edm::ParameterSet>(name, true)) {

        auto nested = ps.getParameterSet(name);
        rewriteParameterSet(nested, changes);
        ps.addParameter(name, nested);
      }

      // untracked ParameterSet
      else if (ps.existsAs<edm::ParameterSet>(name, false)) {

        auto nested = ps.getUntrackedParameterSet(name);
        rewriteParameterSet(nested, changes);
        ps.addUntrackedParameter(name, nested);
      }

      // tracked VParameterSet
      else if (ps.existsAs<edm::VParameterSet>(name, true)) {

        auto vec = ps.getParameterSetVector(name);

        for (auto& p : vec)
          rewriteParameterSet(p, changes);

        ps.addParameter(name, vec);
      }

      // untracked VParameterSet
      else if (ps.existsAs<edm::VParameterSet>(name, false)) {

        auto vec = ps.getUntrackedParameterSetVector(name);

        for (auto& p : vec)
          rewriteParameterSet(p, changes);

        ps.addUntrackedParameter(name, vec);
      }
    }
}

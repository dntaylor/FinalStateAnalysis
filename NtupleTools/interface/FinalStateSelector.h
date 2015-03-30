#ifndef PATFINALSTATESELECTOR_3LY6MT3N
#define PATFINALSTATESELECTOR_3LY6MT3N

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

class FinalState;
typedef std::vector<const FinalState*> FinalStatePtrs;
namespace edm {
  class ParameterSet;
}

class FinalStateSelection : public Selector<FinalStatePtrs> {
  public:
    FinalStateSelection(const edm::ParameterSet& pset);
  private:
}



#endif /* end of include guard: PATFINALSTATESELECTOR_3LY6MT3N */

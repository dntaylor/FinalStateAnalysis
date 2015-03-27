#ifndef FINALSTATEPROXY_7FWRI39L
#define FINALSTATEPROXY_7FWRI39L

/*
 *  A Proxy class which wraps a FinalState object.  Essentially wraps a
 *  boost::shared_ptr and allows
 *
 *  Takes ownership of a FinalState passed via the constructor.
 */

#include <boost/shared_ptr.hpp>
class FinalState;

class FinalStateProxy {
  public:
    FinalStateProxy(FinalState* finalState);
    FinalStateProxy();
    const FinalState* get() const;
    const FinalState* operator->() const;
  private:
    boost::shared_ptr<FinalState> finalState_;
};

#endif /* end of include guard: FINALSTATEPROXY_7FWRI39L */

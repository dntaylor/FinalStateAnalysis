/*
 * FinalStateLS
 *
 * A smart container class that holds lumisection level information and has
 * functionality to process that information.
 *
 * */

#ifndef FINALSTATELS_HT41P3HF
#define FINALSTATELS_HT41P3HF

#include <string>
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "FinalStateAnalysis/DataFormats/interface/FinalStateLSFwd.h"

class FinalStateLS {
  public:
    FinalStateLS();
    FinalStateLS(const edm::LuminosityBlockID& id,
        double integratedLuminosity,
        double instantaneousLumi);

    /// Get the ID with lumi & run number
    const edm::LuminosityBlockID& lsID() const;

    /// Get the inst. luminosity in this lumi
    double instantaneousLumi() const;

    /// Get the recorded int. luminosity in this lumi
    double intLumi() const;

    /// Merge this with another product
    bool mergeProduct(FinalStateLS const& other) {
      integratedLumi_ += other.intLumi();
      instaneousLumi_ += other.instantaneousLumi();
      return true;
    }

  private:
    edm::LuminosityBlockID id_;
    double integratedLumi_;
    double instaneousLumi_;
};

#endif /* end of include guard: FINALSTATELS_HT41P3HF */

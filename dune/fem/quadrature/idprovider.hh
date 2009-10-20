#ifndef DUNE_IDPROVIDER_HH
#define DUNE_IDPROVIDER_HH

#include <cstdlib>

namespace Dune {
  //! Singleton that manages a globally unique identifier.
  class IdProvider {
  public:
    //! Access to the singleton object.
    inline static IdProvider& instance() {
      static IdProvider idProvider; 
      return idProvider;
    }

    //! Return a new identifier.
    //! \note Identifiers are never freed.
    size_t newId() {
      return lowestFreeId_++;
    }

  private:
    //! Constructor (for the singleton object)
    IdProvider() :
      lowestFreeId_(0)
    {}

    IdProvider(const IdProvider&);
    IdProvider& operator=(const IdProvider&);

  private:
    size_t lowestFreeId_;
  };

} // end namespace Dune

// #include "idprovider.cc"

#endif

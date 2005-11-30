#ifndef DUNE_IDPROVIDER_HH
#define DUNE_IDPROVIDER_HH

#include <cstdlib>

namespace Dune {
  //! Singleton that manages a globally unique identifier.
  class IdProvider {
  public:
    //! Access to the singleton object.
    static IdProvider& instance() {
      if (!instance_) {
        IdProvider::instance_ = new IdProvider();
      }
      return *IdProvider::instance_;
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

    //! Destructor
    ~IdProvider() { 
      delete instance_; 
    }

    IdProvider(const IdProvider&);
    IdProvider& operator=(const IdProvider&);

  private:
    static IdProvider* instance_;

    size_t lowestFreeId_;
  };

} // end namespace Dune

#endif

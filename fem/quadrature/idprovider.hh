#ifndef DUNE_IDPROVIDER_HH
#define DUNE_IDPROVIDER_HH

namespace Dune {
  class IdProvider {
  public:
    static IdProvider& instance() {
      if (!instance_) {
        IdProvider::instance_ = new IdProvider();
      }
      return *IdProvider::instance_;
    }

    size_t newId() {
      return lowestFreeId_++;
    }

  private:
    IdProvider() :
      lowestFreeId_(0)
    {}

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

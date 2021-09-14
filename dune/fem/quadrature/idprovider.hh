#ifndef DUNE_FEM_IDPROVIDER_HH
#define DUNE_FEM_IDPROVIDER_HH

#include <cstdlib>

#include <dune/fem/storage/singleton.hh>
#include <dune/fem/misc/threads/threadsafevalue.hh>

namespace Dune
{

  namespace Fem
  {

    //! Singleton that manages a globally unique identifier.
    class IdProvider
    {
    public:
      friend class Dune::Fem::Singleton< IdProvider >;

      //! Access to the singleton object.
      static IdProvider& instance()
      {
        return *(Singleton< ThreadSafeValue<IdProvider> >::instance());
      }

      //! Return a new identifier.
      //! \note Identifiers are never freed.
      size_t newId() { return lowestFreeId_++; }

      //! Constructor (for the singleton object)
      IdProvider() :
        lowestFreeId_(0)
      {}

    private:
      IdProvider(const IdProvider&);
      IdProvider& operator=(const IdProvider&);

    private:
      size_t lowestFreeId_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_IDPROVIDER_HH

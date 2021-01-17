#ifndef DUNE_FEM_COMMON_MEMORY_HH
#define DUNE_FEM_COMMON_MEMORY_HH

#include <type_traits>
#include <memory>

#include <dune/common/shared_ptr.hh>

namespace Dune
{

  namespace Fem
  {

    // referenceToSharedPtr
    // --------------------

    template< class T, std::enable_if_t< !std::is_base_of< std::enable_shared_from_this< std::decay_t< T > >, std::decay_t< T > >::value, int > = 0 >
    inline static std::shared_ptr< T > referenceToSharedPtr ( T &t )
    {
      return std::shared_ptr< T >( &t, Dune::null_deleter< T >() );
    }

    template< class T, std::enable_if_t< std::is_base_of< std::enable_shared_from_this< std::decay_t< T > >, std::decay_t< T > >::value, int > = 0 >
    inline static std::shared_ptr< T > referenceToSharedPtr ( T &t )
    try
    {
      // obtain internal weak pointer of enable_shared_from_this class.
      // if the obtained internal weak_ptr of the enable_shared_from_this object is empty
      // then no previously created shared_ptr exists and we create a phony one using the null_deleter
      if( t.weak_from_this().use_count() == 0 )
      {
        return std::shared_ptr< T >( &t, Dune::null_deleter< T >() );
      }
      else
      {
        // return a shared pointer with the same count from the previously
        // created shared_ptr
        return t.shared_from_this();
      }
    }
    // std::bad_weak_ptr is for example thrown when a shared_ptr
    // is created from and empty weak_ptr
    catch( std::bad_weak_ptr& )
    {
      return std::shared_ptr< T >( &t, Dune::null_deleter< T >() );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_MEMORY_HH

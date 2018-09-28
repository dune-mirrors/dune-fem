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
      return t.shared_from_this();
    }
    catch( std::bad_weak_ptr& )
    {
      return std::shared_ptr< T >( &t, Dune::null_deleter< T >() );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_MEMORY_HH

#ifndef DUNE_FEM_COMMON_GET_HH
#define DUNE_FEM_COMMON_GET_HH

#include <utility>

// allows to call get<I>( v ) also if v is not a std::tuple type
template< std::size_t , typename T >
inline T&& get( T&& v )
{
  return std::forward< T >( v );
}

#endif // #ifndef DUNE_FEM_COMMON_GET_HH

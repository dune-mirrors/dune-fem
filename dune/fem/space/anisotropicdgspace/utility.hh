#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_UTILTIY_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_UTILTIY_HH

// dune-common includes
#include <dune/common/static_assert.hh>
#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>

/**
  @file
  @author Christoph Gersbacher
  @brief Provide a simple template meta programm for constructing tuple types
*/


namespace AnisotropicDG
{

  // MakeTuple
  // ---------

  /**
   * \brief Make tuple of given length from type
   *
   * Thus, the following types are equal: 
\code
    Dune::tuple< T, T, T, T >
    MakeTuple< T, 4 >::Type
\endcode
   *
   * \tparam  T    Type
   * \tparam  len  Tuple length 
   */
  template< class T, int len >
  struct MakeTuple;

  template< class T >
  struct MakeTuple< T, 1 >
  {
    typedef Dune::tuple< T > Type;
  };

  template< class T, int len >
  struct MakeTuple
  {
    dune_static_assert( len > 0, "Negative or zero template parameter." );
    typedef typename Dune::JoinTuples< Dune::tuple< T >, typename MakeTuple< T, len-1 >::Type >::type Type;
  };

} //  namespace AnisotropicDG

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_UTILTIY_HH

#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_UTILTIY_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_UTILTIY_HH

// C++ includes
#include <cassert>
#include <cstddef>

// dune-common includes
#include <dune/common/array.hh>
#include <dune/common/densevector.hh>
#include <dune/common/forloop.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>

// dune-fem includes
#include <dune/fem/space/dgspace/legendredgbasefunctions.hh>

/**
  @file
  @author Christoph Gersbacher
  @brief Please doc me
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
   endcode
   *
   * \tparam  T    Type
   * \tparam  len  Tuple length 
   */
  template< class T, int len >
  struct MakeTuple;

  template< class T >
  struct MakeTuple< T, 1 >
  {
    typename Dune::tuple< T > Type;
  };

  template< class T, int len >
  struct MakeTuple
  {
    dune_static_assert( len > 0, "Negative template parameter." );
    typedef Dune::JoinTuples< Dune::tuple< T >, typename MakeTuple< T, len-1 >::Type > Type;
  };



  // NumShapeFunctions
  // -----------------

  /**
   * \brief Provide number of shape functions for AnisotropicDGSpace
   *
   * \tparam  dimension  grid dimension
   *
   * \tparam  maxOrder   maximum polynomal order
   */
  template < int dimension, int maxOrder >
  class NumShapeFunctions
  {
    // this type
    typedef NumShapeFunctions< dimension, maxOrder > ThisType;

    template< int order >
    struct Initialize
    {
      static void apply ( Dune::array< int, maxOrder+1 > &sizes )
      {
        dune_static_assert( order >= 0 && order <= maxOrder, "Invalid template parameter" );
        sizes[ order ] = Dune::Fem::NumLegendreBaseFunctions< order, dimension >::numBaseFct;
      }
    };

  protected:
    //! \brief constructor
    NumShapeFunctions ()
    {
      Dune::ForLoop< Initialize, 0, maxOrder+1 >::apply( sizes_ );
    }

  public:
    //! \brief get singleton
    static ThisType &instance ()
    {
      static ThisType instance_;
      return instance_;
    }

    //! \brief return max number of shape functions
    std::size_t max () const
    {
      std::size_t min = 1;
      for( int i = 0; i < dimension; ++i )
        min *= sizes_[ maxOrder ];
      return min ;
    }

    //! \brief return min number of shape functions
    std::size_t min () const
    {
      std::size_t min = 1;
      for( int i = 0; i < dimension; ++i )
        min *= sizes_[ 0 ];
      return min ;
    }

    //! \brief return number of shape functions for given multi index
    template< class Implementation >
    std::size_t count ( const Dune::DenseVector< Implementation > &multiIndex ) const
    {
      assert( multiIndex.size() == dimension );
      std::size_t count = 1;
      for( int i = 0; i < dimension; ++i )
        count *= sizes_[ multiIndex[ i ] ];
      return count;
    }

  private:
    //forbid copy constructor
    NumShapeFunctions ( const ThisType &other );
    // forbid assignment operator
    ThisType &operator= ( const ThisType &other );

    Dune::array< std::size_t, maxOrder+1 > sizes_;
  };

} //  namespace AnisotropicDG

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_UTILTIY_HH

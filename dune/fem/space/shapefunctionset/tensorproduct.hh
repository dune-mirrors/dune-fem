#ifndef DUNE_FEM_SHAPEFUNCTIONSET_TENSORPRODUCT_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_TENSORPRODUCT_HH

#include <algorithm>
#include <array>
#include <cstddef>
#include <tuple>
#include <utility>

#include <dune/common/fmatrix.hh>
#include <dune/fem/common/forloop.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/tupleutility.hh>
#include <dune/fem/common/coordinate.hh>

namespace Dune
{

  namespace Fem
  {

    // TensorProductShapeFunctionSet
    // -----------------------------

    template< class FunctionSpace, class ShapeFunctionSetTuple >
    class TensorProductShapeFunctionSet
    {
      typedef TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple > ThisType;

      static_assert( (FunctionSpace::dimDomain == std::tuple_size< ShapeFunctionSetTuple >::value),
                     "dimDomain of FunctionSpace must coincide with length of ShapeFunctionSetTuple." );
      static_assert( (FunctionSpace::dimRange == 1),
                     "FunctionSpace must be scalar (i.e., dimRange = 1)." );

      struct Assign;

      template< int i > struct Size;
      template< int i > struct EvaluateAll;
      template< int i > struct JacobianAll;
      template< int i > struct HessianAll;

      static const int dimension = FunctionSpace::dimDomain;

    public:
      typedef FunctionSpace FunctionSpaceType;
      typedef ShapeFunctionSetTuple ShapeFunctionSetTupleType;

      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      TensorProductShapeFunctionSet () = default;

      explicit TensorProductShapeFunctionSet ( const ShapeFunctionSetTupleType &shapeFunctionSetTuple );
      ~TensorProductShapeFunctionSet ();

      TensorProductShapeFunctionSet ( const ThisType &other );
      const ThisType &operator= ( const ThisType &other );

      int order () const;

      std::size_t size () const;

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const;

    private:
      template< class Functor >
      void doEvaluateEach ( int d, RangeType value, std::size_t &index, const RangeFieldType *buffer, Functor functor ) const;
      template< class Functor >
      void doJacobianEach ( int d, JacobianRangeType jacobian, std::size_t &index, const RangeFieldType *buffer, Functor functor ) const;
      template< class Functor >
      void doHessianEach ( int d, HessianRangeType hessian, std::size_t &index, const RangeFieldType *buffer, Functor functor ) const;

      ShapeFunctionSetTuple shapeFunctionSetTuple_;
      std::array< std::size_t, dimension > sizes_;
      RangeFieldType *buffer_;
    };



    // TensorProductShapeFunctionSet::Assign
    // -------------------------------------

    template< class FunctionSpace, class ShapeFunctionSetTuple >
    struct TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >::Assign
    {
      explicit Assign ( RangeFieldType *buffer ) : buffer_( buffer ) {}

      void operator() ( const std::size_t i, const RangeFieldType &value )
      {
        buffer_[ i ] = value;
      }

      template< class T >
      void operator() ( const std::size_t i, const FieldVector< T, 1 > &value )
      {
        (*this)( i, value[ 0 ] );
      }

      template< class T >
      void operator() ( const std::size_t i, const FieldMatrix< T, 1, 1 > &value )
      {
        (*this)( i, value[ 0 ][ 0 ] );
      }

    private:
      RangeFieldType *buffer_;
    };



    // TensorProductShapeFunctionSet::Size
    // -----------------------------------

    template< class FunctionSpace, class ShapeFunctionSetTuple >
    template< int i >
    struct TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >::Size
    {
      static void apply ( const ShapeFunctionSetTuple &tuple, std::array< std::size_t, FunctionSpace::dimDomain > &size )
      {
        size[ i ] = std::get< i >( tuple ).size();
      }
    };



    // TensorProductShapeFunctionSet::EvaluateAll
    // ------------------------------------------

    template< class FunctionSpace, class ShapeFunctionSetTuple >
    template< int i >
    struct TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >::EvaluateAll
    {
      static void apply ( const ShapeFunctionSetTuple &tuple, const DomainType &x, RangeFieldType *&it )
      {
        Dune::FieldVector< DomainFieldType, 1 > xi( x[ i ] );
        std::get< i >( tuple ).evaluateEach( xi, Assign( it ) );
        it += std::get< i >( tuple ).size();
      }
    };



    // TensorProductShapeFunctionSet::JacobianAll
    // ------------------------------------------

    template< class FunctionSpace, class ShapeFunctionSetTuple >
    template< int i >
    struct TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >::JacobianAll
    {
      static void apply ( const ShapeFunctionSetTuple &tuple, const DomainType &x, RangeFieldType *&it )
      {
        Dune::FieldVector< DomainFieldType, 1 > xi( x[ i ] );
        const std::size_t size = std::get< i >( tuple ).size();
        std::get< i >( tuple ).evaluateEach( xi, Assign( it ) );
        std::get< i >( tuple ).jacobianEach( xi, Assign( it+size ) );
        it += 2*size;
      }
    };



    // TensorProductShapeFunctionSet::HessianAll
    // -----------------------------------------

    template< class FunctionSpace, class ShapeFunctionSetTuple >
    template< int i >
    struct TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >::HessianAll
    {
      static void apply ( const ShapeFunctionSetTuple &tuple, const DomainType &x, RangeFieldType *&it )
      {
        Dune::FieldVector< DomainFieldType, 1 > xi( x[ i ] );
        const std::size_t size = std::get< i >( tuple ).size();
        std::get< i >( tuple ).evaluateEach( xi, Assign( it ) );
        std::get< i >( tuple ).jacobianEach( xi, Assign( it+size ) );
        std::get< i >( tuple ).hessianEach( xi, Assign( it+2*size ) );
        it += 3*size;
      }
    };



    // Implementation of TensorProductShapeFunctionSet
    // -----------------------------------------------

    template< class FunctionSpace, class ShapeFunctionSetTuple >
    inline TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >
      ::TensorProductShapeFunctionSet ( const ShapeFunctionSetTupleType &shapeFunctionSetTuple )
    : shapeFunctionSetTuple_( shapeFunctionSetTuple )
    {
      Fem::ForLoop< Size, 0, dimension-1 >::apply( shapeFunctionSetTuple_, sizes_ );
      std::size_t buffer_size = 0;
      for( int i = 0; i < dimension; ++i )
        buffer_size += sizes_[ i ];
      buffer_ = new RangeFieldType[ 3*buffer_size ];
    }


    template< class FunctionSpace, class ShapeFunctionSetTuple >
    inline TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >
      ::~TensorProductShapeFunctionSet ()
    {
      delete[]( buffer_ );
    }


    template< class FunctionSpace, class ShapeFunctionSetTuple >
    inline TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >
      ::TensorProductShapeFunctionSet ( const ThisType &other )
    : shapeFunctionSetTuple_( other.shapeFunctionSetTuple_ )
    {
      std::size_t buffer_size = 0;
      for( int i = 0; i < dimension; ++i )
      {
        sizes_[ i ] = other.sizes_[ i ];
        buffer_size += sizes_[ i ];
      }
      buffer_ = new RangeFieldType[ 3*buffer_size ];
    }


    template< class FunctionSpace, class ShapeFunctionSetTuple >
    inline const typename TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >::ThisType &
    TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >
      ::operator= ( const ThisType &other )
    {
      if( this == &other )
        return *this;
      delete[]( buffer_ );

      shapeFunctionSetTuple_ = other.shapeFunctionSetTuple_;
      std::size_t buffer_size = 0;
      for( int i = 0; i < dimension; ++i )
      {
        sizes_[ i ] = other.sizes_[ i ];
        buffer_size += sizes_[ i ];
      }
      buffer_ = new RangeFieldType[ 3*buffer_size ];
      return *this;
    }


    template< class FunctionSpace, class ShapeFunctionSetTuple >
    inline int TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >::order () const
    {
      int value( 0 );
      Hybrid::forEach( std::make_index_sequence< std::tuple_size< ShapeFunctionSetTupleType >::value >{},
        [ & ]( auto i ){ value = std::max( value, std::get< i >( shapeFunctionSetTuple_ ).order() ); } );
      return value;
    }


    template< class FunctionSpace, class ShapeFunctionSetTuple >
    inline std::size_t
    TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >::size () const
    {
      std::size_t size( 1 );
      for( int i = 0; i < dimension; ++i )
        size *= sizes_[ i ];
      return size;
    }


    template< class FunctionSpace, class ShapeFunctionSetTuple >
    template< class Point, class Functor >
    inline void TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >
      ::evaluateEach ( const Point &x, Functor functor ) const
    {
      RangeFieldType *it = buffer_;
      Fem::ForLoop< EvaluateAll, 0, dimension-1 >::apply( shapeFunctionSetTuple_, coordinate( x ), it );

      std::size_t index = 0;
      doEvaluateEach( 0, RangeType( RangeFieldType( 1 ) ), index, buffer_, functor );
    }


    template< class FunctionSpace, class ShapeFunctionSetTuple >
    template< class Point, class Functor >
    inline void TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >
      ::jacobianEach ( const Point &x, Functor functor ) const
    {
      RangeFieldType *it = buffer_;
      Fem::ForLoop< JacobianAll, 0, dimension-1 >::apply( shapeFunctionSetTuple_, coordinate( x ), it );

      std::size_t index = 0;
      JacobianRangeType jacobian;
      for( int i = 0; i < dimension; ++i )
        jacobian[ 0 ][ i ] = RangeFieldType( 1 );
      doJacobianEach( 0, jacobian, index, buffer_, functor );
    }


    template< class FunctionSpace, class ShapeFunctionSetTuple >
    template< class Point, class Functor >
    inline void TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >
      ::hessianEach ( const Point &x, Functor functor ) const
    {
      RangeFieldType *it = buffer_;
      Fem::ForLoop< HessianAll, 0, dimension-1 >::apply( shapeFunctionSetTuple_, coordinate( x ), it );

      std::size_t index = 0;
      HessianRangeType hessian;
      for( int i = 0; i < dimension; ++i )
        for( int j = 0; j < dimension; ++j )
          hessian[ 0 ][ i ][ j ] = RangeFieldType( 1 );
      doHessianEach( 0, hessian, index, buffer_, functor );
    }


    template< class FunctionSpace, class ShapeFunctionSetTuple >
    template< class Functor >
    inline void TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >
      ::doEvaluateEach ( int d, RangeType value, std::size_t &index, const RangeFieldType *buffer, Functor functor ) const
    {
      if( d < dimension )
      {
        for( std::size_t i = 0; i < sizes_[ d ]; ++i )
        {
          RangeType v( value );
          v[ 0 ] *= buffer[ i ];
          doEvaluateEach( d+1, v, index, buffer+sizes_[ d ], functor );
        }
      }
      else
        functor( index++, value );
    }


    template< class FunctionSpace, class ShapeFunctionSetTuple >
    template< class Functor >
    inline void TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >
      ::doJacobianEach ( int d, JacobianRangeType jacobian, std::size_t &index, const RangeFieldType *buffer, Functor functor ) const
    {
      if( d < dimension )
      {
        for( std::size_t i = 0; i < sizes_[ d ]; ++i )
        {
          JacobianRangeType j( jacobian );
          j[ 0 ][ d ] *= buffer[ i + sizes_[ d ] ];
          for( int k = 1; k < dimension; ++k )
            j[ 0 ][ (d+k)%dimension ] *= buffer[ i ];
          doJacobianEach( d+1, j, index, buffer+2*sizes_[ d ], functor );
        }
      }
      else
        functor( index++, jacobian );
    }


    template< class FunctionSpace, class ShapeFunctionSetTuple >
    template< class Functor >
    inline void TensorProductShapeFunctionSet< FunctionSpace, ShapeFunctionSetTuple >
      ::doHessianEach ( int d, HessianRangeType hessian, std::size_t &index, const RangeFieldType *buffer, Functor functor ) const
    {
      if( d < dimension )
      {
        for( std::size_t i = 0; i < sizes_[ d ]; ++i )
        {
          HessianRangeType h( hessian );
          h[ 0 ][ d ][ d ] *= buffer[ i + 2*sizes_[ d ] ];
          for( int j = 1; j < dimension; ++j )
          {
            h[ 0 ][ (d+j)%dimension ][ d ] *= buffer[ i * sizes_[ d ] ];
            h[ 0 ][ d ][ (d+j)%dimension ] *= buffer[ i * sizes_[ d ] ];
            for( int k = 1; k < dimension; ++k )
              h[ 0 ][ (d+j)%dimension ][ (d+k)%dimension ] *= buffer[ i ];
          }
          doHessianEach( d+1, h, index, buffer+3*sizes_[ d ], functor );
        }
      }
      else
        functor( index++, hessian );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_TENSORPRODUCT_HH

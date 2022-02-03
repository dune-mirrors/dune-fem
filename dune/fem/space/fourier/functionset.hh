#ifndef DUNE_FEM_SPACE_FOURIER_FUNCTIONSET_HH
#define DUNE_FEM_SPACE_FOURIER_FUNCTIONSET_HH

#include <cassert>
#include <cstddef>
#include <limits>
#include <array>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/version.hh>


namespace Dune
{

  namespace Fem
  {

    // FourierFunctionSetSize
    // ----------------------

    template< int dimension, int Order >
    struct FourierFunctionSetSize
    {
      static constexpr int v = Dune::power( (2*Order+1), dimension );
    };



    // FourierFunctionSet
    // ------------------

    template< class FunctionSpace, int Order >
    class FourierFunctionSet;



    // Template specialization for dimDomain = dimRange = 1
    // ----------------------------------------------------

    template< class DomainFieldType, class RangeFieldType, int Order >
    class FourierFunctionSet< FunctionSpace< DomainFieldType, RangeFieldType, 1, 1 >, Order >
    {
      typedef FourierFunctionSet< FunctionSpace< DomainFieldType, RangeFieldType, 1, 1 >, Order > ThisType;

    public:
      typedef FunctionSpace< DomainFieldType, RangeFieldType, 1, 1 > FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef std::size_t SizeType;

      explicit FourierFunctionSet ( int order ) : order_( order ) {}

      int order () const { return order_; }

      static SizeType size () { return FourierFunctionSetSize< 1, Order >::v; }

      template< class Functor >
      static void evaluateEach ( const DomainType &x, Functor functor )
      {
        using std::sin;
        using std::cos;
        functor( 0, RangeType( RangeFieldType( 1 ) / RangeFieldType( 2 ) ) );
        // use recursion:
        // sin((n+1)*x) = sin(n*x)*cos(x) + cos(n*x)*sin(x)
        // cos((n+1)*x) = cos(n*x)*cos(x) - sin(n*x)*sin(x)
        SizeType basisFunction = 1;
        for( int n = 1; n <= Order; ++n )
        {
          functor( basisFunction++, RangeType( cos( n*x[ 0 ] ) ) );
          functor( basisFunction++, RangeType( sin( n*x[ 0 ] ) ) );
        }
      }

      template< class Functor >
      static void jacobianEach ( const DomainType &x, Functor functor )
      {
        using std::sin;
        using std::cos;
        functor( 0, JacobianRangeType( RangeFieldType( 0 ) ) );
        SizeType basisFunction = 1;
        for( int n = 1; n <= Order; ++n )
        {
          functor( basisFunction++, JacobianRangeType( -n*sin( n*x[ 0 ] ) ) );
          functor( basisFunction++, JacobianRangeType( n*cos( n*x[ 0 ] ) ) );
        }
      }

      template< class Functor >
      static void hessianEach ( const DomainType &x, Functor functor )
      {
        using std::sin;
        using std::cos;
        functor( 0, HessianRangeType( RangeFieldType( 0 ) ) );
        SizeType basisFunction = 1;
        for( int n = 1; n <= Order; ++n )
        {
          functor( basisFunction++, HessianRangeType( -(n*n)*cos( n*x[ 0 ] ) ) );
          functor( basisFunction++, HessianRangeType( -(n*n)*sin( n*x[ 0 ] ) ) );
        }
      }
    private:
      int order_;
    };



    // Template specialization for dimDomain > 1, dimRange = 1
    // -------------------------------------------------------

    template< class DomainFieldType, class RangeFieldType, int dimDomain, int Order >
    class FourierFunctionSet< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, Order >
    {
      typedef FourierFunctionSet< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, Order > ThisType;

    public:
      typedef FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 > FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef std::size_t SizeType;

    private:
      // number of Fourier basis function for dimDomain = 1
      static const int buffer_size = FourierFunctionSetSize< 1, Order >::v;

      // tags used for building cache
      struct Evaluate {};  //< evaluate basis functions
      struct Jacobian {};  //< evaluate basis functions and jacobians
      struct Hessian {};   //< evaluate basis functions, jacobians, and hessians

      struct Assign;

    protected:
      // multi index used for tensor product basis functions
      typedef Dune::FieldVector< int, dimDomain > MultiIndexType;

      // iterator type and methods for accessing multi indices
      struct MultiIndexIterator;
      typedef MultiIndexIterator IteratorType;
      static IteratorType begin () { return IteratorType::begin(); }
      static IteratorType end () { return IteratorType::end(); }

    public:
      explicit FourierFunctionSet ( int order ) : order_( order ) {}

      int order () const { return order_; }

      static SizeType size () { return FourierFunctionSetSize< dimDomain, Order >::v; }

      template< class Functor >
      void evaluateEach ( const DomainType &x, Functor functor ) const
      {
        prepare( Evaluate(), x );
        SizeType index( 0 );
        const IteratorType end = ThisType::end();
        for( IteratorType it = ThisType::begin(); it != end; ++it, ++index )
        {
          RangeType value;
          evaluate( *it, value );
          assert( index == IteratorType::index( *it ) );
          functor( index, value );
        }
        assert( index == size() );
      }

      template< class Functor >
      void jacobianEach ( const DomainType &x, Functor functor ) const
      {
        prepare( Jacobian(), x );
        SizeType index( 0 );
        const IteratorType end = ThisType::end();
        for( IteratorType it = ThisType::begin(); it != end; ++it, ++index )
        {
          JacobianRangeType jacobian;
          evaluate( *it, jacobian );
          assert( index == IteratorType::index( *it ) );
          functor( index, jacobian );
        }
        assert( index == size() );
      }

      template< class Functor >
      void hessianEach ( const DomainType &x, Functor functor ) const
      {
        prepare( Hessian(), x );
        SizeType index( 0 );
        const IteratorType end = ThisType::end();
        for( IteratorType it = ThisType::begin(); it != end; ++it, ++index )
        {
          HessianRangeType hessian;
          evaluate( *it, hessian );
          assert( index == IteratorType::index( *it ) );
          functor( index, hessian );
        }
        assert( index == size() );
      }

    protected:
      // evaluate tensor product basis function
      void evaluate ( const MultiIndexType &multiIndex, RangeType &value ) const
      {
        value = RangeType( RangeFieldType( 1 ) );
        for( SizeType i = 0; i < dimDomain; ++i )
          value *= buffer_[ i ][ multiIndex[ i ] ];
      }

      // evaluate jacobian of tensor product basis function
      void evaluate ( const MultiIndexType &multiIndex, JacobianRangeType &jacobian ) const
      {
        jacobian = JacobianRangeType( 1 );
        for( int k = 0; k < dimDomain; ++k )
        {
          const RangeFieldType phi = buffer_[ k ][ multiIndex[ k ] ];
          const RangeFieldType dphi = buffer_[ k ][ buffer_size + multiIndex[ k ] ];
          for( int i = 0; i < dimDomain; ++i )
            jacobian[ 0 ][ i ] *= (k == i ? dphi : phi);
        }
      }

      // evaluate hessian of tensor product basis function
      void evaluate ( const MultiIndexType &multiIndex, HessianRangeType &hessian ) const
      {
        for( int i = 0; i < dimDomain; ++i )
          for( int j = 0; j < dimDomain; ++j )
            hessian[ 0 ][ i ][ j ] = RangeFieldType( 1 );

        for( int k = 0; k < dimDomain; ++k )
        {
          const RangeFieldType phi = buffer_[ k ][ multiIndex[ k ] ];
          const RangeFieldType dphi = buffer_[ k ][ buffer_size + multiIndex[ k ] ];
          for( int i = 0; i < dimDomain; ++i )
          {
            hessian[ 0 ][ i ][ i ] *= (k == i ? buffer_[ i ][ 2*buffer_size + multiIndex[ i ] ] : phi);
            for( int j = i+1; j < dimDomain; ++j )
            {
              RangeFieldType tmp = ( k == i || k == j ) ? dphi : phi;
              hessian[ 0 ][ i ][ j ] *= tmp;
              hessian[ 0 ][ j ][ i ] *= tmp;
            }
          }
        }
      }

      // methods for building cache
      void prepare ( const Evaluate &, const DomainType &x ) const;
      void prepare ( const Jacobian &, const DomainType &x ) const;
      void prepare ( const Hessian &, const DomainType &x ) const;

    private:
      int order_;
      // cache for evaluation of basis functions
      mutable std::array< std::array< RangeFieldType, 3*buffer_size>, dimDomain > buffer_;
    };



    // Implementation of FourierFunctionSet::Assign
    // --------------------------------------------

    template< class DomainFieldType, class RangeFieldType, int dimDomain, int Order >
    struct FourierFunctionSet< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, Order >::Assign
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



    // Implementation of FourierFunctionSet::MultiIndexIterator
    // --------------------------------------------------------

    template< class DomainFieldType, class RangeFieldType, int dimDomain, int Order >
    struct FourierFunctionSet< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, Order >::MultiIndexIterator
    {
      typedef MultiIndexIterator ThisType;

    protected:
      typedef int IndexType;

      explicit MultiIndexIterator ( IndexType n ) : multiIndex_( n ) {}

      static IndexType invalidIndex () { return std::numeric_limits< IndexType >::max(); }

      static const int N = FourierFunctionSetSize< 1, Order >::v;

    public:
      static ThisType begin () { return ThisType( 0 ); }
      static ThisType end () { return ThisType( invalidIndex() ); }

      ThisType operator++ ()
      {
        // try to increment and leave eventually
        for( int i = 0; i < dimDomain; ++i )
        {
          const int j = dimDomain-i-1;
          if( ++multiIndex_[ j ] < N )
            return *this;
          multiIndex_[ j ] = 0;
        }

        // otherwise, reset this iterator to end iterator
        *this = end();
        return *this;
      }

      bool operator== ( const ThisType &other ) const { return ( multiIndex_ == other.multiIndex_ ); }

      bool operator!= ( const ThisType &other ) const { return !( *this == other ); }

      const MultiIndexType &operator* () const { return multiIndex_; }

      const MultiIndexType *operator-> () const { return &multiIndex_; }

      SizeType index () const { return index( multiIndex_ ); }

      static SizeType index ( const MultiIndexType &multiIndex )
      {
        SizeType index = 0, factor = 1;
        for( int i = dimDomain-1; i >= 0; --i )
        {
          index += multiIndex[ i ]*factor;
          factor *= N;
        }
        assert( index < size() );
        return index;
      }

    private:
      MultiIndexType multiIndex_;
    };



    // Implementation of FourierFunctionSet
    // ------------------------------------

    template< class DomainFieldType, class RangeFieldType, int dimDomain, int Order >
    void FourierFunctionSet< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, Order >
      ::prepare ( const Evaluate &, const DomainType &x ) const
    {
      typedef FourierFunctionSet< FunctionSpace< DomainFieldType, RangeFieldType, 1, 1 >, Order > FunctionSetImp;
      for( SizeType i = 0; i < dimDomain; ++i )
      {
        RangeFieldType *it = &buffer_[ i ][ 0 ];
        Dune::FieldVector< DomainFieldType, 1 > y( x[ i ] );
        FunctionSetImp::evaluateEach( y, Assign( it ) );
      }
    };


    template< class DomainFieldType, class RangeFieldType, int dimDomain, int Order >
    void FourierFunctionSet< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, Order >
      ::prepare ( const Jacobian &, const DomainType &x ) const
    {
      typedef FourierFunctionSet< FunctionSpace< DomainFieldType, RangeFieldType, 1, 1 >, Order > FunctionSetImp;
      for( SizeType i = 0; i < dimDomain; ++i )
      {
        RangeFieldType *it = &buffer_[ i ][ 0 ];
        Dune::FieldVector< DomainFieldType, 1 > y( x[ i ] );
        FunctionSetImp::evaluateEach( y, Assign( it ) );
        FunctionSetImp::jacobianEach( y, Assign( it+buffer_size ) );
      }
    };


    template< class DomainFieldType, class RangeFieldType, int dimDomain, int Order >
    void FourierFunctionSet< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, Order >
      ::prepare ( const Hessian &, const DomainType &x ) const
    {
      typedef FourierFunctionSet< FunctionSpace< DomainFieldType, RangeFieldType, 1, 1 >, Order > FunctionSetImp;
      for( SizeType i = 0; i < dimDomain; ++i )
      {
        RangeFieldType *it = &buffer_[ i ][ 0 ];
        Dune::FieldVector< DomainFieldType, 1 > y( x[ i ] );
        FunctionSetImp::evaluateEach( y, Assign( it ) );
        FunctionSetImp::jacobianEach( y, Assign( it+buffer_size) );
        FunctionSetImp::hessianEach( y, Assign( it+2*buffer_size ) );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FOURIER_FUNCTIONSET_HH

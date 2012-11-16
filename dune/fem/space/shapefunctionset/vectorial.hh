#ifndef DUNE_FEM_SHAPEFUNCTIONSET_VECTORIAL_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_VECTORIAL_HH

// C++ includes
#include <algorithm>
#include <cstddef>

// dune-fem includes
#include <dune/fem/space/common/functionspace.hh>


namespace Dune
{

  namespace Fem 
  {

    // MakeVectorialExpression
    // -----------------------

    template< class Scalar, class Vectorial >
    class MakeVectorialExpression;

    template< class K, int dimR >
    class MakeVectorialExpression< FieldVector< K, 1 >, FieldVector< K, dimR > >
    {
      typedef MakeVectorialExpression< FieldVector< K, 1 >, FieldVector< K, dimR > > ThisType;
      
    public:
      typedef FieldVector< K, 1 > ScalarType;
      typedef FieldVector< K, dimR > VectorialType;

      typedef typename VectorialType::size_type size_type;

      MakeVectorialExpression ( int k, const ScalarType &scalar )
      : k_( k ),
        scalar_( scalar )
      {}

      operator VectorialType () const
      {
        VectorialType vectorial( K( 0 ) );
        vectorial[ k_ ] = scalar_[ 0 ];
        return vectorial;
      }

      const ThisType &operator*= ( const K &s )
      {
        scalar_ *= s;
        return *this;
      }

      const ThisType &operator/= ( const K &s )
      {
        scalar_ /= s;
        return *this;
      }

      bool operator== ( const ThisType &other ) const
      {
        return ((k_ == other.k_) && (scalar_ == other.scalar_));
      }

      bool operator!= ( const ThisType &other ) const
      {
        return ((k_ != other.k_) || (scalar_ != other.scalar_));
      }

      K operator* ( const ThisType &other ) const
      {
        return (k_ == other.k_ ? scalar_ * other.scalar_ : K( 0 ));
      }

      K operator* ( const VectorialType &other ) const
      {
        return (scalar_[ 0 ] * other[ k_ ]);
      }

      K one_norm () const { return scalar_.one_norm(); }
      K two_norm () const { return scalar_.two_norm(); }
      K two_norm2 () const { return scalar_.two_norm2(); }
      K infinity_norm () const { return scalar_.infinity_norm(); }

      size_type size () const { return dimR; }

      friend K operator* ( const VectorialType &a, ThisType &b ) { return b*a; }

      friend void axpy ( const K &a, const ThisType &x, VectorialType &y )
      {
        y[ x.k_ ] += a * x.scalar_[ 0 ];
      }

    protected:
      int k_;
      ScalarType scalar_;
    };

    template< class K, int dimR, int dimD >
    class MakeVectorialExpression< FieldMatrix< K, 1, dimD >, FieldMatrix< K, dimR, dimD > >
    {
      typedef MakeVectorialExpression< FieldMatrix< K, 1, dimD >, FieldMatrix< K, dimR, dimD > > ThisType;
      
    public:
      typedef FieldMatrix< K, 1, dimD > ScalarType;
      typedef FieldMatrix< K, dimR, dimD > VectorialType;

      typedef typename VectorialType::size_type size_type;

      MakeVectorialExpression ( int k, const ScalarType &scalar )
      : k_( k ),
        scalar_( scalar )
      {}

      operator VectorialType () const
      {
        VectorialType vectorial( K( 0 ) );
        vectorial[ k_ ] = scalar_[ 0 ];
        return vectorial;
      }

      const ThisType &operator*= ( const K &s )
      {
        scalar_ *= s;
        return *this;
      }

      const ThisType &operator/= ( const K &s )
      {
        scalar_ /= s;
        return *this;
      }

      bool operator== ( const ThisType &other ) const
      {
        return ((k_ == other.k_) && (scalar_ == other.scalar_));
      }

      bool operator!= ( const ThisType &other ) const
      {
        return ((k_ != other.k_) || (scalar_ != other.scalar_));
      }

      template< class X, class Y >
      void mv ( const X &x, Y &y ) const
      {
        for( size_type i = 0; i < rows(); ++i )
          y[ i ] = K( 0 );
        for( size_type j= 0; j < cols(); ++j )
          y[ k_ ] += scalar_[ k_ ][ j ] * x[ j ];
      }

      template< class X, class Y >
      void mtv ( const X &x, Y &y ) const
      {
        for( size_type i = 0; i < rows(); ++i )
          y[ i ] = scalar_[ i ][ k_ ] * x[ k_ ];
      }

      template< class X, class Y >
      void umv ( const X &x, Y &y ) const
      {
        for( size_type j= 0; j < cols(); ++j )
          y[ k_ ] += scalar_[ k_ ][ j ] * x[ j ];
      }

      template< class X, class Y >
      void umtv ( const X &x, Y &y ) const
      {
        for( size_type i = 0; i < rows(); ++i )
          y[ i ] += scalar_[ i ][ k_ ] * x[ k_ ];
      }

      template< class X, class Y >
      void mmv ( const X &x, Y &y ) const
      {
        for( size_type j= 0; j < cols(); ++j )
          y[ k_ ] -= scalar_[ k_ ][ j ] * x[ j ];
      }

      template< class X, class Y >
      void mmtv ( const X &x, Y &y ) const
      {
        for( size_type i = 0; i < rows(); ++i )
          y[ i ] -= scalar_[ i ][ k_ ] * x[ k_ ];
      }

      K frobenius_norm () const { return scalar_.frobenius_norm(); }
      K frobenius_norm2 () const { return scalar_.frobenius_norm2(); }
      K infinity_norm () const { return scalar_.infinity_norm(); }

      K determinant () const { return (dimR == 1 ? scalar_.determinant() : K( 0 )); }

      size_type N () const { return rows(); }
      size_type M () const { return cols(); }

      size_type rows () const { return dimR; }
      size_type cols () const { return scalar_.cols(); }

      friend void axpy ( const K &a, const ThisType &x, VectorialType &y )
      {
        y[ x.k_ ].axpy( a, x.scalar_[ 0 ] );
      }

    protected:
      int k_;
      ScalarType scalar_;
    };



    // Auxilliary Functions for MakeVectorialExpression
    // ------------------------------------------------

    template< class Scalar, class Vectorial >
    inline bool
    operator== ( const Vectorial &a,
                 const MakeVectorialExpression< Scalar, Vectorial > &b )
    {
      return (a == static_cast< Vectorial >( b ));
    }

    template< class Scalar, class Vectorial >
    inline bool
    operator!= ( const Vectorial &a,
                 const MakeVectorialExpression< Scalar, Vectorial > &b )
    {
      return (a != static_cast< Vectorial >( b ));
    }

    template< class Scalar, class Vectorial >
    inline bool
    operator== ( const MakeVectorialExpression< Scalar, Vectorial > &a,
                 const Vectorial &b )
    {
      return (static_cast< Vectorial >( a ) == b);
    }

    template< class Scalar, class Vectorial >
    inline bool
    operator!= ( const MakeVectorialExpression< Scalar, Vectorial > &a,
                 const Vectorial &b )
    {
      return (static_cast< Vectorial >( a ) != b);
    }



    // MakeVectorial
    // -------------

    template< class ScalarFunctionSpace, class RangeVector >
    struct MakeVectorial;

    template< class DomainField, class RangeField, int dimD, int dimR >
    struct MakeVectorial< FunctionSpace< DomainField, RangeField, dimD, 1 >, FieldVector< RangeField, dimR > >
    {
      typedef FunctionSpace< DomainField, RangeField, dimD, 1 > ScalarFunctionSpaceType;
      typedef FunctionSpace< DomainField, RangeField, dimD, dimR > VectorialFunctionSpaceType;

      static const int dimRangeFactor = dimR;

      static typename VectorialFunctionSpaceType::RangeType
      makeVectorial ( int k, const typename ScalarFunctionSpaceType::RangeType &scalarValue )
      {
        typename VectorialFunctionSpaceType::RangeType vectorialValue( RangeField( 0 ) );
        vectorialValue[ k ] = scalarValue[ 0 ];
        return vectorialValue;
      }

      static typename VectorialFunctionSpaceType::JacobianRangeType
      makeVectorial ( int k, const typename ScalarFunctionSpaceType::JacobianRangeType &scalarValue )
      {
        typename VectorialFunctionSpaceType::JacobianRangeType vectorialValue( RangeField( 0 ) );
        vectorialValue[ k ] = scalarValue[ 0 ];
        return vectorialValue;
      }

      static typename VectorialFunctionSpaceType::HessianRangeType
      makeVectorial ( int k, const typename ScalarFunctionSpaceType::HessianRangeType &scalarValue )
      {
        typename VectorialFunctionSpaceType::HessianRangeType 
          vectorialValue( typename VectorialFunctionSpaceType::HessianRangeType::value_type( RangeField( 0 ) ) );
        vectorialValue[ k ] = scalarValue[ 0 ];
        return vectorialValue;
      }
    };



    // VectorialShapeFunctionSet
    // -------------------------

    template< class ScalarShapeFunctionSet, class RangeVector >
    class VectorialShapeFunctionSet
    {
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector > ThisType;

    public:
      typedef typename MakeVectorial< typename ScalarShapeFunctionSet::FunctionSpaceType, RangeVector >::VectorialFunctionSpaceType FunctionSpaceType;
      typedef ScalarShapeFunctionSet ScalarShapeFunctionSetType;

    protected:
      typedef typename ScalarShapeFunctionSetType::FunctionSpaceType ScalarFunctionSpaceType;
      typedef MakeVectorial< ScalarFunctionSpaceType, RangeVector > MakeVectorialType;

      static const int dimRangeFactor = MakeVectorialType::dimRangeFactor;

      template< class Functor >
      struct VectorialFunctor;

    public:
      explicit VectorialShapeFunctionSet ( const ScalarShapeFunctionSetType &scalarShapeFunctionSet = ScalarShapeFunctionSetType() )
      : scalarShapeFunctionSet_( scalarShapeFunctionSet )
      {}

      const ScalarShapeFunctionSetType &scalarShapeFunctionSet () const { return scalarShapeFunctionSet_; }

      // Shape Function Set Interface Methods
      std::size_t size () const { return dimRangeFactor * scalarShapeFunctionSet().size(); }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor > 
      void jacobianEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor > 
      void hessianEach ( const Point &x, Functor functor ) const;
     
    protected:
      ScalarShapeFunctionSet scalarShapeFunctionSet_;
    };



    // VectorialShapeFunctionSet::VectorialFunctor
    // -------------------------------------------

    template< class ScalarShapeFunctionSet, class RangeVector >
    template< class Functor >
    struct VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector >::VectorialFunctor
    {
      explicit VectorialFunctor ( const Functor &functor )
      : functor_( functor )
      {}

      template< class Value >
      void operator() ( const std::size_t i, const Value &value )
      {
        for( int k = 0; k < dimRangeFactor; ++k )
          functor_( i*dimRangeFactor+k, MakeVectorialType::makeVectorial( k, value ) );
      }

    private:
      Functor functor_;
    };



    // Implementation of VectorialShapeFunctionSet
    // -------------------------------------------

    template< class ScalarShapeFunctionSet, class RangeVector >
    template< class Point, class Functor >
    inline void VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector >
      ::evaluateEach ( const Point &x, Functor functor ) const
    {
      scalarShapeFunctionSet().evaluateEach( x, VectorialFunctor< Functor >( functor ) );
    }
    

    template< class ScalarShapeFunctionSet, class RangeVector >
    template< class Point, class Functor >
    inline void VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector >
      ::jacobianEach ( const Point &x, Functor functor ) const
    {
      scalarShapeFunctionSet().jacobianEach( x, VectorialFunctor< Functor >( functor ) );
    }


    template< class ScalarShapeFunctionSet, class RangeVector >
    template< class Point, class Functor >
    inline void VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector >
      ::hessianEach ( const Point &x, Functor functor ) const
    {
      scalarShapeFunctionSet().hessianEach( x, VectorialFunctor< Functor >( functor ) );
    }

  } // namespace Fem 

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_VECTORIAL_HH

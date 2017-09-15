#ifndef DUNE_FEM_SHAPEFUNCTIONSET_VECTORIAL_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_VECTORIAL_HH

// C++ includes
#include <cstddef>

#include <algorithm>
#include <type_traits>

// dune-fem includes
#include <dune/fem/common/fmatrixcol.hh>
#include <dune/fem/space/basisfunctionset/functor.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/storage/subvector.hh>

namespace Dune
{

  namespace Fem
  {

    // MakeVectorialTraits
    // -------------------

    template< class Scalar, class Vectorial, class = void >
    struct MakeVectorialTraits;

    template< class K, int dimR, int dimNR >
    struct MakeVectorialTraits< FieldVector< K, dimR >, FieldVector< K, dimNR >, std::enable_if_t< dimNR % dimR == 0 > >
    {
      typedef FieldVector< K, dimR > ScalarType;
      typedef FieldVector< K, dimNR > VectorialType;

      typedef typename FieldTraits< VectorialType >::field_type field_type;
      typedef typename VectorialType::size_type IndexType;
      typedef typename VectorialType::size_type size_type;

      typedef SubVector< const VectorialType, StaticOffsetSubMapper< dimR > > ConstComponentType;
      typedef SubVector< VectorialType, StaticOffsetSubMapper< dimR > > ComponentType;

      static const size_type factor = dimNR / dimR;

      static IndexType begin () { return IndexType( 0 ); }
      static IndexType end () { return IndexType( factor ); }

      static VectorialType zeroVectorial () { return VectorialType( K( field_type( 0 ) ) ); }

      static ConstComponentType component ( const VectorialType &x, const IndexType &i )
      {
        return ConstComponentType( x, StaticOffsetSubMapper< dimR >( i*dimR ) );
      }

      static ComponentType component ( VectorialType &x, const IndexType &i )
      {
        return ComponentType( x, StaticOffsetSubMapper< dimR >( i*dimR ) );
      }

      static size_type index ( const IndexType &i ) { return i; }
    };

    template< class K, int dimR, int dimD, int dimNR >
    struct MakeVectorialTraits< FieldMatrix< K, dimR, dimD >, FieldMatrix< K, dimNR, dimD >, std::enable_if_t< dimNR % dimR == 0 > >
    {
      typedef FieldMatrix< K, dimR, dimD > ScalarType;
      typedef FieldMatrix< K, dimNR, dimD > VectorialType;

      typedef typename FieldTraits< VectorialType >::field_type field_type;
      typedef typename VectorialType::size_type IndexType;
      typedef typename VectorialType::size_type size_type;

      typedef SubRowMatrix< const VectorialType, StaticOffsetSubMapper< dimR > > ConstComponentType;
      typedef SubRowMatrix< VectorialType, StaticOffsetSubMapper< dimR > > ComponentType;

      static const size_type factor = dimNR / dimR;

      static IndexType begin () { return IndexType( 0 ); }
      static IndexType end () { return IndexType( factor ); }

      static VectorialType zeroVectorial () { return VectorialType( K( field_type( 0 ) ) ); }

      static ConstComponentType component ( const VectorialType &x, const IndexType &i )
      {
        return ConstComponentType( x, StaticOffsetSubMapper< dimR >( i*dimR ) );
      }

      static ComponentType component ( VectorialType &x, const IndexType &i )
      {
        return ComponentType( x, StaticOffsetSubMapper< dimR >( i*dimR ) );
      }

      static size_type index ( const IndexType &i ) { return i; }
    };



    // MakeVectorialExpression
    // -----------------------

    template< class Scalar, class Vectorial >
    class BasicMakeVectorialExpression
    {
      typedef BasicMakeVectorialExpression< Scalar, Vectorial > ThisType;

      typedef MakeVectorialTraits< Scalar, Vectorial > Traits;

    public:
      typedef typename Traits::ScalarType ScalarType;
      typedef typename Traits::VectorialType VectorialType;

      typedef typename Traits::field_type field_type;
      typedef typename Traits::IndexType IndexType;
      typedef typename Traits::size_type size_type;

      BasicMakeVectorialExpression ( const IndexType &component, const ScalarType &scalar )
      : component_( component ),
        scalar_( scalar )
      {}

      operator VectorialType () const
      {
        VectorialType vectorial = Traits::zeroVectorial();
        Traits::component( vectorial, component() ) = scalar();
        return vectorial;
      }

      const ThisType &operator*= ( const field_type &s )
      {
        scalar() *= s;
        return *this;
      }

      const ThisType &operator/= ( const field_type &s )
      {
        scalar() /= s;
        return *this;
      }

      const IndexType &component () const { return component_; }

      const ScalarType &scalar () const { return scalar_; }
      ScalarType &scalar () { return scalar_; }

    protected:
      IndexType component_;
      ScalarType scalar_;
    };



    // MakeVectorialExpression
    // -----------------------

    template< class Scalar, class Vectorial >
    class MakeVectorialExpression
    : public BasicMakeVectorialExpression< Scalar, Vectorial >
    {
      typedef MakeVectorialExpression< Scalar, Vectorial > ThisType;
      typedef BasicMakeVectorialExpression< Scalar, Vectorial > BaseType;

    public:
      typedef typename BaseType::IndexType IndexType;
      typedef typename BaseType::ScalarType ScalarType;

      MakeVectorialExpression ( const IndexType &component, const ScalarType &scalar )
      : BaseType( component, scalar )
      {}
    };

    template< class K, int dimR >
    class MakeVectorialExpression< FieldVector< K, 1 >, FieldVector< K, dimR > >
    : public BasicMakeVectorialExpression< FieldVector< K, 1 >, FieldVector< K, dimR > >
    {
      typedef MakeVectorialExpression< FieldVector< K, 1 >, FieldVector< K, dimR > > ThisType;
      typedef BasicMakeVectorialExpression< FieldVector< K, 1 >, FieldVector< K, dimR > > BaseType;

    public:
      typedef typename BaseType::ScalarType ScalarType;
      typedef typename BaseType::VectorialType VectorialType;

      typedef typename BaseType::field_type field_type;
      typedef typename BaseType::IndexType IndexType;
      typedef typename BaseType::size_type size_type;

      using BaseType::component;
      using BaseType::scalar;

      MakeVectorialExpression ( const IndexType &component, const ScalarType &scalar )
      : BaseType( component, scalar )
      {}

      field_type operator* ( const ThisType &other ) const
      {
        return (component() == other.component() ? scalar() * other.scalar() : field_type( 0 ));
      }

      field_type operator* ( const VectorialType &other ) const
      {
        return (scalar()[ 0 ] * other[ component() ]);
      }

      field_type one_norm () const { return scalar().one_norm(); }
      field_type two_norm () const { return scalar().two_norm(); }
      field_type two_norm2 () const { return scalar().two_norm2(); }
      field_type infinity_norm () const { return scalar().infinity_norm(); }

      size_type size () const { return dimR; }

      friend field_type operator* ( const VectorialType &a, ThisType &b ) { return b*a; }
    };

    template< class K, int dimR, int dimD >
    class MakeVectorialExpression< FieldMatrix< K, 1, dimD >, FieldMatrix< K, dimR, dimD > >
    : public BasicMakeVectorialExpression< FieldMatrix< K, 1, dimD >, FieldMatrix< K, dimR, dimD > >
    {
      typedef MakeVectorialExpression< FieldMatrix< K, 1, dimD >, FieldMatrix< K, dimR, dimD > > ThisType;
      typedef BasicMakeVectorialExpression< FieldMatrix< K, 1, dimD >, FieldMatrix< K, dimR, dimD > > BaseType;

    public:
      typedef typename BaseType::ScalarType ScalarType;
      typedef typename BaseType::VectorialType VectorialType;

      typedef typename BaseType::field_type field_type;
      typedef typename BaseType::IndexType IndexType;
      typedef typename BaseType::size_type size_type;

      using BaseType::component;
      using BaseType::scalar;

      MakeVectorialExpression ( const IndexType &component, const ScalarType &scalar )
      : BaseType( component, scalar )
      {}

      template< class X, class Y >
      void mv ( const X &x, Y &y ) const
      {
        for( size_type i = 0; i < rows(); ++i )
          y[ i ] = field_type( 0 );
        for( size_type j= 0; j < cols(); ++j )
          y[ component() ] += scalar()[ component() ][ j ] * x[ j ];
      }

      template< class X, class Y >
      void mtv ( const X &x, Y &y ) const
      {
        for( size_type i = 0; i < rows(); ++i )
          y[ i ] = scalar()[ i ][ component() ] * x[ component() ];
      }

      template< class X, class Y >
      void umv ( const X &x, Y &y ) const
      {
        for( size_type j= 0; j < cols(); ++j )
          y[ component() ] += scalar()[ component() ][ j ] * x[ j ];
      }

      template< class X, class Y >
      void umtv ( const X &x, Y &y ) const
      {
        for( size_type i = 0; i < rows(); ++i )
          y[ i ] += scalar()[ i ][ component() ] * x[ component() ];
      }

      template< class X, class Y >
      void mmv ( const X &x, Y &y ) const
      {
        for( size_type j= 0; j < cols(); ++j )
          y[ component() ] -= scalar()[ component() ][ j ] * x[ j ];
      }

      template< class X, class Y >
      void mmtv ( const X &x, Y &y ) const
      {
        for( size_type i = 0; i < rows(); ++i )
          y[ i ] -= scalar()[ i ][ component() ] * x[ component() ];
      }

      field_type frobenius_norm () const { return scalar().frobenius_norm(); }
      field_type frobenius_norm2 () const { return scalar().frobenius_norm2(); }
      field_type infinity_norm () const { return scalar().infinity_norm(); }

      field_type determinant () const { return (dimR == 1 ? scalar().determinant() : field_type( 0 )); }

      size_type N () const { return rows(); }
      size_type M () const { return cols(); }

      size_type rows () const { return dimR; }
      size_type cols () const { return scalar().cols(); }
    };



    // Auxilliary Functions for MakeVectorialExpression
    // ------------------------------------------------

    template< class Scalar, class Vectorial >
    inline bool
    operator== ( const MakeVectorialExpression< Scalar, Vectorial > &a,
                 const MakeVectorialExpression< Scalar, Vectorial > &b )
    {
      return ((a.component() == b.component()) && (a.scalar() == b.scalar()));
    }

    template< class Scalar, class Vectorial >
    inline bool
    operator!= ( const MakeVectorialExpression< Scalar, Vectorial > &a,
                 const MakeVectorialExpression< Scalar, Vectorial > &b )
    {
      return ((a.component() != b.component()) || (a.scalar() != b.scalar()));
    }

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

    template< class Scalar, class Vectorial >
    inline void
    axpy ( const typename MakeVectorialTraits< Scalar, Vectorial >::field_type &a,
           const MakeVectorialExpression< Scalar, Vectorial > &x,
           typename MakeVectorialTraits< Scalar, Vectorial >::VectorialType &y )
    {
      typedef MakeVectorialTraits< Scalar, Vectorial > Traits;
      auto z = Traits::component( y, x.component() );
      axpy( a, x.scalar(), z );
    }



    // jacobianTransformation
    // ----------------------

    template< class GJIT, class K, int ROWS, int NROWS >
    inline std::enable_if_t< NROWS % ROWS == 0 >
    jacobianTransformation ( const GJIT &gjit,
                             const MakeVectorialExpression< FieldMatrix< K, ROWS, GJIT::cols >, FieldMatrix< K, NROWS, GJIT::cols > > &a,
                             FieldMatrix< K, NROWS, GJIT::rows > &b )
    {
      typedef MakeVectorialTraits< FieldMatrix< K, ROWS, GJIT::rows >, FieldMatrix< K, NROWS, GJIT::rows > > Traits;
      b = Traits::zeroVectorial();
      for( int i = 0; i < ROWS; ++i )
        gjit.mv( a.scalar()[ i ], Traits::component( b, a.component() )[ i ] );
    }



    // hessianTransformation
    // ---------------------

    template< class GJIT, class K, int SIZE, int NSIZE >
    inline std::enable_if_t< NSIZE % SIZE == 0 >
    hessianTransformation ( const GJIT &gjit,
                            const MakeVectorialExpression< FieldVector< FieldMatrix< K, GJIT::cols, GJIT::cols >, SIZE >, FieldVector< FieldMatrix< K, GJIT::cols, GJIT::cols >, NSIZE > > &a,
                            FieldVector< FieldMatrix< K, GJIT::rows, GJIT::rows >, NSIZE > &b )
    {
      const int dimLocal = GJIT::cols;
      const int dimGlobal = GJIT::rows;
      typedef MakeVectorialTraits< FieldVector< FieldMatrix< K, dimGlobal, dimGlobal >, SIZE >, FieldVector< FieldMatrix< K, dimGlobal, dimGlobal >, NSIZE > > Traits;

      b = Traits::zeroVectorial();

      // c = J^{-T} a_r^T
      // FieldMatrix< K, dimLocal, dimGlobal > c;
      FieldVector< FieldMatrix< K, dimLocal, dimGlobal >, SIZE > c;
      for( int i = 0; i < dimLocal; ++i )
      {
        for( int j = 0; j < SIZE; ++j )
          gjit.mv( a.scalar()[ j ][ i ], c[ j ][ i ] );
      }

      // b_r = J^{-T} c
      for( int i = 0; i < dimGlobal; ++i )
      {
        for( int j = 0; j < SIZE; ++j )
        {
          FieldMatrixColumn< const FieldMatrix< K, dimLocal, dimGlobal > > cji( c[ j ], i );
          FieldMatrixColumn< FieldMatrix< K, dimGlobal, dimGlobal > > bji( Traits::component( b, a.component() )[ j ], i );
          gjit.umv( cji, bji );
        }
      }
    }



    // scalarProduct
    // -------------

    template< class Scalar, class Vectorial >
    inline typename MakeVectorialTraits< Scalar, Vectorial >::field_type
    scalarProduct ( const MakeVectorialExpression< Scalar, Vectorial > &a,
                    const Vectorial &b )
    {
      return scalarProduct( a.scalar()[ 0 ], b[ a.component() ] );
    }

    template< class Scalar, class Vectorial >
    inline typename MakeVectorialTraits< Scalar, Vectorial >::field_type
    scalarProduct ( const Vectorial &a,
                    const MakeVectorialExpression< Scalar, Vectorial > &b )
    {
      return scalarProduct( a[ b.component() ], b.scalar()[ 0 ] );
    }

    template< class Scalar, class Vectorial >
    inline typename MakeVectorialTraits< Scalar, Vectorial >::field_type
    scalarProduct ( const MakeVectorialExpression< Scalar, Vectorial > &a,
                    const MakeVectorialExpression< Scalar, Vectorial > &b )
    {
      typedef typename MakeVectorialTraits< Scalar, Vectorial >::field_type field_type;
      return (a.component() == b.component() ? scalarProduct( a.scalar(), b.scalar() ) : field_type( 0 ));
    }



    // ToNewRange
    // ----------

    template< class ScalarFunctionSpace, class RangeVector >
    struct ToNewRange;

    template< class DomainField, class RangeField, int dimD, int dimR, int dimNR >
    struct ToNewRange< FunctionSpace< DomainField, RangeField, dimD, dimR >, FieldVector< RangeField, dimNR > >
    {
      typedef FunctionSpace< DomainField, RangeField, dimD, dimNR > Type;
    };



    // VectorialShapeFunctionSet
    // -------------------------

    template< class ScalarShapeFunctionSet, class RangeVector >
    class VectorialShapeFunctionSet
    {
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector > ThisType;

    public:
      typedef ScalarShapeFunctionSet ScalarShapeFunctionSetType;

    protected:
      typedef typename ScalarShapeFunctionSetType::FunctionSpaceType ScalarFunctionSpaceType;

      static const std::size_t dimRangeFactor = MakeVectorialTraits< typename ScalarFunctionSpaceType::RangeType, RangeVector >::factor;

      template< class Functor, class Vectorial >
      struct VectorialFunctor;

    public:
      typedef typename ToNewRange< ScalarFunctionSpaceType, RangeVector >::Type FunctionSpaceType;

      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      template< class ... Args >
      VectorialShapeFunctionSet ( Args&& ... args )
      : scalarShapeFunctionSet_( std::forward< Args > ( args ) ... )
      {}

      explicit VectorialShapeFunctionSet ( const ScalarShapeFunctionSetType &scalarShapeFunctionSet )
      : scalarShapeFunctionSet_( scalarShapeFunctionSet )
      {}

      const ScalarShapeFunctionSetType &scalarShapeFunctionSet () const { return scalarShapeFunctionSet_; }

      int order () const { return scalarShapeFunctionSet().order(); }

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
    template< class Functor, class Vectorial >
    struct VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector >::VectorialFunctor
    {
      explicit VectorialFunctor ( const Functor &functor )
      : functor_( functor )
      {}

      template< class Scalar >
      void operator() ( const std::size_t i, const Scalar &value )
      {
        typedef MakeVectorialTraits< Scalar, Vectorial > Traits;
        typedef MakeVectorialExpression< Scalar, Vectorial > Expression;
        const typename Traits::IndexType end = Traits::end();
        for( typename Traits::IndexType k = Traits::begin(); k != end; ++k )
          functor_( i*Traits::factor + Traits::index( k ), Expression( k, value ) );
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
      typedef typename FunctionSpaceType::RangeType VectorialType;
      scalarShapeFunctionSet().evaluateEach( x, VectorialFunctor< Functor, VectorialType >( functor ) );
    }


    template< class ScalarShapeFunctionSet, class RangeVector >
    template< class Point, class Functor >
    inline void VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector >
      ::jacobianEach ( const Point &x, Functor functor ) const
    {
      typedef typename FunctionSpaceType::JacobianRangeType VectorialType;
      scalarShapeFunctionSet().jacobianEach( x, VectorialFunctor< Functor, VectorialType >( functor ) );
    }


    template< class ScalarShapeFunctionSet, class RangeVector >
    template< class Point, class Functor >
    inline void VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector >
      ::hessianEach ( const Point &x, Functor functor ) const
    {
      typedef typename FunctionSpaceType::HessianRangeType VectorialType;
      scalarShapeFunctionSet().hessianEach( x, VectorialFunctor< Functor, VectorialType >( functor ) );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_VECTORIAL_HH

#ifndef DUNE_FEM_COMBINEDBASISFUNCTIONSET_HH
#define DUNE_FEM_COMBINEDBASISFUNCTIONSET_HH


//- dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>


//- dune-fem includes
#include <dune/fem/space/basisfunctionset/basisfunctionset.hh>
#include <dune/fem/misc/subobjects.hh>

namespace Dune
{

  namespace Fem
  {

    // forward declaration
    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    class CombinedBasisFunctionSet;


    //! Traits class for a combined BasisFunctionSetType
    template <class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    struct CombinedBasisFunctionSetTraits
    {
      typedef CombFunctSpace              FunctionSpaceType;
      typedef BasisSetType1                BasisFunctionSetType1;
      typedef BasisSetType2                BasisFunctionSetType2;

      typedef typename BasisFunctionSetType1 :: FunctionSpaceType    FunctionSpaceType1;
      typedef typename BasisFunctionSetType2 :: FunctionSpaceType    FunctionSpaceType2;

      typedef typename FunctionSpaceType1 :: DomainType         DomianType;
      typedef typename FunctionSpaceType1 :: DomainFieldType    DomainFieldType;

      typedef typename FunctionSpaceType1 :: RangeFieldType     RangeFieldType1;
      typedef typename FunctionSpaceType1 :: RangeType          RangeType1;
      typedef typename FunctionSpaceType1 :: JacobianRangeType  JacobianRangeType1;
      typedef typename FunctionSpaceType1 :: HessianRangeType   HessianRangeType1;

      typedef typename FunctionSpaceType2 :: RangeFieldType     RangeFieldType2;
      typedef typename FunctionSpaceType2 :: RangeType          RangeType2;
      typedef typename FunctionSpaceType2 :: JacobianRangeType  JacobianRangeType2;
      typedef typename FunctionSpaceType2 :: HessianRangeType   HessianRangeType2;

      typedef CombinedBasisFunctionSet< FunctionSpaceType, BasisFunctionSetType1, BasisFunctionSetType2 >
        BasisFunctionSetType;
    };

    //! CombinedBasisFunctionSet
    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    class CombinedBasisFunctionSet
    {
      public:
      typedef CombinedBasisFunctionSetTraits< CombFunctSpace, BasisSetType1, BasisSetType2> Traits;
      typedef CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
        ThisType;

      typedef typename Traits :: FunctionSpaceType        FunctionSpaceType;

      typedef typename FunctionSpaceType :: DomainType        DomainType;
      typedef typename FunctionSpaceType :: DomainFieldType   DomainFieldType;
      typedef typename FunctionSpaceType :: RangeType         RangeType;
      typedef typename FunctionSpaceType :: RangeFieldType    RangeFieldType;
      typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType :: HessianRangeType  HessianRangeType;

      typedef typename Traits :: BasisFunctionSetType1     BasisFunctionSetType1;
      typedef typename Traits :: BasisFunctionSetType2     BasisFunctionSetType2;


      typedef typename BasisFunctionSetType1 :: EntityType EntityType;
      typedef typename BasisFunctionSetType1 :: ReferenceElementType ReferenceElementType;


      typedef typename BasisFunctionSetType1 :: FunctionSpaceType    FunctionSpaceType1;
      typedef typename BasisFunctionSetType2 :: FunctionSpaceType    FunctionSpaceType2;

      typedef typename FunctionSpaceType1 :: RangeFieldType     RangeFieldType1;
      typedef typename FunctionSpaceType1 :: RangeType          RangeType1;
      typedef typename FunctionSpaceType1 :: JacobianRangeType  JacobianRangeType1;
      typedef typename FunctionSpaceType1 :: HessianRangeType   HessianRangeType1;

      typedef typename FunctionSpaceType2 :: RangeFieldType     RangeFieldType2;
      typedef typename FunctionSpaceType2 :: RangeType          RangeType2;
      typedef typename FunctionSpaceType2 :: JacobianRangeType  JacobianRangeType2;
      typedef typename FunctionSpaceType2 :: HessianRangeType   HessianRangeType2;
      //! dimension of domain
      enum { dimDomain = FunctionSpaceType :: dimDomain };
      //! dimension of range
      enum { dimRange  = FunctionSpaceType :: dimRange };
      //! dimensions of the two ranges
      enum { dimRange1 = FunctionSpaceType1 :: dimRange };
      enum { dimRange2 = FunctionSpaceType2 :: dimRange };

      public:
      //! constructor
      CombinedBasisFunctionSet( const BasisFunctionSetType1 &basisSet1, const BasisFunctionSetType2 &basisSet2 )
      : basisSet1_( basisSet1 ),
        basisSet2_( basisSet2 ),
        size1_( basisSet1.size() ),
        size2_( basisSet2.size() ),
        offset_( basisSet1_.size() ),
        phi1_( size1_ ),
        phi2_( size2_ ),
        dPhi1_( size1_ ),
        dPhi2_( size2_ )
      {}

      //! HACK for LocalMatrixDefault interface  !//
      CombinedBasisFunctionSet( )
      : basisSet1_(),
        basisSet2_(),
        size1_( 0 ),
        size2_( 0 ),
        offset_( 0 ),
        phi1_( size1_ ),
        phi2_( size2_ )
      {}

      //! \copydoc BasisFunctionSet::size
      std::size_t size()  const
      {
        return size1_ + size2_;
      }

      //! \copydoc BasisFunctionSet::order
      int order ()  const
      {
        return std::max( basisSet1_.order(), basisSet2_.order() );
      }

      //! \copydoc BasisFunctionSet::type
      Dune::GeometryType type () const
      {
        assert( basisSet1_.type() == basisSet2_.type() );
        return basisSet1_.type();
      }

      //! \copydoc BasisFunctionSet::entity
      const EntityType &entity () const
      {
        assert( basisSet1_.entity() == basisSet2_.entity() );
        return basisSet1_.entity();
      }

      //! \copydoc BasisFunctionSet::entity
      const ReferenceElementType &referenceElement () const
      {
        assert( basisSet1_.referenceElement() == basisSet2_.referenceElement() );
        return basisSet1_.referenceElement();
      }

      //! \copydoc BasisFunctionSet::evaluateAll( x, dofs, value )
      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const;

      //! \copydoc BasisFunctionSet::evaluateAll( x, values )
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const;

      //! \copydoc BasisFunctionSet::evaluateAll( quad, dofs, ranges )
      template< class QuadratureType, class DofVector, class RangeArray >
      void evaluateAll ( const QuadratureType &quad, const DofVector &dofs, RangeArray &ranges ) const;


      //! \copydoc BasisFunctionSet::jacobianAll( x, dofs, jacobian )
      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const;

      //! \copydoc BasisFunctionSet::jacobianAll( x, dofs, jacobians )
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const;

      //! \brief evaluate the jacobian of all basis functions and store the result in the jacobians array
      template< class QuadratureType, class DofVector, class JacobianArray >
      void jacobianAll ( const QuadratureType &quad, const DofVector &dofs, JacobianArray &jacobians ) const;


      //! \copydoc BasisFunctionSet::hessianAll( x, dofs, hessian )
      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const;

      //! \todo please doc me
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const;


      //! \copydoc BasisFunctionSet::axpy( quad, values, dofs )
      template< class QuadratureType, class Vector, class DofVector >
      void axpy ( const QuadratureType &quad, const Vector &values, DofVector &dofs ) const;

      //! \copydoc BasisFunctionSet::axpy( quad, valuesA, valuesB, dofs )
      template< class QuadratureType, class VectorA, class VectorB, class DofVector >
      void axpy ( const QuadratureType &quad, const VectorA &valuesA, const VectorB &valuesB, DofVector &dofs ) const;

      //! \copydoc BasisFunctionSet::axpy( x, valueFactor, dofs )
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const;

      //! \copydoc BasisFunctionSet::axpy( x, jacobianFactor, dofs )
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &valueFactor, DofVector &dofs ) const;

      //! \copydoc BasisFunctionSet::axpy( x, valueFactor, jacobianFactor, dofs )
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const;

      protected:
      BasisFunctionSetType1 basisSet1_;
      BasisFunctionSetType2 basisSet2_;
      std::size_t size1_;
      std::size_t size2_;
      std::size_t offset_;
      mutable std::vector< RangeType1 > phi1_;
      mutable std::vector< RangeType2 > phi2_;
      mutable std::vector< JacobianRangeType1 > dPhi1_;
      mutable std::vector< JacobianRangeType2 > dPhi2_;
    };


    // CombinedBasisFunctionSet::evaluateAll
    // -------------------------------------

    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    template< class Point, class DofVector >
    inline void CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
    :: evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
    {
      assert( offset() == basisSet1().size() );
      typedef typename DofVector :: DofType DofType;
      SubDofVector<const DofVector, DofType > dofs1( dofs, size1_, 0  );
      SubDofVector<const DofVector, DofType > dofs2( dofs, size2_, offset_ );
      RangeType1 value1;
      RangeType2 value2;
      basisSet1_.evaluateAll( x, dofs1, value1 );
      basisSet2_.evaluateAll( x, dofs2, value2 );

      for( int r=0;r<dimRange1;++r)
        value[r] = value1[r];

      for( int r=0;r<dimRange2;++r)
        value[ r + dimRange1 ] = value2[ r ];
    }


    // CombinedBasisFunctionSet::evaluateAll
    // -------------------------------------

    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    template< class Point, class RangeArray >
    inline void CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
    :: evaluateAll ( const Point &x, RangeArray &values ) const
    {
      assert( offset_ == basisSet1_.size() );

      phi1_.resize( size1_ );
      phi2_.resize( size2_ );
      basisSet1_.evaluateAll( x, phi1_ );
      basisSet2_.evaluateAll( x, phi2_ );

      const int size = size1_ + size2_;

      values.clear();
      values.resize( size, RangeType(0) );
      for( size_t i=0;i<size1_;++i)
      {
        for(int r=0;r<dimRange1;++r)
          values[ i ][ r ] = phi1_[ i ][ r ];
      }
      for( size_t i=0;i< size2_;++i)
      {
        for(int r=0;r<dimRange2;++r)
          values[ i+ offset_ ][ r +dimRange1 ] = phi2_[ i ][ r ];
      }
    }

    
    // CombinedBasisFunctionSet::evaluateAll
    // -------------------------------------
   
    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    template< class Quadrature, class DofVector, class RangeArray >
    inline void CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
    :: evaluateAll ( const Quadrature &quad, const DofVector &dofs, RangeArray &ranges ) const
    {
      const int nop = quad.nop();
      for( int qp = 0; qp< nop; ++qp )
        evaluateAll( quad[ qp ], dofs, ranges[ qp ] );
    }



    // CombinedBasisFunctionSet::jacobianAll
    // -------------------------------------

    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    template< class Point, class DofVector >
    inline void CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
    :: jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
    {
      typedef typename DofVector :: DofType DofType;
      SubDofVector<const DofVector, DofType > dofs1( dofs, size1_, 0  );
      SubDofVector<const DofVector, DofType > dofs2( dofs, size2_, offset_ );

      JacobianRangeType1 jacobian1;
      JacobianRangeType2 jacobian2;

      basisSet1_.jacobianAll( x, dofs1, jacobian1 );
      basisSet2_.jacobianAll( x, dofs2, jacobian2 );

      for( int r=0;r<dimRange1;++r)
        jacobian[r] = jacobian1[ r ];

      for( int r=0;r<dimRange2;++r)
        jacobian[ r + dimRange1 ] = jacobian2[ r ];
    }


    // CombinedBasisFunctionSet::jacobianAll
    // -------------------------------------

    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    template< class Point, class GlobalJacobianRangeArray >
    inline void CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
    :: jacobianAll ( const Point &x, GlobalJacobianRangeArray &jacobians ) const
    {
      dPhi1_.resize( size1_, JacobianRangeType1( 0.) );
      basisSet1_.jacobianAll( x, dPhi1_ );

      dPhi2_.resize( size2_, JacobianRangeType2( 0.) );
      basisSet2_.jacobianAll( x, dPhi2_ );
      const int size = size1_ + size2_;

      jacobians.resize( size, JacobianRangeType( 0. ) );

      for( size_t i=0;i<size1_;++i)
      {
        for(int r=0;r<dimRange1;++r)
          jacobians[ i ][ r ] =  dPhi1_[i][r];
      }

      for( size_t i=0;i<size2_;++i)
      {
        for(int r=0;r<dimRange2;++r)
          jacobians[ i +offset_ ][ r +dimRange1 ]  = dPhi2_[i][ r ];
      }
    }


    // CombinedBasisFunctionSet::jacobianAll
    // -------------------------------------

    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    template< class Quadrature, class DofVector, class JacobianArray >
    inline void CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
    :: jacobianAll ( const Quadrature &quad, const DofVector &dofs, JacobianArray &jacobians ) const
    {
      const int nop = quad.nop();
      for( int qp = 0; qp< nop; ++qp )
        jacobianAll( quad[ qp ], dofs, jacobians[ qp ] );
    }



    // CombinedBasisFunctionSet::hessianAll
    // -------------------------------------

    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    template< class Point, class DofVector >
    inline void CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
    :: hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
    {
      typedef typename DofVector :: DofType DofType;
      SubDofVector<const DofVector, DofType > dofs1( dofs, size1_, 0  );
      SubDofVector<const DofVector, DofType > dofs2( dofs, size2_, offset_ );

      HessianRangeType1 hessian1;
      HessianRangeType2 hessian2;

      basisSet1_.hessianAll( x, dofs1, hessian1 );
      basisSet2_.hessianAll( x, dofs2, hessian2 );

      for( int r=0;r<dimRange1;++r)
        hessian[ r ] = hessian1[ r ];

      for( int r=0;r<dimRange2;++r)
        hessian[ r + dimRange1 ] = hessian2[ r ];
    }


    // CombinedBasisFunctionSet::hessianAll
    // -------------------------------------

    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    template< class Point, class HessianRangeArray >
    inline void CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
    :: hessianAll ( const Point &x, HessianRangeArray &hessians ) const
    {
      std::vector< HessianRangeType1 > phi1( size1_ );
      basisSet1_.hessianAll( x, phi1 );

      std::vector< HessianRangeType2 > phi2( size2_ );
      basisSet2_.hessianAll( x, phi2 );

      const int size = size1_ + size2_;

      hessians.resize( size );

      for( size_t i=0;i<size1_;++i)
      {
        for(int r=0;r<dimRange1;++r)
          hessians[ i ][ r ] =  phi1[i][r];
      }

      for( size_t i=0;i<size2_;++i)
      {
        for(int r=0;r<dimRange2;++r)
          hessians[ i +offset_ ][ r +dimRange1 ]  = phi2[i][ r ];
      }
    }



    // CombinedBasisFunctionSet::axpy
    // ------------------------------

    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    template< class Point, class DofVector >
    inline void CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
    :: axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
    {
      typedef typename DofVector :: DofType DofType;
      SubDofVector<DofVector, DofType > dofs1( dofs, size1_, 0  );
      SubDofVector<DofVector, DofType > dofs2( dofs, size2_, offset_ );
      SubObject< const RangeType, const RangeType1, 0 > valueFactor1( valueFactor );
      SubObject< const RangeType, const RangeType2, dimRange1 > valueFactor2( valueFactor );

      basisSet1_.axpy(x, (RangeType1) valueFactor1, dofs1);
      basisSet2_.axpy(x, (RangeType2) valueFactor2, dofs2);
    }


    // CombinedBasisFunctionSet::axpy
    // ------------------------------

    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    template< class Point, class DofVector >
    inline void CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
    :: axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
    {
      typedef typename DofVector :: DofType DofType;
      SubDofVector<DofVector, DofType > dofs1( dofs, size1_, 0  );
      SubDofVector<DofVector, DofType > dofs2( dofs, size2_, offset_ );
      SubObject< const JacobianRangeType, const JacobianRangeType1, 0 > jacobianFactor1( jacobianFactor );
      SubObject< const JacobianRangeType, const JacobianRangeType2, dimRange1 > jacobianFactor2( jacobianFactor );

      basisSet1_.axpy(x, (JacobianRangeType1) jacobianFactor1, dofs1);
      basisSet2_.axpy(x, (JacobianRangeType2) jacobianFactor2, dofs2);
    }


    // CombinedBasisFunctionSet::axpy
    // ------------------------------

    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    template< class Point, class DofVector >
    inline void CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
    :: axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
              DofVector &dofs ) const
    {
      typedef typename DofVector :: DofType DofType;
      SubDofVector<DofVector, DofType > dofs1( dofs, size1_, 0  );
      SubDofVector<DofVector, DofType > dofs2( dofs, size2_, offset_ );
      SubObject< const RangeType, const RangeType1, 0 > valueFactor1( valueFactor );
      SubObject< const RangeType, const RangeType2, dimRange1 > valueFactor2( valueFactor );
      SubObject< const JacobianRangeType, const JacobianRangeType1, 0 > jacobianFactor1( jacobianFactor );
      SubObject< const JacobianRangeType, const JacobianRangeType2, dimRange1 > jacobianFactor2( jacobianFactor );

      basisSet1_.axpy(x, (RangeType1) valueFactor1, (JacobianRangeType1) jacobianFactor1, dofs1);
      basisSet2_.axpy(x, (RangeType1) valueFactor2, (JacobianRangeType2) jacobianFactor2, dofs2);
    }


    // CombinedBasisFunctionSet::axpy
    // ------------------------------

    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    template< class Quadrature, class Vector, class DofVector >  
    inline void CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
    :: axpy ( const Quadrature &quad, const Vector &values, DofVector &dofs ) const
    {
      const int nop = quad.nop();
      for( int qp = 0; qp< nop; ++qp )
        axpy( quad[ qp ], values[ qp ], dofs );
    }


    // CombinedBasisFunctionSet::axpy
    // ------------------------------

    template<class CombFunctSpace, class BasisSetType1, class BasisSetType2>
    template< class Quadrature, class VectorA, class VectorB, class DofVector >  
    inline void CombinedBasisFunctionSet< CombFunctSpace, BasisSetType1, BasisSetType2>
    :: axpy ( const Quadrature &quad, const VectorA &valuesA, const VectorB &valuesB, DofVector &dofs ) const
    {
      const int nop = quad.nop();
      for( int qp = 0; qp< nop; ++qp )
      {
        axpy( quad[ qp ], valuesA[ qp ], dofs );
        axpy( quad[ qp ], valuesB[ qp ], dofs );
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMBINEDBASISFUNCTIONSET_HH

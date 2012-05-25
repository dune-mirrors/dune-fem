#ifndef COMBINEDBASEFUNCTIONS_HH
#define COMBINEDBASEFUNCTIONS_HH

#include <dune/fem/space/basefunctions/basefunctionsetinterface.hh>

#include <dune/fem/misc/subobjects.hh>

namespace Dune
{

  // forward declaration
  template<class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  class CombinedBaseFunctionSet;


  //! Traits class for a combined BaseFunctionSetType
  template <class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  struct CombinedBaseFunctionSetTraits
  { 
    typedef CombFunctSpace              FunctionSpaceType;
    typedef BaseSetType1                BaseFunctionSetType1;
    typedef BaseSetType2                BaseFunctionSetType2;

    typedef typename BaseFunctionSetType1 :: FunctionSpaceType    FunctionSpaceType1;
    typedef typename BaseFunctionSetType2 :: FunctionSpaceType    FunctionSpaceType2;

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
      
    typedef CombinedBaseFunctionSet< FunctionSpaceType, BaseFunctionSetType1, BaseFunctionSetType2 > 
      BaseFunctionSetType;
  };

  //! CombinedBaseFunctionSet
  template<class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  class CombinedBaseFunctionSet 
  : public BaseFunctionSetDefault< CombinedBaseFunctionSetTraits< CombFunctSpace, BaseSetType1, BaseSetType2> >
  {
    public:
    typedef CombinedBaseFunctionSetTraits< CombFunctSpace, BaseSetType1, BaseSetType2> Traits;
    typedef CombinedBaseFunctionSet< CombFunctSpace, BaseSetType1, BaseSetType2>
      ThisType;

    typedef typename Traits :: FunctionSpaceType        FunctionSpaceType;

    typedef typename FunctionSpaceType :: DomainType        DomainType;
    typedef typename FunctionSpaceType :: DomainFieldType   DomainFieldType;
    typedef typename FunctionSpaceType :: RangeType         RangeType;
    typedef typename FunctionSpaceType :: RangeFieldType    RangeFieldType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType :: HessianRangeType  HessianRangeType;
    
    typedef typename Traits :: BaseFunctionSetType1     BaseFunctionSetType1;
    typedef typename Traits :: BaseFunctionSetType2     BaseFunctionSetType2;


    typedef typename BaseFunctionSetType1 :: FunctionSpaceType    FunctionSpaceType1;
    typedef typename BaseFunctionSetType2 :: FunctionSpaceType    FunctionSpaceType2;

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
    CombinedBaseFunctionSet( const BaseFunctionSetType1 &baseSet1, const BaseFunctionSetType2 &baseSet2 ) 
      : baseSet1_( baseSet1 ),
        baseSet2_( baseSet2 ),
        size1_( baseSet1.size() ),
        size2_( baseSet2.size() ), 
        offset_( baseSet1_.size() )
    {}

    //! HACK for LocalMatrixDefault interface  !//
    CombinedBaseFunctionSet( )
      : baseSet1_(), 
        baseSet2_(),
        size1_( 0 ),
        size2_( 0 ),
        offset_( 0 )
    {}

    size_t size()  const
    {
      return size1_ + size2_; 
    }

    GeometryType geometryType () const
    {
      assert( baseSet1_.geometryType() == baseSet2_.geometryType() );
      return baseSet1_.geometryType();
    }

    template< int diffOrd, class PointType >
    DUNE_DEPRECATED
    void evaluate ( const int baseFunction,
                    const FieldVector< int, diffOrd> &diffVariable,
                    const PointType &x,
                    RangeType &phi ) const;

    template< class Point >
    DUNE_DEPRECATED
    void evaluate ( const int baseFunction, const Point &x, RangeType &value ) const;

    template< int diffOrder, class Point, class DofVector >
    void evaluateAll ( const FieldVector< int, diffOrder > &diffVariable,
                       const Point &x, const DofVector &dofs, RangeType &value ) const;
    template< class Point, class DofVector >
    void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const;

    template< class Point, class RangeArray >
    void evaluateAll ( const Point &x, RangeArray &values ) const;

    // jacobian methods
    template< class Point >
    DUNE_DEPRECATED
    void jacobian ( const int baseFunction, const Point &x, JacobianRangeType &jacobian ) const;

    template< class Point, class GeometryJacobianInverse, class DofVector, class GlobalJacobianRange >
    void jacobianAll ( const Point &x, const GeometryJacobianInverse &gjit,
                       const DofVector &dofs, GlobalJacobianRange &jacobian ) const;

    template< class Point, class GeometryJacobianInverse, class GlobalJacobianRangeArray >
    void jacobianAll ( const Point &x, const GeometryJacobianInverse &gjit,
                       GlobalJacobianRangeArray &jacobians ) const;


    // axpy methods
    template< class Point, class DofVector >
    void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const;

    template< class Point, class GeometryJacobianInverse, class GlobalJacobianRange, class DofVector >
    void axpy ( const Point &x, const GeometryJacobianInverse &gjit,
                const GlobalJacobianRange &jacobianFactor, DofVector &dofs ) const;
    template< class Point, class GeometryJacobianInverse, class GlobalJacobianRange, class DofVector >
    void axpy ( const Point &x, const GeometryJacobianInverse &gjit,
                const RangeType &valueFactor, const GlobalJacobianRange &jacobianFactor,
                DofVector &dofs ) const;
      
    protected:
    BaseFunctionSetType1 baseSet1_;
    BaseFunctionSetType2 baseSet2_;
    std::size_t size1_;
    std::size_t size2_;
    std::size_t offset_;   
  };
    

  template<class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  template< int diffOrd, class PointType >
  DUNE_DEPRECATED
  inline void CombinedBaseFunctionSet< CombFunctSpace, BaseSetType1, BaseSetType2>
  :: evaluate ( const int baseFunction,
                const FieldVector< int, diffOrd> &diffVariable,
                const PointType &x,
                RangeType &phi ) const
  {
    assert( offset_ == size1_ );

    if( baseFunction - offset_ < 0 )
    {
      RangeType1 phi1;
      baseSet1_.evaluate( baseFunction, diffVariable, x, phi1);     

      for(int r=0;r<dimRange1;++r)
        phi[ r ] = phi1[ r ];
      for(int r=0;r<dimRange2;++r)
        phi[ r +dimRange1 ]  = 0;
    }
    else
    {
      RangeType2 phi2;
      baseSet2_.evaluate( baseFunction - offset_,  diffVariable, x, phi2);

      for(int r=0;r<dimRange1;++r)
        phi[ r ] = 0;
      for(int r=0;r<dimRange2;++r)
        phi[ r +dimRange1 ]  = phi2[ r ];
    }
  }

  template<class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  template< class PointType >
  DUNE_DEPRECATED
  inline void CombinedBaseFunctionSet< CombFunctSpace, BaseSetType1, BaseSetType2>
  :: evaluate ( const int baseFunction,
                const PointType &x,
                RangeType &phi ) const
  {
    assert( offset_ == baseSet1_.size() );
    
    if( baseFunction - offset_ < 0 )
    {
      RangeType1 phi1;
      baseSet1_.evaluate( baseFunction, x, phi1);        
      for(int r=0;r<dimRange1;++r)
        phi[ r ] = phi1[ r ];
      for(int r=0;r<dimRange2;++r)
        phi[ r +dimRange1 ]  = 0;
    }
    else
    {
      RangeType2 phi2;
      baseSet2_.evaluate( baseFunction - offset_, x, phi2);
      for(int r=0;r<dimRange1;++r)
        phi[ r ] = 0;
      for(int r=0;r<dimRange2;++r)
        phi[ r +dimRange1 ]  = phi2[ r ];
    }
  }

  template<class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  template< int diffOrder, class Point, class DofVector >
  inline void CombinedBaseFunctionSet< CombFunctSpace, BaseSetType1, BaseSetType2>
  :: evaluateAll ( const FieldVector< int, diffOrder > &diffVariable,
                   const Point &x, const DofVector &dofs, RangeType &value ) const
  {
    SubDofVector<const DofVector, double > dofs1( dofs, size1_, 0  );
    SubDofVector<const DofVector, double > dofs2( dofs, size2_, offset_ );
    RangeType1 value1;
    RangeType2 value2;
    baseSet1_.evaluateAll( diffVariable, x, dofs1, value1);        
    baseSet2_.evaluateAll( diffVariable, x, dofs2, value2);       

    for( int r=0;r<dimRange1;++r)
      value[r] = value1[r];

    for( int r=0;r<dimRange2;++r)
      value[ r + dimRange1 ] = value2[ r ];
  }

  template<class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  template< class Point, class RangeArray >
  inline void CombinedBaseFunctionSet< CombFunctSpace, BaseSetType1, BaseSetType2>
  :: evaluateAll ( const Point &x, RangeArray &values ) const
  {
    std::vector< RangeType1 > phi1;
    baseSet1_.evaluateAll( x, phi1 );        

    std::vector< RangeType2 > phi2;
    baseSet2_.evaluateAll( x, phi2 );      

    const int size = size1_ + size2_;
    values.resize( size );
    for( size_t i=0;i<size1_;++i)
    {
      for(int r=0;r<dimRange1;++r)
        values[ i ][ r ] =  phi1[i][r];
    }
    for( size_t i=0;i< size2_;++i)
    {
      for(int r=0;r<dimRange2;++r)
        values[ i+ offset_ ][ r +dimRange1 ]  = phi2[i][ r ];
    }
  }

  template<class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  template< class Point, class DofVector >
  inline void CombinedBaseFunctionSet< CombFunctSpace, BaseSetType1, BaseSetType2>
  :: evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
  {
    SubDofVector< const DofVector, double > dofs1( dofs, size1_, 0  );
    SubDofVector< const DofVector, double > dofs2( dofs, size2_, offset_ );
    RangeType1 value1;
    RangeType2 value2;
    baseSet1_.evaluateAll( x, dofs1, value1);        
    baseSet2_.evaluateAll( x, dofs2, value2);       

    for( int r=0;r<dimRange1;++r)
      value[r] = value1[r];

    for( int r=0;r<dimRange2;++r)
      value[ r + dimRange1 ] = value2[ r ];
  }


  template<class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  template< class PointType >
  DUNE_DEPRECATED
  inline void CombinedBaseFunctionSet< CombFunctSpace, BaseSetType1, BaseSetType2>
  :: jacobian ( const int baseFunction,
                const PointType &x,
                JacobianRangeType &phi ) const
  {
    assert( offset_ == baseSet1_.size() );

    if( baseFunction - offset_ < 0 )
    {
      JacobianRangeType1 phi1;
      baseSet1_.jacobian( baseFunction, x, phi1);
      for(int r=0;r<dimRange1;++r)
        phi[ r ] = phi1[ r ];
      for(int r=0;r<dimRange2;++r)
        phi[ r +dimRange1 ] = 0;
    }
    else
    {
      JacobianRangeType2 phi2;
      baseSet2_.jacobian( baseFunction - offset_, x, phi2);
      for(int r=0;r<dimRange1;++r)
        phi[ r ] = 0;      
      for(int r=0;r<dimRange2;++r)
        phi[ r +dimRange1 ] = phi2[ r ];
    }
  }

  /** \brief \todo please doc me */
  template<class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  template< class Point, class GeometryJacobianInverse, class DofVector, class GlobalJacobianRange >
  inline void CombinedBaseFunctionSet< CombFunctSpace, BaseSetType1, BaseSetType2>
  :: jacobianAll ( const Point &x, const GeometryJacobianInverse &gjit,
                   const DofVector &dofs, GlobalJacobianRange &jacobian ) const
  {
    SubDofVector<const DofVector, double > dofs1( dofs, size1_, 0  );
    SubDofVector<const DofVector, double > dofs2( dofs, size2_, offset_ );

    JacobianRangeType1 jacobian1; 
    JacobianRangeType2 jacobian2; 

    baseSet1_.jacobianAll( x, gjit, dofs1, jacobian1);        
    baseSet2_.jacobianAll( x, gjit, dofs2, jacobian2);       

    for( int r=0;r<dimRange1;++r)
      jacobian[r] = jacobian1[ r];

    for( int r=0;r<dimRange2;++r)
      jacobian[ r + dimRange1 ] = jacobian2[ r ];
  }

  /** \brief \todo please doc me */
  template<class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  template< class Point, class GeometryJacobianInverse, class GlobalJacobianRangeArray >
  inline void CombinedBaseFunctionSet< CombFunctSpace, BaseSetType1, BaseSetType2>
  :: jacobianAll ( const Point &x, const GeometryJacobianInverse &gjit,
                   GlobalJacobianRangeArray &jacobians ) const
  {
    std::vector< JacobianRangeType1 > phi1;
    baseSet1_.jacobianAll( x, gjit, phi1 );        

    std::vector< JacobianRangeType2 > phi2;
    baseSet2_.jacobianAll( x, gjit, phi2 );       
    const int size = size1_ + size2_;

    jacobians.resize( size );

    for( size_t i=0;i<size1_;++i)
    {
      for(int r=0;r<dimRange1;++r)
        jacobians[ i ][ r ] =  phi1[i][r];
    }

    for( size_t i=0;i<size2_;++i)
    {        
      for(int r=0;r<dimRange2;++r)
        jacobians[ i +offset_ ][ r +dimRange1 ]  = phi2[i][ r ];
    }
  }


  /** \todo please doc me */
  template<class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  template< class Point, class DofVector >
  inline void CombinedBaseFunctionSet< CombFunctSpace, BaseSetType1, BaseSetType2>
  :: axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
  {
    SubDofVector<DofVector, double > dofs1( dofs, size1_, 0  );
    SubDofVector<DofVector, double > dofs2( dofs, size2_, offset_ );
    SubObject< const RangeType, const RangeType1, 0 > valueFactor1( valueFactor ); 
    SubObject< const RangeType, const RangeType2, dimRange1 > valueFactor2( valueFactor );

    baseSet1_.axpy(x, valueFactor1, dofs1);

    baseSet2_.axpy(x, valueFactor2, dofs2);
  }

  /** \todo please doc me */
  template<class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  template< class Point, class GeometryJacobianInverse, class GlobalJacobianRange, class DofVector >
  inline void CombinedBaseFunctionSet< CombFunctSpace, BaseSetType1, BaseSetType2>
  :: axpy ( const Point &x, const GeometryJacobianInverse &gjit,
            const GlobalJacobianRange &jacobianFactor, DofVector &dofs ) const
  {

    SubDofVector<DofVector, double > dofs1( dofs, size1_, 0  );
    SubDofVector<DofVector, double > dofs2( dofs, size2_, offset_ );
    SubObject< const GlobalJacobianRange, const JacobianRangeType1, 0 > jacobianFactor1( jacobianFactor );
    SubObject< const GlobalJacobianRange, const JacobianRangeType2, dimRange1 > jacobianFactor2( jacobianFactor ); 

    baseSet1_.axpy(x, gjit, jacobianFactor1, dofs1);
    baseSet2_.axpy(x, gjit, jacobianFactor2, dofs2);
  }

  /** \todo please doc me */
  template<class CombFunctSpace, class BaseSetType1, class BaseSetType2>
  template< class Point, class GeometryJacobianInverse, class GlobalJacobianRange, class DofVector >
  inline void CombinedBaseFunctionSet< CombFunctSpace, BaseSetType1, BaseSetType2>
  :: axpy ( const Point &x, const GeometryJacobianInverse &gjit,
            const RangeType &valueFactor, const GlobalJacobianRange &jacobianFactor,
            DofVector &dofs ) const
  {
    SubDofVector<DofVector, double > dofs1( dofs, size1_, 0  );
    SubDofVector<DofVector, double > dofs2( dofs, size2_, offset_ );
    SubObject< const RangeType, const RangeType1, 0 > valueFactor1( valueFactor ); 
    SubObject< const RangeType, const RangeType2, dimRange1 > valueFactor2( valueFactor );
    SubObject< const GlobalJacobianRange, const JacobianRangeType1, 0 > jacobianFactor1( jacobianFactor );
    SubObject< const GlobalJacobianRange, const JacobianRangeType2, dimRange1 > jacobianFactor2( jacobianFactor ); 

    baseSet1_.axpy(x, gjit, valueFactor1, jacobianFactor1, dofs1);
    baseSet2_.axpy(x, gjit, valueFactor2, jacobianFactor2, dofs2);
  }

}
#endif // header guards

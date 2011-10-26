#ifndef COMBINEDBASEFUNCTIONS_HH
#define COMBINEDBASEFUNCTIONS_HH

#include <dune/fem/space/basefunctions/basefunctionsetinterface.hh>

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
        offSet_( baseSet1_.numBaseFunctions() )
    {}

    //! HACK for LocalMatrixDefault interface  !//
    CombinedBaseFunctionSet( )
      : baseSet1_(), 
        baseSet2_(),
        offSet_( 0 )
    {}

#if 0
    //! copy constructor
    CombinedBaseFunctionSet( const ThisType &other )
      : baseSet1_( other.baseSet1_ ),
        baseSet2_( other.baseSet2_ ),
        offSet_( other.offSet_ )
    {}
#endif

    int numBaseFunctions()  const
    {
      return baseSet1_.numBaseFunctions() + baseSet2_.numBaseFunctions();            
    }
    
    GeometryType geometryType () const
    {
      assert( baseSet1_.geometryType() == baseSet2_.geometryType() );
      return baseSet1_.geometryType();
    }

    template< int diffOrd, class PointType >
    void evaluate ( const int baseFunction,
                    const FieldVector< int, diffOrd > &diffVariable,
                    const PointType &x,
                    RangeType &phi ) const
    {
      assert( offSet_ == baseSet1_.numBaseFunctions() );

      if( baseFunction - offSet_ < 0 )
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
        baseSet2_.evaluate( baseFunction - offSet_,  diffVariable, x, phi2);

        for(int r=0;r<dimRange1;++r)
          phi[ r ] = 0;
        for(int r=0;r<dimRange2;++r)
          phi[ r +dimRange1 ]  = phi2[ r ];
      }
    }


    template< class PointType >
    void evaluate ( const int baseFunction,
                    const PointType &x,
                    RangeType &phi ) const
    {
      assert( offSet_ == baseSet1_.numBaseFunctions() );
      
      if( baseFunction - offSet_ < 0 )
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
        baseSet2_.evaluate( baseFunction - offSet_, x, phi2);
        for(int r=0;r<dimRange1;++r)
          phi[ r ] = 0;
        for(int r=0;r<dimRange2;++r)
          phi[ r +dimRange1 ]  = phi2[ r ];
      }
    }


    template< class PointType >
    void jacobian ( const int baseFunction,
                    const PointType &x,
                    JacobianRangeType &phi ) const
    {
      assert( offSet_ == baseSet1_.numBaseFunctions() );

      if( baseFunction - offSet_ < 0 )
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
        baseSet2_.jacobian( baseFunction - offSet_, x, phi2);
        for(int r=0;r<dimRange1;++r)
          phi[ r ] = 0;      
        for(int r=0;r<dimRange2;++r)
          phi[ r +dimRange1 ] = phi2[ r ];
      }
    }


    protected:
    BaseFunctionSetType1 baseSet1_;
    BaseFunctionSetType2 baseSet2_;
    int offSet_;

  };

}
#endif // header guards

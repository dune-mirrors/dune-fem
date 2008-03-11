#ifndef DUNE_LAGRANGESPACE_BASEFUNCTIONS_HH
#define DUNE_LAGRANGESPACE_BASEFUNCTIONS_HH

#include "genericbasefunctions.hh"
#include "lagrangepoints.hh"

namespace Dune
{
   
  template< class FunctionSpace, GeometryType :: BasicType type,
            unsigned int dim,  unsigned int pOrder >
  class LagrangeBaseFunction
  : public BaseFunctionInterface< FunctionSpace >
  {
    typedef LagrangeBaseFunction< FunctionSpace, type, dim, pOrder > ThisType;
    typedef BaseFunctionInterface< FunctionSpace > BaseType;

    typedef typename GeometryWrapper< type, dim > :: GenericGeometryType
      GenericGeometryType;
    typedef GenericLagrangeBaseFunction
      < FunctionSpace, GenericGeometryType, pOrder >
      GenericBaseFunctionType;
      
  public:
    enum { dimension = GenericGeometryType :: dimension };
    
    enum { polynomialOrder = GenericBaseFunctionType :: polynomialOrder };

    enum { numBaseFunctions = GenericBaseFunctionType :: numBaseFunctions };

    typedef typename GenericBaseFunctionType :: DomainType DomainType;
    typedef typename GenericBaseFunctionType :: RangeType RangeType;

    typedef typename GenericBaseFunctionType :: DomainFieldType DomainFieldType;
    typedef typename GenericBaseFunctionType :: RangeFieldType RangeFieldType;

    typedef LagrangePoint< type, dim, polynomialOrder > LagrangePointType;
    typedef LagrangePointListImplementation
      < DomainFieldType, type, dim, polynomialOrder >
      LagrangePointListType;

  protected:
    GenericBaseFunctionType baseFunction_;

  public:
    inline LagrangeBaseFunction( unsigned int baseNum )
    : BaseType(),
      baseFunction_( baseNum )
    {}

    virtual void evaluate ( const FieldVector< deriType, 0 > &diffVariable,
                            const DomainType &x,
                            RangeType &phi ) const
    {
      baseFunction_.evaluate( diffVariable, x, phi );
    }
    
    virtual void evaluate ( const FieldVector< deriType, 1 > &diffVariable,
                            const DomainType &x,
                            RangeType &phi ) const
    {
      baseFunction_.evaluate( diffVariable, x, phi );
    }

    virtual void evaluate ( const FieldVector< deriType, 2 > &diffVariable,
                            const DomainType &x,
                            RangeType &phi ) const
    {
      baseFunction_.evaluate( diffVariable, x, phi );
    }

    virtual int order () const
    {
      return polynomialOrder;
    }
  };



  template< class ScalarFunctionSpaceImp, unsigned int dim, unsigned int order >
  class LagrangeBaseFunctionFactory
  : public BaseFunctionFactory< ScalarFunctionSpaceImp >
  {
  public:
    typedef ScalarFunctionSpaceImp ScalarFunctionSpaceType;

    typedef BaseFunctionInterface< ScalarFunctionSpaceType > BaseFunctionType;

    enum { dimension = dim };

    enum { polynomialOrder = order };

  private:
    typedef BaseFunctionFactory< ScalarFunctionSpaceType > BaseType;
    typedef LagrangeBaseFunctionFactory
      < ScalarFunctionSpaceType, dimension, polynomialOrder >
      ThisType;

  public:
    LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseType( geometry )
    {
      assert( this->geometry().dim() == dimension );
    }

    virtual ~LagrangeBaseFunctionFactory ()
    {
    }

    virtual BaseFunctionType* baseFunction ( int i ) const
    {
      const GeometryType :: BasicType basicType = this->geometry().basicType();

      switch( basicType ) {
      case GeometryType :: simplex:
        return new LagrangeBaseFunction< ScalarFunctionSpaceType,
                                         GeometryType :: simplex, dimension,
                                         polynomialOrder >
                                       ( i );

      case GeometryType :: cube:
        return new LagrangeBaseFunction< ScalarFunctionSpaceType,
                                         GeometryType :: cube, dimension,
                                         polynomialOrder >
                                       ( i );

      default:
        DUNE_THROW( NotImplemented, "No such geometry type implemented." );
      }
    }
    
    virtual int numBaseFunctions () const
    {
      const GeometryType :: BasicType basicType = this->geometry().basicType();

      switch( basicType ) {
      case GeometryType :: simplex:
        return LagrangeBaseFunction< ScalarFunctionSpaceType,
                                     GeometryType :: simplex, dimension,
                                     polynomialOrder >
                 :: numBaseFunctions;

      case GeometryType :: cube:
        return LagrangeBaseFunction< ScalarFunctionSpaceType,
                                     GeometryType :: cube, dimension,
                                     polynomialOrder >
                 :: numBaseFunctions;

      default:
        DUNE_THROW( NotImplemented, "No such geometry type implemented." );
      }
    }
  };



  template< class ScalarFunctionSpaceImp, unsigned int polOrder >
  class LagrangeBaseFunctionFactory< ScalarFunctionSpaceImp, 3, polOrder >
  : public BaseFunctionFactory< ScalarFunctionSpaceImp >
  {
  public:
    typedef ScalarFunctionSpaceImp ScalarFunctionSpaceType;

    typedef BaseFunctionInterface< ScalarFunctionSpaceType > BaseFunctionType;

    enum { dimension = 3 };

    enum { polynomialOrder = polOrder };

  private:
    typedef BaseFunctionFactory< ScalarFunctionSpaceType > BaseType;
    typedef LagrangeBaseFunctionFactory
      < ScalarFunctionSpaceType, dimension, polynomialOrder >
      ThisType;

  public:
    LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseType( geometry )
    {
      assert( this->geometry().dim() == dimension );
    }

    virtual ~LagrangeBaseFunctionFactory ()
    {
    }

    virtual BaseFunctionType* baseFunction ( int i ) const
    {
      const GeometryType :: BasicType basicType = this->geometry().basicType();

      switch( basicType ) {
      case GeometryType :: simplex:
        return new LagrangeBaseFunction< ScalarFunctionSpaceType,
                                         GeometryType :: simplex, dimension,
                                         polynomialOrder >
                                       ( i );

      case GeometryType :: cube:
        return new LagrangeBaseFunction< ScalarFunctionSpaceType,
                                         GeometryType :: cube, dimension,
                                         polynomialOrder >
                                       ( i );

      case GeometryType :: pyramid:
        return new LagrangeBaseFunction< ScalarFunctionSpaceType,
                                         GeometryType :: pyramid, dimension,
                                         polynomialOrder >
                                       ( i );

      case GeometryType :: prism:
        return new LagrangeBaseFunction< ScalarFunctionSpaceType,
                                         GeometryType :: prism, dimension,
                                         polynomialOrder >
                                       ( i );

      default:
        DUNE_THROW( NotImplemented, "No such geometry type implemented." );
      }
    }
    
    virtual int numBaseFunctions () const
    {
      const GeometryType :: BasicType basicType = this->geometry().basicType();

      switch( basicType ) {
      case GeometryType :: simplex:
        return LagrangeBaseFunction< ScalarFunctionSpaceType,
                                     GeometryType :: simplex, dimension,
                                     polynomialOrder >
                 :: numBaseFunctions;

      case GeometryType :: cube:
        return LagrangeBaseFunction< ScalarFunctionSpaceType,
                                     GeometryType :: cube, dimension,
                                     polynomialOrder >
                 :: numBaseFunctions;

      case GeometryType :: pyramid:
        return LagrangeBaseFunction< ScalarFunctionSpaceType,
                                     GeometryType :: pyramid, dimension,
                                     polynomialOrder >
                 :: numBaseFunctions;

      case GeometryType :: prism:
        return LagrangeBaseFunction< ScalarFunctionSpaceType,
                                     GeometryType :: prism, dimension,
                                     polynomialOrder >
                 :: numBaseFunctions;

      default:
        DUNE_THROW( NotImplemented, "No such geometry type implemented." );
      }
    }
  };


  
}

#endif

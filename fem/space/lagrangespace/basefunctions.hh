#ifndef DUNE_LAGRANGESPACE_BASEFUNCTIONS_HH
#define DUNE_LAGRANGESPACE_BASEFUNCTIONS_HH

#include "genericbasefunctions.hh"
#include "lagrangepoints.hh"

namespace Dune
{
   
  template< class FunctionSpaceType,
            GeometryType :: BasicType type,
            unsigned int dim,
            unsigned int order >
  class LagrangeBaseFunction
  : public GenericLagrangeBaseFunction
    < FunctionSpaceType,
      typename GeometryWrapper< type, dim > :: GenericGeometryType,
      order >
  {
  private:
    typedef GenericLagrangeBaseFunction
      < FunctionSpaceType,
        typename GeometryWrapper< type, dim > :: GenericGeometryType,
        order >
      BaseType;
    typedef LagrangeBaseFunction< FunctionSpaceType, type, dim, order >
      ThisType;

  public:
    enum { dimension = BaseType :: GeometryType :: dimension };
    
    enum { polynomialOrder = BaseType :: polynomialOrder };

    enum { numBaseFunctions = BaseType :: numBaseFunctions };

    typedef typename BaseType :: DomainType DomainType;
    typedef typename BaseType :: RangeType RangeType;

    typedef typename BaseType :: DomainFieldType DomainFieldType;
    typedef typename BaseType :: RangeFieldType RangeFieldType;

    typedef LagrangePoint< type, dim, polynomialOrder > LagrangePointType;
    typedef LagrangePointListImplementation< DomainFieldType, type, dim, polynomialOrder >
      LagrangePointListType;

  public:
    inline LagrangeBaseFunction( unsigned int baseNum )
    : BaseType( baseNum )
    {
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

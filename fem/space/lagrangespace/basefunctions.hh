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
    LagrangeBaseFunction( unsigned int baseNum );

    virtual void evaluate ( const FieldVector< deriType, 0 > &diffVariable,
                            const DomainType &x,
                            RangeType &phi ) const;
    
    virtual void evaluate ( const FieldVector< deriType, 1 > &diffVariable,
                            const DomainType &x,
                            RangeType &phi ) const;

    virtual void evaluate ( const FieldVector< deriType, 2 > &diffVariable,
                            const DomainType &x,
                            RangeType &phi ) const;

    virtual int order () const;
  };



  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  class LagrangeBaseFunctionFactory
  : public BaseFunctionFactory< ScalarFunctionSpace >
  {
    typedef LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
      ThisType;
    typedef BaseFunctionFactory< ScalarFunctionSpace > BaseType;

  public:
    typedef ScalarFunctionSpace ScalarFunctionSpaceType;

    typedef BaseFunctionInterface< ScalarFunctionSpaceType > BaseFunctionType;

    enum { dimension = dim };

    enum { polynomialOrder = pOrder };

  public:
    LagrangeBaseFunctionFactory ( GeometryType geometry );

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionType *baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };



  template< class ScalarFunctionSpace, unsigned int pOrder >
  class LagrangeBaseFunctionFactory< ScalarFunctionSpace, 3, pOrder >
  : public BaseFunctionFactory< ScalarFunctionSpace >
  {
    typedef LagrangeBaseFunctionFactory< ScalarFunctionSpace, 3, pOrder >
      ThisType;
    typedef BaseFunctionFactory< ScalarFunctionSpace > BaseType;

  public:
    typedef ScalarFunctionSpace ScalarFunctionSpaceType;

    typedef BaseFunctionInterface< ScalarFunctionSpaceType > BaseFunctionType;

    enum { dimension = 3 };

    enum { polynomialOrder = pOrder };

  public:
    LagrangeBaseFunctionFactory ( GeometryType geometry );

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionType *baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };
  
}

#include "basefunction_inline.hh"

#endif

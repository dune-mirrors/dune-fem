#ifndef DUNE_LAGRANGESPACE_BASEFUNCTIONS_HH
#define DUNE_LAGRANGESPACE_BASEFUNCTIONS_HH

#include "genericbasefunctions.hh"

namespace Dune
{
   
  template< class FunctionSpace, GeometryType :: BasicType type,
            unsigned int dim, unsigned int pOrder >
  class LagrangeBaseFunction
  : public BaseFunctionInterface< FunctionSpace >
  {
    typedef LagrangeBaseFunction< FunctionSpace, type, dim, pOrder > ThisType;
    typedef BaseFunctionInterface< FunctionSpace > BaseType;

  public:
    typedef typename GeometryWrapper< type, dim > :: GenericGeometryType
      GenericGeometryType;
    typedef GenericLagrangeBaseFunction
      < FunctionSpace, GenericGeometryType, pOrder >
      GenericBaseFunctionType;
      
    typedef typename GenericBaseFunctionType :: DomainType DomainType;
    typedef typename GenericBaseFunctionType :: RangeType RangeType;

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
    LagrangeBaseFunctionFactory ( GeometryType geometry );

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

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
    LagrangeBaseFunctionFactory ( GeometryType geometry );

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };
  
}

#include "basefunctions_inline.hh"

#endif

#ifndef DUNE_LAGRANGESPACE_BASEFUNCTIONS_INLINE_HH
#define DUNE_LAGRANGESPACE_BASEFUNCTIONS_INLINE_HH

#include "basefunctions.hh"

namespace Dune
{

  // LagrangeBaseFunction
  // --------------------
  
  template< class FunctionSpace, GeometryType :: BasicType type,
            unsigned int dim, unsigned int pOrder >
  LagrangeBaseFunction< FunctionSpace, type, dim, pOrder >
    :: LagrangeBaseFunction ( unsigned int baseNum )
  : BaseType(),
    baseFunction_( baseNum )
  {}
  
  template< class FunctionSpace, GeometryType :: BasicType type,
            unsigned int dim, unsigned int pOrder >
  void LagrangeBaseFunction< FunctionSpace, type, dim, pOrder >
    :: evaluate ( const FieldVector< deriType, 0 > &diffVariable,
                  const DomainType &x,
                  RangeType &phi ) const
  {
    return baseFunction_.evaluate( diffVariable, x, phi );
  }
  
  
  template< class FunctionSpace, GeometryType :: BasicType type,
            unsigned int dim, unsigned int pOrder >
  void LagrangeBaseFunction< FunctionSpace, type, dim, pOrder >
    :: evaluate ( const FieldVector< deriType, 1 > &diffVariable,
                  const DomainType &x,
                  RangeType &phi ) const
  {
    return baseFunction_.evaluate( diffVariable, x, phi );
  }
  
  
  template< class FunctionSpace, GeometryType :: BasicType type,
            unsigned int dim, unsigned int pOrder >
  void LagrangeBaseFunction< FunctionSpace, type, dim, pOrder >
    :: evaluate ( const FieldVector< deriType, 2 > &diffVariable,
                  const DomainType &x,
                  RangeType &phi ) const
  {
    return baseFunction_.evaluate( diffVariable, x, phi );
  }

  
  template< class FunctionSpace, GeometryType :: BasicType type,
            unsigned int dim, unsigned int pOrder >
  int LagrangeBaseFunction< FunctionSpace, type, dim, pOrder >
    :: order () const
  {
    return pOrder;
  }



  // LagrangeBaseFunctionFactory
  // ---------------------------

  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
    :: LagrangeBaseFunctionFactory ( GeometryType geometry )
  : BaseType( geometry )
  {
    assert( this->geometry().dim() == dim );
  }


  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
    :: ~LagrangeBaseFunctionFactory ()
  {}


  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  BaseFunctionInterface< ScalarFunctionSpace > *
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
   :: baseFunction ( int i ) const
  {
    const GeometryType :: BasicType basicType = this->geometry().basicType();

    switch( basicType )
    {
    case GeometryType :: simplex:
      return new LagrangeBaseFunction
        < ScalarFunctionSpace, GeometryType :: simplex, dim, pOrder >( i );

    case GeometryType :: cube:
      return new LagrangeBaseFunction
        < ScalarFunctionSpace, GeometryType :: cube, dim, pOrder >( i );

    default:
      DUNE_THROW( NotImplemented, "No such geometry type implemented." );
    }
  }
 

  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  int LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
    :: numBaseFunctions () const
  {
    const GeometryType :: BasicType basicType = this->geometry().basicType();

    switch( basicType )
    {
    case GeometryType :: simplex:
      return LagrangeBaseFunction
               < ScalarFunctionSpace, GeometryType :: simplex, dim, pOrder >
               :: GenericBaseFunctionType :: numBaseFunctions;

    case GeometryType :: cube:
      return LagrangeBaseFunction
               < ScalarFunctionSpace, GeometryType :: cube, dim, pOrder >
               :: GenericBaseFunctionType :: numBaseFunctions;

    default:
      DUNE_THROW( NotImplemented, "No such geometry type implemented." );
    }
  }



  template< class ScalarFunctionSpace, unsigned int pOrder >
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, 3, pOrder >
    :: LagrangeBaseFunctionFactory ( GeometryType geometry )
  : BaseType( geometry )
  {
    assert( this->geometry().dim() == 3 );
  }


  template< class ScalarFunctionSpace, unsigned int pOrder >
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, 3, pOrder >
    :: ~LagrangeBaseFunctionFactory ()
  {}

  
  template< class ScalarFunctionSpace, unsigned int pOrder >
  BaseFunctionInterface< ScalarFunctionSpace > *
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, 3, pOrder >
    :: baseFunction ( int i ) const
  {
    const GeometryType :: BasicType basicType = this->geometry().basicType();

    switch( basicType )
    {
    case GeometryType :: simplex:
      return new LagrangeBaseFunction
        < ScalarFunctionSpace, GeometryType :: simplex, 3, pOrder >( i );

    case GeometryType :: cube:
      return new LagrangeBaseFunction
        < ScalarFunctionSpace, GeometryType :: cube, 3, pOrder >( i );

    case GeometryType :: pyramid:
      return new LagrangeBaseFunction
        < ScalarFunctionSpace, GeometryType :: pyramid, 3, pOrder >( i );

    case GeometryType :: prism:
      return new LagrangeBaseFunction
        < ScalarFunctionSpace, GeometryType :: prism, 3, pOrder >( i );

    default:
      DUNE_THROW( NotImplemented, "No such geometry type implemented." );
    }
  }

  
  template< class ScalarFunctionSpace, unsigned int pOrder >
  int LagrangeBaseFunctionFactory< ScalarFunctionSpace, 3, pOrder >
    :: numBaseFunctions () const
  {
    const GeometryType :: BasicType basicType = this->geometry().basicType();

    switch( basicType )
    {
    case GeometryType :: simplex:
      return LagrangeBaseFunction
               < ScalarFunctionSpace, GeometryType :: simplex, 3, pOrder >
               :: GenericBaseFunctionType :: numBaseFunctions;

    case GeometryType :: cube:
      return LagrangeBaseFunction
               < ScalarFunctionSpace, GeometryType :: cube, 3, pOrder >
               :: GenericBaseFunctionType :: numBaseFunctions;

    case GeometryType :: pyramid:
      return LagrangeBaseFunction
               < ScalarFunctionSpace, GeometryType :: pyramid, 3, pOrder >
               :: GenericBaseFunctionType :: numBaseFunctions;

    case GeometryType :: prism:
      return LagrangeBaseFunction
               < ScalarFunctionSpace, GeometryType :: prism, 3, pOrder >
               :: GenericBaseFunctionType :: numBaseFunctions;

    default:
      DUNE_THROW( NotImplemented, "No such geometry type implemented." );
    }
  }

}

#endif

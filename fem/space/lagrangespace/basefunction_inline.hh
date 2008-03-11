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
    return polynomialOrder;
  }



  // LagrangeBaseFunctionFactory
  // ---------------------------

  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
    :: LagrangeBaseFunctionFactory ( GeometryType geometry )
  : BaseType( geometry )
  {
    assert( this->geometry().dim() == dimension );
  }


  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
    :: ~LagrangeBaseFunctionFactory ()
  {}


  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  typename LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
    :: BaseFunctionType *
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
   :: baseFunction ( int i ) const
  {
    const GeometryType :: BasicType basicType = this->geometry().basicType();

    switch( basicType )
    {
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
 

  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  int LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
    :: numBaseFunctions () const
  {
    const GeometryType :: BasicType basicType = this->geometry().basicType();

    switch( basicType )
    {
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



  template< class ScalarFunctionSpace, unsigned int pOrder >
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, 3, pOrder >
    :: LagrangeBaseFunctionFactory ( GeometryType geometry )
  : BaseType( geometry )
  {
    assert( this->geometry().dim() == dimension );
  }


  template< class ScalarFunctionSpace, unsigned int pOrder >
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, 3, pOrder >
    :: ~LagrangeBaseFunctionFactory ()
  {}

  
  template< class ScalarFunctionSpace, unsigned int pOrder >
  typename LagrangeBaseFunctionFactory< ScalarFunctionSpace, 3, pOrder >
    :: BaseFunctionType *
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, 3, pOrder >
    :: baseFunction ( int i ) const
  {
    const GeometryType :: BasicType basicType = this->geometry().basicType();

    switch( basicType )
    {
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

  
  template< class ScalarFunctionSpace, unsigned int pOrder >
  int LagrangeBaseFunctionFactory< ScalarFunctionSpace, 3, pOrder >
    :: numBaseFunctions () const
  {
    const GeometryType :: BasicType basicType = this->geometry().basicType();

    switch( basicType )
    {
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

}

#ifndef DUNE_FEM_SPACE_PADAPTIVE_LAGRANGEBASEFUNCTIONS_HH
#define DUNE_FEM_SPACE_PADAPTIVE_LAGRANGEBASEFUNCTIONS_HH

#include <dune/fem/space/lagrange/genericbasefunctions.hh>


namespace Dune
{

  namespace Fem 
  {

    // LagrangeBaseFunction
    // --------------------
     
    template< class FunctionSpace, unsigned int topologyId,
              unsigned int dim, unsigned int pOrder >
    class LagrangeBaseFunction
    : public BaseFunctionInterface< FunctionSpace >
    {
      typedef LagrangeBaseFunction< FunctionSpace, topologyId, dim, pOrder > ThisType;
      typedef BaseFunctionInterface< FunctionSpace > BaseType;

    public:
      typedef typename GeometryWrapper< topologyId, dim >::GenericGeometryType
        GenericGeometryType;
      typedef GenericLagrangeBaseFunction< FunctionSpace, GenericGeometryType, pOrder >
        GenericBaseFunctionType;

      enum { numBaseFunctions = GenericBaseFunctionType::numBaseFunctions };

      typedef typename GenericBaseFunctionType :: DomainType DomainType;
      typedef typename GenericBaseFunctionType :: RangeType RangeType;

      LagrangeBaseFunction( unsigned int baseNum );

      virtual void evaluate ( const FieldVector< int, 0 > &diffVariable,
                              const DomainType &x,
                              RangeType &phi ) const;
      
      virtual void evaluate ( const FieldVector< int, 1 > &diffVariable,
                              const DomainType &x,
                              RangeType &phi ) const;

      virtual void evaluate ( const FieldVector< int, 2 > &diffVariable,
                              const DomainType &x,
                              RangeType &phi ) const;

      virtual int order () const;

    protected:
      GenericBaseFunctionType baseFunction_;
    };



    // LagrangeBaseFunctionFactory
    // ---------------------------

    template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
    class LagrangeBaseFunctionFactory
    : public BaseFunctionFactory< ScalarFunctionSpace >
    {
      typedef LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
        ThisType;
      typedef BaseFunctionFactory< ScalarFunctionSpace > BaseType;

      template< class Topology >
      struct Switcher;

    public:
      using BaseType::geometry;

      explicit LagrangeBaseFunctionFactory ( GeometryType geometry );

      virtual ~LagrangeBaseFunctionFactory ();

      virtual BaseFunctionInterface< ScalarFunctionSpace > *
      baseFunction ( int i ) const;

      virtual int numBaseFunctions () const;
    };

  } // namespace Fem 

} // namespace Dune

#include "lagrangebasefunctions_inline.hh"

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVE_LAGRANGEBASEFUNCTIONS_HH

#ifndef DUNE_LAGRANGESPACE_BASEFUNCTIONS_INLINE_HH
#define DUNE_LAGRANGESPACE_BASEFUNCTIONS_INLINE_HH

#include <dune/common/version.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include "basefunctions.hh"

namespace Dune
{

  // LagrangeBaseFunction
  // --------------------
  
  template< class FunctionSpace, unsigned int topologyId,
            unsigned int dim, unsigned int pOrder >
  LagrangeBaseFunction< FunctionSpace, topologyId, dim, pOrder >
    :: LagrangeBaseFunction ( unsigned int baseNum )
  : BaseType(),
    baseFunction_( baseNum )
  {}
  
  template< class FunctionSpace, unsigned int topologyId,
            unsigned int dim, unsigned int pOrder >
  void LagrangeBaseFunction< FunctionSpace, topologyId, dim, pOrder >
    :: evaluate ( const FieldVector< int, 0 > &diffVariable,
                  const DomainType &x,
                  RangeType &phi ) const
  {
    return baseFunction_.evaluate( diffVariable, x, phi );
  }
  
  
  template< class FunctionSpace, unsigned int topologyId,
            unsigned int dim, unsigned int pOrder >
  void LagrangeBaseFunction< FunctionSpace, topologyId, dim, pOrder >
    :: evaluate ( const FieldVector< int, 1 > &diffVariable,
                  const DomainType &x,
                  RangeType &phi ) const
  {
    return baseFunction_.evaluate( diffVariable, x, phi );
  }
  
  
  template< class FunctionSpace, unsigned int topologyId,
            unsigned int dim, unsigned int pOrder >
  void LagrangeBaseFunction< FunctionSpace, topologyId, dim, pOrder >
    :: evaluate ( const FieldVector< int, 2 > &diffVariable,
                  const DomainType &x,
                  RangeType &phi ) const
  {
    return baseFunction_.evaluate( diffVariable, x, phi );
  }

  
  template< class FunctionSpace, unsigned int topologyId,
            unsigned int dim, unsigned int pOrder >
  int LagrangeBaseFunction< FunctionSpace, topologyId, dim, pOrder >
    :: order () const
  {
    return pOrder;
  }



  // LagrangeBaseFunctionFactory::TopologyId
  // ---------------------------------------

  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  template< class Topology >
  struct LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >::Switcher
  {
    typedef LagrangeBaseFunction< ScalarFunctionSpace, Topology::id, dim, pOrder > BaseFunction;

    static void apply( int &numBaseFunctions )
    {
      numBaseFunctions = BaseFunction::GenericBaseFunctionType::numBaseFunctions; 
    }

    static void apply ( const int &i, BaseFunctionInterface< ScalarFunctionSpace > *&baseFunction )
    {
      baseFunction = new BaseFunction( i );
    }
  };


  // LagrangeBaseFunctionFactory
  // ---------------------------

  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
    ::LagrangeBaseFunctionFactory ( GeometryType geometry )
  : BaseType( geometry )
  {
    assert( geometry.dim() == dim );
  }


  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
    :: ~LagrangeBaseFunctionFactory ()
  {}


  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  BaseFunctionInterface< ScalarFunctionSpace > *
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >::baseFunction ( int i ) const
  {
    BaseFunctionInterface< ScalarFunctionSpace > *baseFunction;
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,1,0)
    const unsigned int topologyId = geometry().id();
#else
    const unsigned int topologyId = GenericGeometry::topologyId( geometry() );
#endif
    GenericGeometry::IfTopology< Switcher, dim >::apply( topologyId, i, baseFunction );
    return baseFunction;
  }
 

  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  int
  LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >::numBaseFunctions () const
  {
    int numBaseFunctions;
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,1,0)
    const unsigned int topologyId = geometry().id();
#else
    const unsigned int topologyId = GenericGeometry::topologyId( geometry() );
#endif
    GenericGeometry::IfTopology< Switcher, dim >::apply( topologyId, numBaseFunctions );
    return numBaseFunctions;
  }

}

#endif // #ifndef DUNE_LAGRANGESPACE_BASEFUNCTIONS_INLINE_HH

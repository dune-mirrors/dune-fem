#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LAGRANGE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LAGRANGE_HH

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/lagrange/genericbasefunctions.hh>
#include <dune/fem/space/lagrange/shapefunctionset.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>

// local includes
#include "default.hh"


namespace Dune
{

  namespace Fem
  {

    // LagrangeDiscontinuousGalerkinSpaceTraits
    // ----------------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    struct LagrangeDiscontinuousGalerkinSpaceTraits
    {
      typedef LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = 0;
      static const int polynomialOrder = polOrder;

    private:
      static const int dimLocal = GridPartType::dimension;
      typedef typename ToLocalFunctionSpace< FunctionSpaceType, dimLocal >::Type ShapeFunctionSpaceType;

      template< class GridPartType, bool hasSingleGeometryType = Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPartType >::v >
      struct HaveSingleGeometryType;

      template< class GridPartType >
      class HaveSingleGeometryType< GridPartType, true >
      {
        typedef typename GeometryWrapper<
            Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPartType >::topologyId, dimLocal
          >::GenericGeometryType GenericGeometryType;
        typedef GenericLagrangeBaseFunction<
            typename FunctionSpaceType::ScalarFunctionSpaceType, GenericGeometryType, polOrder
          > GenericBaseFunctionType;

      public:
        static const int localBlockSize = FunctionSpaceType::dimRange * GenericBaseFunctionType::numBaseFunctions;
      };

    public:
      static const int localBlockSize = HaveSingleGeometryType< GridPartType >::localBlockSize;

      typedef LagrangeShapeFunctionSet< ShapeFunctionSpaceType, polOrder > LagrangeShapeFunctionSetType;
      typedef SelectCachingShapeFunctionSet< LagrangeShapeFunctionSetType, Storage > ScalarShapeFunctionSetType;

      struct ScalarShapeFunctionSetFactory
      {
        static ScalarShapeFunctionSetType *createObject ( const GeometryType &type )
        {
          return new ScalarShapeFunctionSetType( type, LagrangeShapeFunctionSetType( type ) );
        }

        static void deleteObject ( ScalarShapeFunctionSetType *object ) { delete object; }
      };
      typedef ScalarShapeFunctionSetFactory ScalarShapeFunctionSetFactoryType;
    };



    // LagrangeDiscontinuousGalerkinSpace
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    class LagrangeDiscontinuousGalerkinSpace
    : public DiscontinuousGalerkinSpaceDefault< 
        DiscontinuousGalerkinSpaceDefaultTraits< LagrangeDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >, Storage 
      >
    {
      typedef LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscontinuousGalerkinSpaceDefault< 
          DiscontinuousGalerkinSpaceDefaultTraits< LagrangeDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >, Storage 
        > BaseType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      LagrangeDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                           const InterfaceType commInterface = BaseType::defaultInterface,
                                           const CommunicationDirection commDirection = BaseType::defaultDirection )
      : BaseType( gridPart, commInterface, commDirection )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LAGRANGE_HH

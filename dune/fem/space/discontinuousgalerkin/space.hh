#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-common includes
#include <dune/common/static_assert.hh>

// dune-fem includes
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>

// local includes
#include <dune/fem/space/discontinuousgalerkin/capabilities.hh>
#include <dune/fem/space/discontinuousgalerkin/default.hh>
#include <dune/fem/space/discontinuousgalerkin/shapefunctionset.hh>


namespace Dune
{

  namespace Fem
  {

    // DiscontinuousGalerkinSpaceTraits
    // --------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    struct DiscontinuousGalerkinSpaceTraits
    {
      dune_static_assert( (GridPart::dimensionworld <= 3), "Use Legendre spaces for higher spatial dimensions." );

      typedef DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = 0;
      static const int polynomialOrder = polOrder;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      static const int dimLocal = GridPartType::dimension;
      typedef typename ToLocalFunctionSpace< ScalarFunctionSpaceType, dimLocal >::Type ScalarShapeFunctionSpaceType;

    public:
      static const int localBlockSize = FunctionSpaceType::dimRange * OrthonormalShapeFunctionSetSize< ScalarShapeFunctionSpaceType, polOrder >::v;

      typedef OrthonormalShapeFunctionSet< ScalarShapeFunctionSpaceType, polOrder > OrthonormalShapeFunctionSetType;
      typedef SelectCachingShapeFunctionSet< OrthonormalShapeFunctionSetType, Storage > ScalarShapeFunctionSetType;

      struct ScalarShapeFunctionSetFactory
      {
        static ScalarShapeFunctionSetType *createObject ( const GeometryType &type )
        {
          return new ScalarShapeFunctionSetType( type, OrthonormalShapeFunctionSetType( type ) );
        }

        static void deleteObject ( ScalarShapeFunctionSetType *object ) { delete object; }
      };
      typedef ScalarShapeFunctionSetFactory ScalarShapeFunctionSetFactoryType;
    };



    // DiscontinuousGalerkinSpace
    // --------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    class DiscontinuousGalerkinSpace
    : public DiscontinuousGalerkinSpaceDefault< 
        DiscontinuousGalerkinSpaceDefaultTraits< DiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
      >
    {
      typedef DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType; 
      typedef DiscontinuousGalerkinSpaceDefault< 
          DiscontinuousGalerkinSpaceDefaultTraits< DiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
        > BaseType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      DiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                   const InterfaceType commInterface = BaseType::defaultInterface,
                                   const CommunicationDirection commDirection = BaseType::defaultDirection )
      : BaseType( gridPart, commInterface, commDirection )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH

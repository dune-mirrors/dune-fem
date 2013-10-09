#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HLEGENDRE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HLEGENDRE_HH

// dune-common includes
#include <dune/common/power.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/legendre.hh>
#include <dune/fem/space/shapefunctionset/mapped.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>

// local includes
#include "default.hh"
#include "hierarchiclegendremap.hh"

namespace Dune
{

  namespace Fem
  {

    // LegendreDiscontinuousGalerkinSpaceTraits
    // ----------------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    struct HierarchicLegendreDiscontinuousGalerkinSpaceTraits
    {
      typedef  HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;
      
      static const int codimension = 0;
      static const int polynomialOrder = polOrder;

    private:
      template <int p, int dim>
      struct NumLegendreShapeFunctions
      {
        static const int v = StaticPower< p+1, dim >::power;
      };

      static const int dimLocal = GridPartType::dimension;

    public:
      static const size_t numShapeFunctions = NumLegendreShapeFunctions< polOrder, dimLocal >::v;

      static const int localBlockSize = FunctionSpaceType::dimRange * NumLegendreShapeFunctions< polOrder, dimLocal >::v;

      typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      typedef LegendreShapeFunctionSet< typename ToNewDimDomainFunctionSpace< ScalarFunctionSpaceType, dimLocal >::Type > OriginalShapeFunctionSetType;
      
      typedef HierarchicLegendreMap<polOrder,dimLocal> MappingType;
     
      typedef MappedShapeFunctionSet<OriginalShapeFunctionSetType, MappingType> LegendreShapeFunctionSetType;
  

      typedef SelectCachingShapeFunctionSet< LegendreShapeFunctionSetType, Storage > ScalarShapeFunctionSetType;
     
      struct ScalarShapeFunctionSetFactory
      {
         
        static ScalarShapeFunctionSetType *createObject ( const GeometryType &type )
        {
          return new ScalarShapeFunctionSetType( type, 
                                                 LegendreShapeFunctionSetType( OriginalShapeFunctionSetType( polOrder ) ) );
        }

        static void deleteObject ( ScalarShapeFunctionSetType *object ) { delete object; }
      
       
     
      };
      typedef ScalarShapeFunctionSetFactory ScalarShapeFunctionSetFactoryType;
    };



    // LegendreDiscontinuousGalerkinSpace
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    class HierarchicLegendreDiscontinuousGalerkinSpace
    : public DiscontinuousGalerkinSpaceDefault< 
        DiscontinuousGalerkinSpaceDefaultTraits< HierarchicLegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
      >
    {
      typedef HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscontinuousGalerkinSpaceDefault< 
          DiscontinuousGalerkinSpaceDefaultTraits< HierarchicLegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
        > BaseType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      HierarchicLegendreDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                           const InterfaceType commInterface = BaseType::defaultInterface,
                                           const CommunicationDirection commDirection = BaseType::defaultDirection )
      : BaseType( gridPart, commInterface, commDirection )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HLEGENDRE_HH

#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HLEGENDRE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HLEGENDRE_HH

#include <cassert>

#include <dune/geometry/type.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/legendre.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>

#include "legendre.hh"

namespace Dune
{

  namespace Fem
  {

    // HierarchicalLegendreDiscontinuousGalerkinSpace
    // ----------------------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage >
    class HierarchicLegendreDiscontinuousGalerkinSpace
    : public LegendreDiscontinuousGalerkinSpaceBase< FunctionSpace, GridPart, polOrder, Storage, true >
    {
      // hierarchicalOrdering = true
      typedef LegendreDiscontinuousGalerkinSpaceBase< FunctionSpace, GridPart, polOrder, Storage, true > BaseType;
      typedef HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;

    public:
      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType   EntityType;
      typedef DiscontinuousGalerkinLocalInterpolation< ThisType > InterpolationType;
      typedef InterpolationType  InterpolationImplType;

      explicit HierarchicLegendreDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                              const InterfaceType commInterface = InteriorBorder_All_Interface,
                                                              const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection )
      {}

      InterpolationType interpolation () const
      {
        return InterpolationType( *this );
      }

      [[deprecated("Use LocalInterpolation( space ) instead!")]]
      InterpolationType interpolation ( const EntityType &entity ) const
      {
        return interpolation();
      }

      InterpolationType localInterpolation () const
      {
        return interpolation();
      }


    };

    namespace Capabilities
    {

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasFixedPolynomialOrder< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasStaticPolynomialOrder< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isContinuous< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isLocalized< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isAdaptive< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct threadSafe< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct viewThreadSafe< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isHierarchic< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HLEGENDRE_HH

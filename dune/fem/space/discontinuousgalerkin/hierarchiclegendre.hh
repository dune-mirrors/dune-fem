#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HLEGENDRE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HLEGENDRE_HH

#include <dune/common/power.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/legendre.hh>
#include <dune/fem/space/shapefunctionset/mapped.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>

#include "basisfunctionsets.hh"
#include "generic.hh"
#include "hierarchiclegendremap.hh"
#include "interpolation.hh"
#include "shapefunctionsets.hh"

namespace Dune
{

  namespace Fem
  {

    // HierarchicLegendreDiscontinuousGalerkinSpaceTraits
    // --------------------------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    struct HierarchicLegendreDiscontinuousGalerkinSpaceTraits
    {
      typedef HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = 0;

      typedef Dune::Fem::FunctionSpace<
          typename FunctionSpace::DomainFieldType, typename FunctionSpace::RangeFieldType,
           GridPartType::dimension, 1
        > ScalarShapeFunctionSpaceType;

      typedef LegendreShapeFunctionSet< ScalarShapeFunctionSpaceType > HostShapeFunctionSetType;
      typedef HierarchicLegendreMap< polOrder, GridPartType::dimension> MappingType;
      struct ScalarShapeFunctionSet
      : public MappedShapeFunctionSet< HostShapeFunctionSetType, MappingType >
      {
        typedef MappedShapeFunctionSet< HostShapeFunctionSetType, MappingType > BaseType;

      public:
        explicit ScalarShapeFunctionSet ( Dune::GeometryType type )
          : BaseType( HostShapeFunctionSetType( polOrder ) )
        {
          assert( type.isCube() );
        }
      };

      typedef SelectCachingShapeFunctionSets< GridPartType, ScalarShapeFunctionSet, Storage > ScalarShapeFunctionSetsType;
      typedef VectorialShapeFunctionSets< ScalarShapeFunctionSetsType, typename FunctionSpaceType::RangeType > ShapeFunctionSetsType;

      typedef DefaultBasisFunctionSets< GridPartType, ShapeFunctionSetsType > BasisFunctionSetsType;
      typedef typename BasisFunctionSetsType::BasisFunctionSetType BasisFunctionSetType;

      typedef CodimensionMapper< GridPartType, codimension > BlockMapperType;
      static const int localBlockSize
        = FunctionSpaceType::dimRange*StaticPower< polOrder+1, GridPartType::dimension >::power;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef Operation OperationType;
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      };
    };



    // HierarchicLegendreDiscontinuousGalerkinSpace
    // --------------------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    class HierarchicLegendreDiscontinuousGalerkinSpace
    : public GenericDiscontinuousGalerkinSpace< HierarchicLegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef GenericDiscontinuousGalerkinSpace< HierarchicLegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

    public:
      using BaseType::basisFunctionSet;

      static const int polynomialOrder = polOrder;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::BasisFunctionSetsType BasisFunctionSetsType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef DiscontinuousGalerkinLocalL2Projection< GridPartType, BasisFunctionSetType > InterpolationType;

      explicit HierarchicLegendreDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                              const InterfaceType commInterface = InteriorBorder_All_Interface,
                                                              const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, basisFunctionSets( gridPart ), commInterface, commDirection )
      {}

      static DFSpaceIdentifier type () { return HierarchicLegendreDGSpace_id; }

      InterpolationType interpolation ( const EntityType &entity ) const
      {
        return InterpolationType( basisFunctionSet( entity ) );
      }

    private:
      static BasisFunctionSetsType basisFunctionSets ( const GridPartType &gridPart )
      {
        typedef typename BasisFunctionSetsType::ShapeFunctionSetsType ShapeFunctionSetsType;
        ShapeFunctionSetsType shapeFunctionSets( gridPart );
        return BasisFunctionSetsType( std::move( shapeFunctionSets ) );
      }
    };



    namespace Capabilities
    {

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct hasFixedPolynomialOrder< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct hasStaticPolynomialOrder< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isContinuous< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isLocalized< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isParallel< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::isParallel< GridPart >::v;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isAdaptive< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct threadSafe< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct viewThreadSafe< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HLEGENDRE_HH

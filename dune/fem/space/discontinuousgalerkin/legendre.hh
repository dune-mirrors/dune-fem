#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH

#include <cassert>

#include <dune/common/math.hh>

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

#include "basisfunctionsets.hh"
#include "declaration.hh"
#include "generic.hh"
#include "localinterpolation.hh"
#include "shapefunctionsets.hh"

namespace Dune
{

  namespace Fem
  {
    // Forward declaration
    template< class FunctionSpace, class GridPart, int polOrder, class Storage = CachingStorage >
    class LegendreDiscontinuousGalerkinSpace;

    template< class FunctionSpace, class GridPart, int polOrder, class Storage = CachingStorage >
    class HierarchicLegendreDiscontinuousGalerkinSpace;

    // LegendreDiscontinuousGalerkinSpaceTraits
    // ----------------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage, bool hierarchicalOrdering  >
    struct LegendreDiscontinuousGalerkinSpaceTraits
    {
      // select space implementation depending on basis function ordering
      typedef typename std::conditional< hierarchicalOrdering,
          HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage >,
          LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >::type  DiscreteFunctionSpaceType;

      typedef GridPart GridPartType;
      typedef GridFunctionSpace< GridPartType, FunctionSpace > FunctionSpaceType;

      static const int codimension = 0;

      typedef Dune::Fem::FunctionSpace<
          typename FunctionSpace::DomainFieldType, typename FunctionSpace::RangeFieldType,
           GridPartType::dimension, 1
        > ScalarShapeFunctionSpaceType;

      struct ScalarShapeFunctionSet
        : public Dune::Fem::LegendreShapeFunctionSet< ScalarShapeFunctionSpaceType, hierarchicalOrdering >
      {
        typedef Dune::Fem::LegendreShapeFunctionSet< ScalarShapeFunctionSpaceType, hierarchicalOrdering > BaseType;
        static const int numberShapeFunctions = Dune::power( int(polOrder+1), int(ScalarShapeFunctionSpaceType::dimDomain) );

      public:
        explicit ScalarShapeFunctionSet ( Dune::GeometryType type )
          : BaseType( polOrder )
        {
          assert( type.isCube() );
          assert( size() == BaseType::size() );
        }

        // overload size method because it's a static value
        static constexpr unsigned int size() { return numberShapeFunctions; }
      };

      typedef SelectCachingShapeFunctionSets< GridPartType, ScalarShapeFunctionSet, Storage > ScalarShapeFunctionSetsType;
      typedef VectorialShapeFunctionSets< ScalarShapeFunctionSetsType, typename FunctionSpaceType::RangeType > ShapeFunctionSetsType;

      typedef DefaultBasisFunctionSets< GridPartType, ShapeFunctionSetsType > BasisFunctionSetsType;
      typedef typename BasisFunctionSetsType::BasisFunctionSetType BasisFunctionSetType;

      typedef CodimensionMapper< GridPartType, codimension > BlockMapperType;

      typedef Hybrid::IndexRange< int, FunctionSpaceType::dimRange * ScalarShapeFunctionSet::numberShapeFunctions > LocalBlockIndices;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef Operation OperationType;
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      };
    };


    // LegendreDiscontinuousGalerkinSpaceBase
    // --------------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage, bool hierarchicalOrdering >
    class LegendreDiscontinuousGalerkinSpaceBase
    : public GenericDiscontinuousGalerkinSpace< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage, hierarchicalOrdering > >
    {
      typedef GenericDiscontinuousGalerkinSpace< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage, hierarchicalOrdering > > BaseType;
      typedef LegendreDiscontinuousGalerkinSpaceBase< FunctionSpace, GridPart, polOrder, Storage, hierarchicalOrdering > ThisType;

    public:
      using BaseType::basisFunctionSet;

      static const int polynomialOrder = polOrder;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::BasisFunctionSetsType BasisFunctionSetsType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      explicit LegendreDiscontinuousGalerkinSpaceBase ( GridPartType &gridPart,
                                                        const InterfaceType commInterface = InteriorBorder_All_Interface,
                                                        const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, makeBasisFunctionSets( gridPart ), commInterface, commDirection )
      {}

      static DFSpaceIdentifier type () { return LegendreDGSpace_id; }

    private:
      static BasisFunctionSetsType makeBasisFunctionSets ( const GridPartType &gridPart )
      {
        typedef typename BasisFunctionSetsType::ShapeFunctionSetsType ShapeFunctionSetsType;
        ShapeFunctionSetsType shapeFunctionSets( gridPart );
        return BasisFunctionSetsType( std::move( shapeFunctionSets ) );
      }
    };

    // LegendreDiscontinuousGalerkinSpace
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage >
    class LegendreDiscontinuousGalerkinSpace
    : public LegendreDiscontinuousGalerkinSpaceBase< FunctionSpace, GridPart, polOrder, Storage, false >
    {
      // hierarchicalOrdering = false
      typedef LegendreDiscontinuousGalerkinSpaceBase< FunctionSpace, GridPart, polOrder, Storage, false > BaseType;
      typedef LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;

    public:
      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType   EntityType;
      typedef DiscontinuousGalerkinLocalInterpolation< ThisType > InterpolationType;
      typedef InterpolationType  InterpolationImplType;

      explicit LegendreDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                    const InterfaceType commInterface = InteriorBorder_All_Interface,
                                                    const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection )
      {}

      InterpolationType interpolation () const
      {
        return InterpolationType( *this );
      }

      [[deprecated]]
      InterpolationType interpolation ( const EntityType &entity ) const
      {
        return interpolation();
      }

      InterpolationType localInterpolation ( const EntityType &entity ) const
      {
        return interpolation();
      }
    };


    namespace Capabilities
    {

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasFixedPolynomialOrder< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasStaticPolynomialOrder< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isContinuous< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isLocalized< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isAdaptive< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct threadSafe< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct viewThreadSafe< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH

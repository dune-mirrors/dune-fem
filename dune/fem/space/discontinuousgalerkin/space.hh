#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH

#include <utility>

#include <dune/geometry/type.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/orthonormal.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>

#include "basisfunctionsets.hh"
#include "declaration.hh"
#include "generic.hh"
#include "interpolation.hh"
#include "shapefunctionsets.hh"

namespace Dune
{

  namespace Fem
  {

    // DiscontinuousGalerkinSpaceTraits
    // --------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    struct DiscontinuousGalerkinSpaceTraits
    {
      typedef DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = 0;

      typedef Dune::Fem::FunctionSpace<
          typename FunctionSpace::DomainFieldType, typename FunctionSpace::RangeFieldType,
           GridPartType::dimension, 1
        > ScalarShapeFunctionSpaceType;

      typedef SelectCachingShapeFunctionSets< GridPartType, OrthonormalShapeFunctionSet< ScalarShapeFunctionSpaceType, polOrder >, Storage > ScalarShapeFunctionSetsType;
      typedef VectorialShapeFunctionSets< ScalarShapeFunctionSetsType, typename FunctionSpaceType::RangeType > ShapeFunctionSetsType;

      typedef DefaultBasisFunctionSets< GridPartType, ShapeFunctionSetsType > BasisFunctionSetsType;
      typedef typename BasisFunctionSetsType::BasisFunctionSetType BasisFunctionSetType;

      typedef CodimensionMapper< GridPartType, codimension > BlockMapperType;
      static const int localBlockSize
        = FunctionSpaceType::dimRange*OrthonormalShapeFunctionSetSize< ScalarShapeFunctionSpaceType, polOrder >::v;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef Operation OperationType;
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      };
    };



    // DiscontinuousGalerkinSpace
    // --------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    class DiscontinuousGalerkinSpace
    : public GenericDiscontinuousGalerkinSpace< DiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef GenericDiscontinuousGalerkinSpace< DiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

    public:
      using BaseType::basisFunctionSet;

      static const int polynomialOrder = polOrder;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::BasisFunctionSetsType BasisFunctionSetsType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef DiscontinuousGalerkinLocalL2Projection< GridPartType, BasisFunctionSetType > InterpolationType;

      explicit DiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                            const InterfaceType commInterface = InteriorBorder_All_Interface,
                                            const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, basisFunctionSets( gridPart ), commInterface, commDirection )
      {}

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
      struct hasFixedPolynomialOrder< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct hasStaticPolynomialOrder< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isContinuous< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isLocalized< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isAdaptive< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct threadSafe< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct viewThreadSafe< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isHierarchic< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH

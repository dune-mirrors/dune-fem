#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LAGRANGE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LAGRANGE_HH

#include <dune/geometry/type.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/common/hybrid.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/operator/projection/local/l2projection.hh>
#include <dune/fem/operator/projection/local/riesz.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/localinterpolation.hh>
#include <dune/fem/space/lagrange/genericbasefunctions.hh>
#include <dune/fem/space/lagrange/shapefunctionset.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>

#include "basisfunctionsets.hh"
#include "declaration.hh"
#include "generic.hh"
#include "shapefunctionsets.hh"
#include "localinterpolation.hh"

namespace Dune
{

  namespace Fem
  {

    // LagrangeDiscontinuousGalerkinSpaceTraits
    // ----------------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage >
    struct LagrangeDiscontinuousGalerkinSpaceTraits
    {
      typedef LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      typedef GridPart GridPartType;
      typedef GridFunctionSpace< GridPartType, FunctionSpace > FunctionSpaceType;

      static const int codimension = 0;

      typedef Dune::Fem::FunctionSpace<
          typename FunctionSpace::DomainFieldType, typename FunctionSpace::RangeFieldType,
           GridPartType::dimension, 1
        > ScalarShapeFunctionSpaceType;

      typedef SelectCachingShapeFunctionSets< GridPartType, LagrangeShapeFunctionSet< ScalarShapeFunctionSpaceType, polOrder >, Storage > ScalarShapeFunctionSetsType;
      typedef VectorialShapeFunctionSets< ScalarShapeFunctionSetsType, typename FunctionSpaceType::RangeType > ShapeFunctionSetsType;

      typedef DefaultBasisFunctionSets< GridPartType, ShapeFunctionSetsType > BasisFunctionSetsType;
      typedef typename BasisFunctionSetsType::BasisFunctionSetType BasisFunctionSetType;

      typedef CodimensionMapper< GridPartType, codimension > BlockMapperType;

      typedef typename GeometryWrapper<
          Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPartType >::topologyId, GridPartType::dimension
        >::ImplType ImplType;
      typedef GenericLagrangeBaseFunction<
          typename FunctionSpaceType::ScalarFunctionSpaceType, ImplType, polOrder
        > GenericBaseFunctionType;

      typedef Hybrid::IndexRange< int, FunctionSpaceType::dimRange * GenericBaseFunctionType::numBaseFunctions > LocalBlockIndices;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef Operation OperationType;
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      };
    };



    // LagrangeDiscontinuousGalerkinSpace
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage = CachingStorage >
    class LagrangeDiscontinuousGalerkinSpace
    : public GenericDiscontinuousGalerkinSpace< LagrangeDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef GenericDiscontinuousGalerkinSpace< LagrangeDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;
      typedef LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;

    public:
      using BaseType::basisFunctionSet;

      static const int polynomialOrder = polOrder;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::BasisFunctionSetsType BasisFunctionSetsType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

    private:
      typedef CachingQuadrature< GridPart, EntityType::codimension > QuadratureType;
      typedef DenseLocalRieszProjection< BasisFunctionSetType, QuadratureType > LocalRieszProjectionType;

    public:
      typedef DefaultLocalL2Projection< LocalRieszProjectionType, QuadratureType > InterpolationImplType;
      typedef LocalInterpolationWrapper< ThisType > InterpolationType;

      explicit LagrangeDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                    const InterfaceType commInterface = InteriorBorder_All_Interface,
                                                    const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, makeBasisFunctionSets( gridPart ), commInterface, commDirection )
      {}

      static DFSpaceIdentifier type () { return LagrangeDGSpace_id; }

      InterpolationType interpolation () const
      {
        return InterpolationType( *this );
      }

      [[deprecated]]
      InterpolationImplType interpolation ( const EntityType &entity ) const
      {
        return InterpolationImplType( basisFunctionSet( entity ) );
      }

      InterpolationImplType localInterpolation ( const EntityType &entity ) const
      {
        return InterpolationImplType( basisFunctionSet( entity ) );
      }


      //const InterpolationType& interpolation ( const EntityType &entity ) const
      //{
      //  return interpolation_;
      //}

    private:
      //InterpolationType interpolation_;

      static BasisFunctionSetsType makeBasisFunctionSets ( const GridPartType &gridPart )
      {
        typedef typename BasisFunctionSetsType::ShapeFunctionSetsType ShapeFunctionSetsType;
        ShapeFunctionSetsType shapeFunctionSets( gridPart );
        return BasisFunctionSetsType( std::move( shapeFunctionSets ) );
      }
    };



    namespace Capabilities
    {

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasFixedPolynomialOrder< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasStaticPolynomialOrder< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isContinuous< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isLocalized< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isAdaptive< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct threadSafe< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct viewThreadSafe< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LAGRANGE_HH

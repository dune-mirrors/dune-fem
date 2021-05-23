#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH

#include <utility>

#include <dune/geometry/type.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/common/hybrid.hh>
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
#include "localinterpolation.hh"
#include "interpolation.hh"
#include "shapefunctionsets.hh"

namespace Dune
{

  namespace Fem
  {

    // DiscontinuousGalerkinSpaceTraits
    // --------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage >
    struct DiscontinuousGalerkinSpaceTraits
    {
      typedef DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      typedef GridPart GridPartType;
      typedef GridFunctionSpace< GridPartType, FunctionSpace > FunctionSpaceType;

      static const int codimension = 0;

      typedef Dune::Fem::FunctionSpace<
          typename FunctionSpace::DomainFieldType, typename FunctionSpace::RangeFieldType,
           GridPartType::dimension, 1
        > ScalarShapeFunctionSpaceType;

      struct ScalarShapeFunctionSet
        : public Dune::Fem::OrthonormalShapeFunctionSet< ScalarShapeFunctionSpaceType >
      {
        typedef Dune::Fem::OrthonormalShapeFunctionSet< ScalarShapeFunctionSpaceType >   BaseType;

        static constexpr int numberShapeFunctions =
              OrthonormalShapeFunctions< ScalarShapeFunctionSpaceType::dimDomain >::size(polOrder);
      public:
        explicit ScalarShapeFunctionSet ( Dune::GeometryType type )
          : BaseType( type, polOrder )
        {
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



    // DiscontinuousGalerkinSpace
    // --------------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage = CachingStorage >
    class DiscontinuousGalerkinSpace
    : public GenericDiscontinuousGalerkinSpace< DiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef GenericDiscontinuousGalerkinSpace< DiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;
      typedef DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;

    public:
      using BaseType::basisFunctionSet;

      static const int polynomialOrder = polOrder;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::BasisFunctionSetsType BasisFunctionSetsType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef DiscontinuousGalerkinLocalInterpolation< ThisType > InterpolationType;
      typedef InterpolationType  InterpolationImplType;


      explicit DiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                            const InterfaceType commInterface = InteriorBorder_All_Interface,
                                            const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, makeBasisFunctionSets( gridPart ), commInterface, commDirection )
      {}

      InterpolationType interpolation () const
      {
        return InterpolationType( *this );
      }

      [[deprecated]]
      InterpolationType interpolation ( const EntityType &entity ) const
      {
        return interpolation ();
      }

      InterpolationType localInterpolation ( const EntityType &entity ) const
      {
        return interpolation ();
      }

    private:
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
      struct hasFixedPolynomialOrder< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasStaticPolynomialOrder< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isContinuous< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isLocalized< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isAdaptive< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct threadSafe< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct viewThreadSafe< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isHierarchic< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH

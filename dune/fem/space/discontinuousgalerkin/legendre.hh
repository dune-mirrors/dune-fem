#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH

// C++ includes
#include <cassert>

// dune-geometry includes
#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>

// dune-fem includes
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/legendre.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

// local includes
#include "declaration.hh"
#include "default.hh"
#include "legendredgbasefunctions.hh"


namespace Dune
{

  namespace Fem
  {

    // LegendreDiscontinuousGalerkinSpaceTraits
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class LegendreDiscontinuousGalerkinSpaceTraits
    : public DiscontinuousGalerkinSpaceTraitsBase< FunctionSpace, GridPart, polOrder, Storage >
    {
      typedef DiscontinuousGalerkinSpaceTraitsBase< FunctionSpace, GridPart, polOrder, Storage > BaseType;

      typedef typename BaseType::EntityType EntityType;

      typedef typename FunctionSpace::RangeType RangeType;

    public:
      typedef LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      static const int localBlockSize = BaseType::dimRange * NumLegendreBaseFunctions< polOrder, BaseType::dimLocal >::numBaseFct;
      typedef NonBlockMapper< typename BaseType::BlockMapperType, localBlockSize > MapperType;

      typedef LegendreShapeFunctionSet< typename ToLocalFunctionSpace< FunctionSpace, BaseType::dimLocal >::Type > ScalarShapeFunctionSetType;
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetType, RangeType > ShapeFunctionSetType;

      typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;
    };



    // LegendreDiscontinuousGalerkinSpace
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class LegendreDiscontinuousGalerkinSpace
    : public DiscontinuousGalerkinSpaceDefault< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscontinuousGalerkinSpaceDefault< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

    public:
      using BaseType::blockMapper;

      typedef typename BaseType::Traits Traits;
      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::MapperType MapperType;
      typedef typename BaseType::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename BaseType::EntityType EntityType;

      LegendreDiscontinuousGalerkinSpace ( const GridPartType &gridPart,
                                           const InterfaceType commInterface = BaseType::defaultInterface,
                                           const CommunicationDirection commDirection = BaseType::defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        mapper_( blockMapper() )
      {}

      ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity ) const
      {
        assert( entity.type() == GeometryType( GenericGeometry::CubeTopology< Traits::dimLocal >::type ) );
        return shapeFunctionSet( entity.type() );
      }

      ShapeFunctionSetType shapeFunctionSet ( const GeometryType &type) const;

      DUNE_VERSION_DEPRECATED(1,4,remove)
      MapperType &mapper () const { return mapper_; }

    private:
      mutable MapperType mapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH

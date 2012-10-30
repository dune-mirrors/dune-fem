#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LAGRANGE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LAGRANGE_HH

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/lagrange/genericbasefunctions.hh>

// local includes
#include "declaration.hh"
#include "default.hh"


namespace Dune
{

  namespace Fem
  {

    // LagrangeDiscontinuousGalerkinSpaceTraits
    // ----------------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class LagrangeDiscontinuousGalerkinSpaceTraits
    : public DiscontinuousGalerkinSpaceTraitsBase< FunctionSpace, GridPart, polOrder, Storage >
    {
      typedef DiscontinuousGalerkinSpaceTraitsBase< FunctionSpace, GridPart, polOrder, Storage > BaseType;

      typedef typename GeometryWrapper<
          Dune::Fem::GridPartCapabilities::hasSingleGeometryType< typename BaseType::GridPartType >::topologyId, BaseType::dimLocal
        >::GenericGeometryType GenericGeometryType;

      typedef GenericLagrangeBaseFunction<
          typename FunctionSpace::ScalarFunctionSpaceType, GenericGeometryType, polOrder
        > GenericBaseFunctionType;

    public:
      typedef LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      static const int localBlockSize = BaseType::dimRange * GenericBaseFunctionType::numBaseFunctions;
      typedef NonBlockMapper< typename BaseType::BlockMapperType, localBlockSize > MapperType;

      // ShapeFunctionSetType
      // BasisFunctionSetType
    };



    // LagrangeDiscontinuousGalerkinSpace
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class LagrangeDiscontinuousGalerkinSpace
    : public DiscontinuousGalerkinSpaceDefault< LagrangeDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {

      typedef LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscontinuousGalerkinSpaceDefault< LagrangeDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

    public:
      using BaseType::blockMapper;

      typedef typename BaseType::Traits Traits;
      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::MapperType MapperType;
      typedef typename BaseType::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename BaseType::EntityType EntityType;

      LagrangeDiscontinuousGalerkinSpace ( const GridPartType &gridPart,
                                           const InterfaceType commInterface = BaseType::defaultInterface,
                                           const CommunicationDirection commDirection = BaseType::defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        mapper_( blockMapper() )
      {}

      ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity ) const
      {
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

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LAGRANGE_HH

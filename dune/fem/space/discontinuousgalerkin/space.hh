#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/misc/bartonnackmaninterface.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/version.hh>

// local includes
#include "declaration.hh"
#include "default.hh"


namespace Dune
{

  namespace Fem
  {

    // DiscontinuousGalerkinSpaceTraits
    // --------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class DiscontinuousGalerkinSpaceTraits
    : public DiscontinuousGalerkinSpaceTraitsBase< FunctionSpace, GridPart, polOrder, Storage >
    {
      typedef DiscontinuousGalerkinSpaceTraitsBase< FunctionSpace, GridPart, polOrder, Storage > BaseType;

    public:
      typedef DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      static const int localBlockSize = BaseType::dimRange * DGNumberOfBaseFunctions< polOrder, BaseType::dimLocal >::numBaseFunctions;
      typedef NonBlockMapper< typename BaseType::BlockMapperType, localBlockSize > MapperType;

      // ShapeFunctionSetType
      // BasisFunctionSetType
    };



    // DiscontinuousGalerkinSpace
    // --------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class DiscontinuousGalerkinSpace
    : public DiscontinuousGalerkinSpaceDefault< DiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType; 
      typedef DiscontinuousGalerkinSpaceDefault< DiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

    public:
      using BaseType::blockMapper;

      typedef typename BaseType::Traits Traits;
      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::MapperType MapperType;
      typedef typename BaseType::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename BaseType::EntityType EntityType;

      DiscontinuousGalerkinSpace ( const GridPartType &gridPart,
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
      MapperType mapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH

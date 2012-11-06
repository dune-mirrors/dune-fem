#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
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

      /** \brief return shape function set for given entity
       *
       * \param[in]  entity  entity (of codim 0) for which shape function set 
       *                     is requested
       *
       * \returns  ShapeFunctionSetType  shape function set                     
       */
      ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity ) const
      {
        return shapeFunctionSet( entity.type() );
      }

      /** \brief return shape function set for geometry type 
       *
       * \param[in]  type  geometry type for which shape function set
       *                   is requested
       *
       * \returns  ShapeFunctionSetType  shape function set                     
       */
      ShapeFunctionSetType shapeFunctionSet ( const GeometryType &type) const;

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::mapper */
      DUNE_VERSION_DEPRECATED(1,4,remove)
      MapperType &mapper () const { return mapper_; }

    private:
      MapperType mapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH

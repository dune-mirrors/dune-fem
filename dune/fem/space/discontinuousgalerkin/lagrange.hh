#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LAGRANGE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LAGRANGE_HH

// C++ includes
#include <vector>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>
#include <dune/fem/space/lagrange/genericbasefunctions.hh>
#include <dune/fem/space/lagrange/shapefunctionset.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>

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

    public:
      typedef LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
      typedef typename BaseType::GridPartType GridPartType;

    private:
      typedef typename GridPartType::template Codim< BaseType::codimension >::EntityType EntityType;

      static const int dimLocal = GridPartType::dimension;
      typedef typename ToLocalFunctionSpace< FunctionSpaceType, dimLocal >::Type ShapeFunctionSpaceType;

    public:
      typedef LagrangeShapeFunctionSet< ShapeFunctionSpaceType, polOrder > ShapeFunctionSetImp;
      typedef ShapeFunctionSetProxy< ShapeFunctionSetImp > ShapeFunctionSetType;
      typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

    private:
      template< class GridPartType, bool hasSingleGeometryType = Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPartType >::v >
      struct HaveSingleGeometryType;

      template< class GridPartType >
      class HaveSingleGeometryType< GridPartType, true >
      {
        typedef typename GeometryWrapper<
            Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPartType >::topologyId, BaseType::dimLocal
          >::GenericGeometryType GenericGeometryType;
        typedef GenericLagrangeBaseFunction<
            typename FunctionSpaceType::ScalarFunctionSpaceType, GenericGeometryType, polOrder
          > GenericBaseFunctionType;

      public:
        static const int localBlockSize = BaseType::dimRange * GenericBaseFunctionType::numBaseFunctions;
        typedef NonBlockMapper< typename BaseType::BlockMapperType, localBlockSize > MapperType;
      };

    public:
      static const int localBlockSize = HaveSingleGeometryType< GridPartType >::localBlockSize;
      typedef typename HaveSingleGeometryType< GridPartType >::MapperType MapperType;
    };



    // LagrangeDiscontinuousGalerkinSpace
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    class LagrangeDiscontinuousGalerkinSpace
    : public DiscontinuousGalerkinSpaceDefault< LagrangeDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage >, Storage >
    {

      typedef LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscontinuousGalerkinSpaceDefault< LagrangeDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage >, Storage > BaseType;

    public:
      using BaseType::blockMapper;

      typedef typename BaseType::Traits Traits;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::MapperType MapperType;

      typedef typename Traits::ShapeFunctionSetType ShapeFunctionSetType;

    private:
      // shape function set is a proxy, get underlying type
      typedef typename ShapeFunctionSetType::ImplementationType ShapeFunctionSetImp;
      typedef SingletonList< const GeometryType, ShapeFunctionSetImp > SingletonProviderType;
      typedef BaseSetLocalKeyStorage< ShapeFunctionSetImp > ShapeFunctionSetStorageType;

    public:
      LagrangeDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                           const InterfaceType commInterface = BaseType::defaultInterface,
                                           const CommunicationDirection commDirection = BaseType::defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        mapper_( blockMapper() )
      {
        // get geometry types
        std::vector< GeometryType > geomTypes = AllGeomTypes< IndexSetType, GridType >( gridPart.indexSet()).geomTypes( BaseType::codimension );

        // store shape function sets per type
        const typename std::vector< GeometryType >::const_iterator end = geomTypes.end();
        for( typename std::vector< GeometryType >::const_iterator it = geomTypes.begin(); it != end; ++it )
        {
          const GeometryType &type = *it;
          shapeFunctionSets_.template insert< SingletonProviderType >( type );
        }
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet */
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity, shapeFunctionSet( entity ) );
      }

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
      ShapeFunctionSetType shapeFunctionSet ( const GeometryType &type) const
      {
        // wrap proxy around implementation
        return ShapeFunctionSetType( &shapeFunctionSets_[ type ] );
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::mapper */
      DUNE_VERSION_DEPRECATED(1,4,remove)
      MapperType &mapper () const { return mapper_; }

    private:
      mutable MapperType mapper_;
      ShapeFunctionSetStorageType shapeFunctionSets_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LAGRANGE_HH

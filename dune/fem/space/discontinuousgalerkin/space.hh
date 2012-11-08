#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH

// C++ includes
#include <vector>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/version.hh>

// local includes
#include "declaration.hh"
#include "default.hh"
#include "shapefunctionset.hh"


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

      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
      typedef typename BaseType::GridPartType GridPartType;

    private:
      typedef typename GridPartType::template Codim< BaseType::codimension >::EntityType EntityType;

      typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      typedef typename ToLocalFunctionSpace< ScalarFunctionSpaceType, BaseType::dimLocal >::Type ScalarShapeFunctionSpaceType;

    public:
      typedef OrthonormalShapeFunctionSet< ScalarShapeFunctionSpaceType, polOrder > ScalarShapeFunctionSetType;
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetType, typename FunctionSpaceType::RangeType > ShapeFunctionSetImp;
      typedef ShapeFunctionSetProxy< ShapeFunctionSetImp > ShapeFunctionSetType;

      static const int localBlockSize = BaseType::dimRange * OrthonormalShapeFunctionSetSize< ScalarShapeFunctionSpaceType, polOrder >::v;
      typedef NonBlockMapper< typename BaseType::BlockMapperType, localBlockSize > MapperType;

      typedef DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;
    };



    // DiscontinuousGalerkinSpace
    // --------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = SimpleStorage >
    class DiscontinuousGalerkinSpace
    : public DiscontinuousGalerkinSpaceDefault< DiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType; 
      typedef DiscontinuousGalerkinSpaceDefault< DiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

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
      DiscontinuousGalerkinSpace ( GridPartType &gridPart,
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

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH

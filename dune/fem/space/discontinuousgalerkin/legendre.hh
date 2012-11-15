#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH

// C++ includes
#include <cassert>

// dune-common includes
#include <dune/common/power.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>

// dune-fem includes
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/legendre.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

// local includes
#include "default.hh"


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

    public:
      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
      typedef typename BaseType::GridPartType GridPartType;

    private:
      typedef typename GridPartType::template Codim< BaseType::codimension >::EntityType EntityType;

      typedef typename FunctionSpaceType::RangeType RangeType;

      template <int p, int dim>
      struct NumLegendreShapeFunctions
      {
        static const int v = StaticPower< p+1, dim >::power;
      };

    public:
      typedef LegendreDiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder, Storage > DiscreteFunctionSpaceType;

      static const int localBlockSize = BaseType::dimRange * NumLegendreShapeFunctions< polOrder, BaseType::dimLocal >::v;
      typedef NonBlockMapper< typename BaseType::BlockMapperType, localBlockSize > MapperType;

      typedef LegendreShapeFunctionSet< typename ToLocalFunctionSpace< FunctionSpace, BaseType::dimLocal >::Type > LegendreShapeFunctionSetType;
      typedef SelectCachingShapeFunctionSet< LegendreShapeFunctionSetType, Storage > ShapeFunctionSetImp;

      struct ShapeFunctionSetFactory
      {
        static ShapeFunctionSetImp *createObject ( const GeometryType &type )
        {
          return new ShapeFunctionSetImp( type, LegendreShapeFunctionSetType( polOrder ) );
        }

        static void deleteObject ( ShapeFunctionSetImp *object ) { delete object; }
      };

      typedef ShapeFunctionSetProxy< ShapeFunctionSetImp > ScalarShapeFunctionSetType;
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetType, RangeType > ShapeFunctionSetType;
      typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;
    };



    // LegendreDiscontinuousGalerkinSpace
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    class LegendreDiscontinuousGalerkinSpace
    : public DiscontinuousGalerkinSpaceDefault< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage >, Storage >
    {
      typedef LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscontinuousGalerkinSpaceDefault< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage >, Storage > BaseType;

    public:
      using BaseType::blockMapper;

      static const int polynomialOrder DUNE_DEPRECATED = polOrder;

      typedef typename BaseType::Traits Traits;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::MapperType MapperType;

      typedef typename Traits::ShapeFunctionSetType ShapeFunctionSetType;

    private:
      typedef typename ShapeFunctionSetType::ImplementationType ShapeFunctionSetImp;
      typedef SingletonList< const GeometryType, ShapeFunctionSetImp, typename Traits::ShapeFunctionSetFactory > SingletonProviderType;
      typedef BaseSetLocalKeyStorage< ShapeFunctionSetImp > ShapeFunctionSetStorageType;

    public:
      LegendreDiscontinuousGalerkinSpace ( GridPartType &gridPart,
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

      /** \brief return shape unique function set for geometry type 
       *
       * \param[in]  type  geometry type (must be a cube) for which 
       *                   shape function set is requested
       *
       * \returns  ShapeFunctionSetType  shape function set                     
       */
      ShapeFunctionSetType shapeFunctionSet ( const GeometryType &type) const
      {
        return ShapeFunctionSetType( &shapeFunctionSets_ );
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

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH

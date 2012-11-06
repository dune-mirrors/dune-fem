#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH

// C++ includes
#include <cassert>

// dune-common includes
#include <dune/common/power.hh>
#include <dune/common/typetraits.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>

// dune-fem includes
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/legendre.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

// local includes
#include "declaration.hh"
#include "default.hh"


namespace
{
  template< template< class > class Storage >
  struct ShowWarning;

  template<>
  struct ShowWarning< Dune::Fem::CachingStorage >
  {
    static const bool v = true;
  };

  template<>
  struct ShowWarning< Dune::Fem::SimpleStorage >
  {
    static const bool v = false;
  };

} // namespace



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

      template <int p, int dim>
      struct NumLegendreShapeFunctions
      {
        static const int v = StaticPower< p+1, dim >::power;
      };

    public:
      typedef LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      static const int localBlockSize = BaseType::dimRange * NumLegendreShapeFunctions< polOrder, BaseType::dimLocal >::v;
      typedef NonBlockMapper< typename BaseType::BlockMapperType, localBlockSize > MapperType;

      typedef LegendreShapeFunctionSet< typename ToLocalFunctionSpace< FunctionSpace, BaseType::dimLocal >::Type > ScalarShapeFunctionSetType;
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetType, RangeType > ShapeFunctionSetImp;
      typedef ShapeFunctionSetProxy< ShapeFunctionSetImp > ShapeFunctionSetType;

      typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;
    };



    // LegendreDiscontinuousGalerkinSpace
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    class LegendreDiscontinuousGalerkinSpace
    : public DiscontinuousGalerkinSpaceDefault< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscontinuousGalerkinSpaceDefault< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

    public:
      using BaseType::blockMapper;

      static const int polynomialOrder DUNE_DEPRECATED = polOrder;

      typedef typename BaseType::Traits Traits;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::MapperType MapperType;
      typedef typename BaseType::ShapeFunctionSetType ShapeFunctionSetType;

    private:
      typedef typename ShapeFunctionSetType::ImplementationType ShapeFunctionSetImp;

    public:
      LegendreDiscontinuousGalerkinSpace ( const GridPartType &gridPart,
                                           const InterfaceType commInterface = BaseType::defaultInterface,
                                           const CommunicationDirection commDirection = BaseType::defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        mapper_( blockMapper() ),
        shapeFunctionSet_( polOrder )
      {
        deprecationWarning( Dune::integral_constant< bool, ShowWarning< Storage >::v >() );
      }

      /** \brief return unique shape function set for given entity
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
        assert( entity.type() == GeometryType( GenericGeometry::CubeTopology< Traits::dimLocal >::type ) );
        return ShapeFunctionSetType( &shapeFunctionSet_ );
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::mapper */
      DUNE_VERSION_DEPRECATED(1,4,remove)
      MapperType &mapper () const { return mapper_; }

    private:
      void DUNE_DEPRECATED_MSG( "Caching disabled for FiniteVolumeSpace." )
      deprecationWarning ( Dune::integral_constant< bool, true > ) {}
      void
      deprecationWarning ( Dune::integral_constant< bool, false > ) {}

      mutable MapperType mapper_;
      ShapeFunctionSetImp shapeFunctionSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH

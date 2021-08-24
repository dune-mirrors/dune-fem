#ifndef DUNE_FEM_SPACE_LAGRANGE_SPACE_HH
#define DUNE_FEM_SPACE_LAGRANGE_SPACE_HH

// C++ includes
#include <vector>

// dune-common includes
#include <dune/common/exceptions.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/common/hybrid.hh>

#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/localinterpolation.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/mapper/indexsetdofmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>
#include <dune/fem/storage/singletonlist.hh>

// local includes
#include "adaptmanager.hh"
#include "capabilities.hh"
#include "interpolation.hh"
#include "lagrangepoints.hh"
#include "shapefunctionset.hh"
#include "storage.hh"


namespace Dune
{

  namespace Fem
  {

    // Forward declaration
    // -------------------

    template< class FunctionSpace, class GridPart, int maxPolOrder, class Storage = CachingStorage >
    class LagrangeDiscreteFunctionSpace;



    // LagrangeDiscreteFunctionSpaceTraits
    // -----------------------------------

    template< class FunctionSpace, class GridPart, int maxPolOrder, class Storage >
    struct LagrangeDiscreteFunctionSpaceTraits
    {
      typedef LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, maxPolOrder, Storage > DiscreteFunctionSpaceType;

      // we take 6 as the maximal polynomial order for the dynamic space
      static const int maxPolynomialOrder = ( maxPolOrder < 0 ) ? -maxPolOrder : maxPolOrder;
      static const int minPolynomialOrder = ( maxPolOrder < 0 ) ? 1 : maxPolOrder;

      typedef GridPart GridPartType;
      typedef GridFunctionSpace< GridPartType, FunctionSpace > FunctionSpaceType;

      typedef IndexSetDofMapper< GridPartType,
                                 std::conditional_t< Dune::Capabilities::isCartesian< typename GridPart::GridType >::v,
                                                     DefaultLocalDofMapping< GridPart >,
                                                     LagrangeLocalDofMapping< GridPart > > > BlockMapperType;

      typedef Hybrid::IndexRange< int, FunctionSpaceType::dimRange > LocalBlockIndices;

      static const int codimension = 0;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;
      static const int dimLocal = GridPartType::dimension;
      typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      typedef typename ToNewDimDomainFunctionSpace< ScalarFunctionSpaceType, dimLocal >::Type ShapeFunctionSpaceType;

    public:
      typedef LagrangeShapeFunctionSet< ShapeFunctionSpaceType, maxPolynomialOrder > LagrangeShapeFunctionSetType;
      typedef SelectCachingShapeFunctionSet< LagrangeShapeFunctionSetType, Storage > ScalarShapeFunctionSetType;

      typedef ShapeFunctionSetProxy< ScalarShapeFunctionSetType > ScalarShapeFunctionSetProxyType;
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetProxyType, typename FunctionSpaceType::RangeType > ShapeFunctionSetType;

      typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

      template< class DiscreteFunction, class Operation = Dune::Fem::DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        // type of data handle
        typedef Dune::Fem::DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        // type of operation to perform on scatter
        typedef Operation OperationType;
      };
    };



    /** \addtogroup LagrangeDiscreteFunctionSpace
     *
     *  Provides access to bse function sets for different element types in
     *  one grid and size of function space and maps from local to global dof
     *  number.
     *
     *  \note This space can only be used with special index sets. If you want
     *  to use the LagrangeDiscreteFunctionSpace with an index set only
     *  supporting the index set interface you will have to use the
     *  IndexSetWrapper class to provide the required functionality.
     *
     *  \note For adaptive calculations one has to use index sets that are
     *  capable of adaption (i.e. the method adaptive returns true). See also
     *  AdaptiveLeafIndexSet.
     */

    // LagrangeDiscreteFunctionSpace
    // -----------------------------

    /** \class   LagrangeDiscreteFunctionSpace
     *  \ingroup LagrangeDiscreteFunctionSpace
     *  \brief   Lagrange discrete function space
     */

    template< class FunctionSpace, class GridPart, int maxPolOrder, class Storage >
    class LagrangeDiscreteFunctionSpace
    : public DiscreteFunctionSpaceDefault< LagrangeDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, maxPolOrder, Storage > >
    {

      typedef LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, maxPolOrder, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< LagrangeDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, maxPolOrder, Storage > > BaseType;

    public:
      typedef typename BaseType::Traits Traits;
      static const int maxPolynomialOrder = Traits :: maxPolynomialOrder;
      static const int minPolynomialOrder = Traits :: minPolynomialOrder;

      static_assert( (maxPolynomialOrder > 0), "LagrangeDiscreteFunctionSpace only defined for polOrder > 0" );

      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename BaseType::Traits::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::BlockMapperType BlockMapperType;

      typedef LagrangePointSet< GridPartType, maxPolynomialOrder > LagrangePointSetType;
      typedef LagrangeLocalInterpolation< GridPartType, maxPolynomialOrder, BasisFunctionSetType >  LocalInterpolationType;
      typedef LocalInterpolationType InterpolationImplType;

      typedef LocalInterpolationWrapper< ThisType >  InterpolationType;

    private:
      typedef typename Traits::ScalarShapeFunctionSetType ScalarShapeFunctionSetType;
      typedef BaseSetLocalKeyStorage< ScalarShapeFunctionSetType > ScalarShapeFunctionSetStorageType;

      typedef CompiledLocalKeyContainer< LagrangePointSetType, minPolynomialOrder, maxPolynomialOrder > LagrangePointSetContainerType;
      typedef typename LagrangePointSetContainerType::LocalKeyStorageType LocalKeyStorageType;

      typedef LagrangeMapperSingletonKey< GridPartType, LocalKeyStorageType >
        MapperSingletonKeyType;
      typedef LagrangeMapperSingletonFactory< MapperSingletonKeyType, BlockMapperType >
        BlockMapperSingletonFactoryType;
      typedef SingletonList< MapperSingletonKeyType, BlockMapperType, BlockMapperSingletonFactoryType > BlockMapperProviderType;

      // static const InterfaceType defaultInterface = InteriorBorder_InteriorBorder_Interface;
      static const InterfaceType defaultInterface = GridPart::indexSetInterfaceType;
      static const CommunicationDirection defaultDirection = ForwardCommunication;


      template < int pOrd >
      struct Initialize
      {
        struct ScalarShapeFunctionSetFactory
        {
        public:
          static ScalarShapeFunctionSetType *createObject ( const GeometryType &type )
          {
            typedef typename Traits :: LagrangeShapeFunctionSetType  LagrangeShapeFunctionSetType;
            return new ScalarShapeFunctionSetType( type, LagrangeShapeFunctionSetType( type, pOrd ) );
          }

          static void deleteObject ( ScalarShapeFunctionSetType *object ) { delete object; }
        };

        static void apply ( ScalarShapeFunctionSetStorageType& scalarShapeFunctionSets,
                            const int polynomialOrder,
                            const GeometryType &type )
        {
          if( polynomialOrder == pOrd )
          {
            typedef SingletonList< const GeometryType, ScalarShapeFunctionSetType, ScalarShapeFunctionSetFactory > SingletonProviderType;
            scalarShapeFunctionSets.template insert< SingletonProviderType >( type );
          }
        }
      };


    public:
      ///////////////////////
      // Interface methods //
      ///////////////////////

      using BaseType::order;

      explicit LagrangeDiscreteFunctionSpace ( GridPartType &gridPart,
                                               const InterfaceType commInterface,
                                               const CommunicationDirection commDirection )
      : LagrangeDiscreteFunctionSpace( gridPart, minPolynomialOrder, commInterface, commDirection )
      {}

      explicit LagrangeDiscreteFunctionSpace ( GridPartType &gridPart,
                                               const int polOrder = minPolynomialOrder,
                                               const InterfaceType commInterface = defaultInterface,
                                               const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        blockMapper_(),
        lagrangePointSetContainer_( gridPart ),
        // when min == max polynomial order the dynamic choice is off
        polynomialOrder_( ( minPolynomialOrder == maxPolynomialOrder ) ? minPolynomialOrder : polOrder )
      {
        const IndexSetType &indexSet = gridPart.indexSet();

        AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );
        const std::vector< GeometryType > &geometryTypes = allGeometryTypes.geomTypes( 0 );

        // create shape function sets
        for( unsigned int i = 0; i < geometryTypes.size(); ++i )
        {
          Fem::ForLoop< Initialize, minPolynomialOrder, maxPolynomialOrder >::
            apply( scalarShapeFunctionSets_, polynomialOrder_, geometryTypes[ i ] );
        }

        MapperSingletonKeyType key( gridPart, lagrangePointSetContainer_.compiledLocalKeys( polynomialOrder_ ), polynomialOrder_ );
        blockMapper_.reset( &BlockMapperProviderType::getObject( key ) );
        assert( blockMapper_ );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      DFSpaceIdentifier type () const
      {
        return LagrangeSpace_id;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet */
      const BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity, shapeFunctionSet( entity ) );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous () const
      {
        return true;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType &intersection ) const
      {
        return intersection.conforming();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const
      {
        return polynomialOrder_;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const
      {
        assert( blockMapper_ );
        return *blockMapper_;
      }

      ///////////////////////////
      // Non-interface methods //
      ///////////////////////////

      /** \brief return interpolation object
       */
      InterpolationType interpolation () const
      {
        return InterpolationType( *this );
      }

      /** \brief return local interpolation for given entity
       *
       *  \param[in]  entity  grid part entity
       */
      [[deprecated]]
      LocalInterpolationType interpolation ( const EntityType &entity ) const
      {
        return LocalInterpolationType( lagrangePointSetContainer_.compiledLocalKey( entity.type(), polynomialOrder_ ), basisFunctionSet( entity ) );
      }

      /** \brief return local interpolation for given entity
       *
       *  \param[in]  entity  grid part entity
       */
      LocalInterpolationType localInterpolation ( const EntityType &entity ) const
      {
        return LocalInterpolationType( lagrangePointSetContainer_.compiledLocalKey( entity.type(), polynomialOrder_ ), basisFunctionSet( entity ) );
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
      ShapeFunctionSetType shapeFunctionSet ( const GeometryType &type ) const
      {
        return ShapeFunctionSetType( &scalarShapeFunctionSets_[ type ] );
      }

      LagrangeDiscreteFunctionSpace ( const ThisType & ) = delete;
      ThisType &operator= ( const ThisType & ) = delete;

    protected:
      std::unique_ptr< BlockMapperType, typename BlockMapperProviderType::Deleter> blockMapper_;
      ScalarShapeFunctionSetStorageType scalarShapeFunctionSets_;
      LagrangePointSetContainerType lagrangePointSetContainer_;
      const int polynomialOrder_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_SPACE_HH

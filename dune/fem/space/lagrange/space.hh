#ifndef DUNE_FEM_SPACE_LAGRANGE_SPACE_HH
#define DUNE_FEM_SPACE_LAGRANGE_SPACE_HH

// C++ includes
#include <vector>

// dune-common includes
#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/nullptr.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/basisfunctionset/default.hh>
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

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    class LagrangeDiscreteFunctionSpace;



    // LagrangeDiscreteFunctionSpaceTraits
    // -----------------------------------

    template< class FunctionSpace, class GridPart, unsigned int polOrder, template< class > class Storage >
    struct LagrangeDiscreteFunctionSpaceTraits
    {
      typedef LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int polynomialOrder = polOrder;

      static const int localBlockSize = FunctionSpaceType::dimRange;
      typedef IndexSetDofMapper< GridPartType > BlockMapperType;

      static const int codimension = 0;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;
      static const int dimLocal = GridPartType::dimension;
      typedef typename FunctionSpace::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      typedef typename ToNewDimDomainFunctionSpace< ScalarFunctionSpaceType, dimLocal >::Type ShapeFunctionSpaceType;

    public:
      typedef LagrangeShapeFunctionSet< ShapeFunctionSpaceType, polynomialOrder > LagrangeShapeFunctionSetType;
      typedef SelectCachingShapeFunctionSet< LagrangeShapeFunctionSetType, Storage > ScalarShapeFunctionSetType;

      struct ScalarShapeFunctionSetFactory
      {
        static ScalarShapeFunctionSetType *createObject ( const GeometryType &type )
        {
          return new ScalarShapeFunctionSetType( type, LagrangeShapeFunctionSetType( type ) );
        }

        static void deleteObject ( ScalarShapeFunctionSetType *object ) { delete object; }
      };
      typedef ScalarShapeFunctionSetFactory ScalarShapeFunctionSetFactoryType;

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

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class LagrangeDiscreteFunctionSpace
    : public DiscreteFunctionSpaceDefault< LagrangeDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      static_assert( (polOrder > 0), "LagrangeDiscreteFunctionSpace only defined for polOrder > 0" );

      typedef LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< LagrangeDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

    public:
      typedef typename BaseType::Traits Traits;
      static const int polynomialOrder = polOrder;

      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename BaseType::Traits::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::BlockMapperType BlockMapperType;

      typedef LagrangePointSet< GridPartType, polynomialOrder > LagrangePointSetType;
      typedef LagrangeLocalInterpolation< GridPartType, polynomialOrder, BasisFunctionSetType > InterpolationType;

    private:
      typedef typename Traits::ScalarShapeFunctionSetType ScalarShapeFunctionSetType;
      typedef SingletonList< GeometryType, ScalarShapeFunctionSetType, typename Traits::ScalarShapeFunctionSetFactoryType > SingletonProviderType;
      typedef BaseSetLocalKeyStorage< ScalarShapeFunctionSetType > ScalarShapeFunctionSetStorageType;

      typedef CompiledLocalKeyContainer< LagrangePointSetType, polynomialOrder, polynomialOrder > LagrangePointSetContainerType;
      typedef typename LagrangePointSetContainerType::LocalKeyStorageType LocalKeyStorageType;

      typedef LagrangeMapperSingletonKey< GridPartType, LocalKeyStorageType >
        MapperSingletonKeyType;
      typedef LagrangeMapperSingletonFactory< MapperSingletonKeyType, BlockMapperType >
        BlockMapperSingletonFactoryType;
      typedef SingletonList< MapperSingletonKeyType, BlockMapperType, BlockMapperSingletonFactoryType > BlockMapperProviderType;

      // static const InterfaceType defaultInterface = InteriorBorder_InteriorBorder_Interface;
      static const InterfaceType defaultInterface = GridPart::indexSetInterfaceType;
      static const CommunicationDirection defaultDirection = ForwardCommunication;

    public:
      ///////////////////////
      // Interface methods //
      ///////////////////////

      using BaseType::order;

      explicit LagrangeDiscreteFunctionSpace ( GridPartType &gridPart,
                                               const InterfaceType commInterface = defaultInterface,
                                               const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        blockMapper_( nullptr ),
        lagrangePointSetContainer_( gridPart )
      {
        const IndexSetType &indexSet = gridPart.indexSet();

        AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );
        const std::vector< GeometryType > &geometryTypes = allGeometryTypes.geomTypes( 0 );
        for( unsigned int i = 0; i < geometryTypes.size(); ++i )
        {
          const GeometryType &type = geometryTypes[ i ];
          scalarShapeFunctionSets_.template insert< SingletonProviderType >( type );
        }

        MapperSingletonKeyType key( gridPart, lagrangePointSetContainer_.compiledLocalKeys( polynomialOrder ), polynomialOrder );
        blockMapper_ = &BlockMapperProviderType::getObject( key );
        assert( blockMapper_ );
      }

      ~LagrangeDiscreteFunctionSpace ()
      {
        BlockMapperProviderType::removeObject( *blockMapper_ );
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
        return polOrder;
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

      /** \brief return local interpolation for given entity
       *
       *  \param[in]  entity  grid part entity
       */
      InterpolationType interpolation ( const EntityType &entity ) const
      {
        return InterpolationType( lagrangePointSet( entity ), basisFunctionSet( entity ) );
      }

      /** \brief interpolate a function locally
       *
       *  \param[in]  loalFunction  local function to interpolate
       *  \param[out]  localDofVector local degrees of freedom of the interpolion
       */
      template< class LocalFunction, class LocalDofVector >
      DUNE_DEPRECATED
      void interpolate ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        const EntityType &entity = localFunction.entity();
        const auto interpolation = this->interpolation( entity );
        interpolation( localFunction, localDofVector );
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

      /** \brief provide access to the Lagrange point set for an entity
       *
       *  \note This method is not part of the DiscreteFunctionSpaceInterface. It
       *        is unique to the LagrangeDiscreteFunctionSpace.
       *
       *  \param[in]  entity  entity the Lagrange point set is requested for
       *
       *  \returns LagrangePointSet
       */
      template< class EntityType >
      const LagrangePointSetType &lagrangePointSet ( const EntityType &entity ) const
      {
        return lagrangePointSet( entity.type() );
      }

      /** \brief provide access to the Lagrange point set for a geometry type
       *
       *  \note This method is not part of the DiscreteFunctionSpaceInterface. It
       *        is unique to the LagrangeDiscreteFunctionSpace.
       *
       *  \param[in]  type  type of geometry the Lagrange point set is requested for
       *
       *  \returns LagrangePointSetType
       */
      const LagrangePointSetType &lagrangePointSet ( const GeometryType &type ) const
      {
        return lagrangePointSetContainer_.compiledLocalKey( type, polynomialOrder );
      }

    private:
      // forbid copying
      LagrangeDiscreteFunctionSpace ( const ThisType & );
      // forbid assignment
      ThisType &operator= ( const ThisType & );

      BlockMapperType *blockMapper_;
      ScalarShapeFunctionSetStorageType scalarShapeFunctionSets_;
      LagrangePointSetContainerType lagrangePointSetContainer_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_SPACE_HH

#ifndef DUNE_FEM_SPACE_LAGRANGE_SPACE_HH
#define DUNE_FEM_SPACE_LAGRANGE_SPACE_HH

// C++ includes
#include <vector>

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/nullptr.hh>
#include <dune/common/static_assert.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/dofmapper/indexsetdofmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>
#include <dune/fem/storage/singletonlist.hh>

// local includes
#include "lagrangepoints.hh"
#include "adaptmanager.hh"
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
      // type of discrete function space
      typedef LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      // function space type
      typedef FunctionSpace FunctionSpaceType;
      // grid part type
      typedef GridPart GridPartType;
     
      // export template parameter polynomial order
      static const int polynomialOrder = polOrder;

      // local block size
      static const int localBlockSize = FunctionSpaceType::dimRange;
      // block mapper type
      typedef IndexSetDofMapper< GridPartType > BlockMapperType;
      // mapper type
      typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;

      // codimension
      static const int codimension = 0;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;
      static const int dimLocal = GridPartType::dimension;
      typedef typename FunctionSpace::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      typedef typename ToLocalFunctionSpace< ScalarFunctionSpaceType, dimLocal >::Type ShapeFunctionSpaceType;

    public:
      // scalar shape function set type
      typedef LagrangeShapeFunctionSet< ShapeFunctionSpaceType, polynomialOrder > ScalarShapeFunctionSetType;
      // shape function set
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetType, typename FunctionSpaceType::RangeType > ShapeFunctionSetType;
      // proxy
      typedef ShapeFunctionSetProxy< ShapeFunctionSetType > ShapeFunctionSetProxyType;
      // basis function set type
      typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ShapeFunctionSetProxyType > BasisFunctionSetType;

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
      dune_static_assert( (polOrder > 0), "LagrangeDiscreteFunctionSpace only defined for polOrder > 0" );

      typedef LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< LagrangeDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

    public:
      static const int polynomialOrder = polOrder;

      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
      typedef typename BaseType::IntersectionType IntersectionType;
      typedef typename BaseType::MapperType MapperType;
      typedef typename BaseType::BlockMapperType BlockMapperType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::EntityType EntityType;

      typedef LagrangePointSet< GridPartType, polynomialOrder > LagrangePointSetType;

    private:
      typedef typename BaseType::Traits::ShapeFunctionSetType ShapeFunctionSetType;
      typedef Fem::BaseSetLocalKeyStorage< ShapeFunctionSetType > ShapeSetStorageType;
      typedef SingletonList< GeometryType, ShapeFunctionSetType > ShapeFunctionSetSingletonProviderType;

      typedef CompiledLocalKeyContainer< LagrangePointSetType, polynomialOrder, polynomialOrder > LagrangePointSetContainerType;
      typedef typename LagrangePointSetContainerType::LocalKeyStorageType LocalKeyStorageType;

      typedef LagrangeMapperSingletonKey< GridPartType, LocalKeyStorageType >
        MapperSingletonKeyType;
      typedef LagrangeMapperSingletonFactory< MapperSingletonKeyType, BlockMapperType >
        BlockMapperSingletonFactoryType;
      typedef SingletonList< MapperSingletonKeyType, BlockMapperType, BlockMapperSingletonFactoryType > BlockMapperProviderType;

      static const InterfaceType defaultInterface = InteriorBorder_InteriorBorder_Interface;
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
        lagrangePointSetContainer_( gridPart ),
        blockMapper_( nullptr ),
        mapper_( nullptr )
      {
        const IndexSetType &indexSet = gridPart.indexSet();

        AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );
        const std::vector< GeometryType > &geometryTypes = allGeometryTypes.geomTypes( 0 );
        for( unsigned int i = 0; i < geometryTypes.size(); ++i )
        {
          const GeometryType &type = geometryTypes[ i ];
          shapeFunctionSets_.template insert< ShapeFunctionSetSingletonProviderType >( type );
        }

        MapperSingletonKeyType key( gridPart, lagrangePointSetContainer_.compiledLocalKeys( polynomialOrder ), polynomialOrder );
        blockMapper_ = &BlockMapperProviderType::getObject( key );
        assert( blockMapper_ );
        
        mapper_ = new MapperType( *blockMapper_ );
        assert( mapper_ );
      }

      ~LagrangeDiscreteFunctionSpace ()
      {
        if( mapper_ )
          delete mapper_;
        mapper_ = nullptr;
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
        return BasisFunctionSetType( entity, &shapeFunctionSets_[ entity.type() ] );
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

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::mapper */
      MapperType &mapper () const
      {
        assert( mapper_ );
        return *mapper_;
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

      ShapeSetStorageType shapeFunctionSets_;
      LagrangePointSetContainerType lagrangePointSetContainer_;
      BlockMapperType *blockMapper_;
      MapperType *mapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_SPACE_HH

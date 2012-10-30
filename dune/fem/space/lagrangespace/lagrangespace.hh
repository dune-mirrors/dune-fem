#ifndef DUNE_FEM_LAGRANGESPACE_LAGRANGESPACE_HH
#define DUNE_FEM_LAGRANGESPACE_LAGRANGESPACE_HH

#include <algorithm>
#include <vector>

#include <dune/common/misc.hh>
#include <dune/common/nullptr.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/space/basefunctions/basefunctionproxy.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/dofmapper/indexsetdofmapper.hh>
#include <dune/fem/space/lagrangespace/basefunctions.hh>
#include <dune/fem/space/lagrangespace/dofmappercode.hh>
#include <dune/fem/space/lagrangespace/lagrangepoints.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

namespace Dune
{

  namespace Fem 
  {

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    class LagrangeDiscreteFunctionSpace;


    template< class FunctionSpace, class GridPart, unsigned int polOrder,
              template< class > class BaseFunctionStorage = CachingStorage >
    struct LagrangeDiscreteFunctionSpaceTraits
    {
      static const int codimension = 0; 

      dune_static_assert((polOrder > 0), "LagrangeSpace only defined for polOrder > 0" );
      
      typedef FunctionSpace FunctionSpaceType;
      typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType :: DomainType DomainType;
      typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
      typedef typename FunctionSpaceType :: RangeType RangeType;
      typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

      static const int dimRange = FunctionSpaceType::dimRange;
      
      typedef GridPart GridPartType;
      typedef typename GridPartType::GridType GridType;
      typedef typename GridPartType::IndexSetType IndexSetType;

      typedef typename GridPartType::template Codim< codimension >::IteratorType IteratorType;
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      // get dimension of local coordinate 
      static const int dimLocal = GridPartType::dimension;

      typedef typename ToLocalFunctionSpace< FunctionSpaceType, dimLocal > :: Type 
        BaseFunctionSpaceType;

      static const int polynomialOrder = polOrder;
      
      typedef LagrangeDiscreteFunctionSpace
        < FunctionSpaceType, GridPartType, polynomialOrder, BaseFunctionStorage >
        DiscreteFunctionSpaceType;

      static const int localBlockSize = dimRange;

      // mapper for block
      typedef IndexSetDofMapper< GridPartType > BlockMapperType;
      typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;
      
      // implementation of basefunction set 
      typedef VectorialBaseFunctionSet< BaseFunctionSpaceType, BaseFunctionStorage >
          ShapeFunctionSetType;

      // exported type 
      typedef SimpleBaseFunctionProxy< ShapeFunctionSetType >  BaseFunctionSetType;
      typedef BaseFunctionSetType BasisFunctionSetType;

      /** \brief defines type of communication data handle for this type of space
       */
      template< class DiscreteFunction,
                class Operation = DFCommunicationOperation :: Add >
      struct CommDataHandle
      {
        //! type of data handle 
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        //! type of operatation to perform on scatter 
        typedef Operation OperationType;
      };
    };


    // LagrangeMapperSingletonKey
    // --------------------------

    template< class GridPart, class LagrangePointSetContainer >
    class LagrangeMapperSingletonKey 
    {
      typedef LagrangeMapperSingletonKey< GridPart, LagrangePointSetContainer > ThisType;

    public:
      typedef GridPart GridPartType;
      typedef LagrangePointSetContainer LagrangePointSetContainerType;

      LagrangeMapperSingletonKey ( const GridPartType &gridPart,
                                   const LagrangePointSetContainerType &pointSet,
                                   const int polOrd )
      : gridPart_( gridPart ),
        pointSet_( pointSet ),
        polOrd_( polOrd )
      {}

      bool operator== ( const ThisType &other ) const
      {
        return ((&indexSet() == &other.indexSet()) && (polOrd_ == other.polOrd_));
      }

      bool operator!= ( const ThisType &other ) const
      {
        return ((&indexSet() != &other.indexSet()) || (polOrd_ != other.polOrd_));
      }

      const GridPartType &gridPart () const { return gridPart_; }

      const typename GridPartType::IndexSetType &indexSet () const { return gridPart_.indexSet(); }

      const LagrangePointSetContainerType &pointSet () const { return pointSet_; }

    private:
      const GridPartType &gridPart_; 
      const LagrangePointSetContainerType &pointSet_;
      const int polOrd_;
    };

    
    // LagrangeMapperSingletonFactory
    // ------------------------------

    template< class Key, class Object >
    struct LagrangeMapperSingletonFactory
    {
      typedef typename Key::LagrangePointSetContainerType LagrangePointSetContainerType;

      static Object *createObject ( const Key &key )
      {
        LagrangeDofMapperCodeFactory< LagrangePointSetContainerType > codeFactory( key.pointSet() );
        return new Object( key.gridPart(), codeFactory );
      }
      
      static void deleteObject ( Object *obj )
      {
        delete obj;
      }
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
      typedef LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< LagrangeDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

    public:
      //! traits for the discrete function space
      typedef LagrangeDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > Traits;

      //! type of the discrete function space
      typedef ThisType LagrangeDiscreteFunctionSpaceType;

      typedef typename Traits::GridPartType GridPartType;
      typedef typename Traits::GridType GridType;
      typedef typename Traits::IndexSetType IndexSetType;
      typedef typename Traits::IteratorType IteratorType;

      //! type of intersections
      typedef typename BaseType ::IntersectionType IntersectionType;

      //! dimension of the grid (not the world)
      static const int dimension = GridPartType::dimension;

      typedef typename Traits::FunctionSpaceType FunctionSpaceType;
      //! field type for function space's domain
      typedef typename Traits::DomainFieldType DomainFieldType;
      //! type for function space's domain
      typedef typename Traits::DomainType DomainType;
      //! field type for function space's range
      typedef typename Traits::RangeFieldType RangeFieldType;
      //! type for function space's range
      typedef typename Traits::RangeType RangeType;
      //! dimension of function space's range
      static const int dimRange = Traits::dimRange;
      //! type of scalar function space
      typedef typename Traits::BaseFunctionSpaceType BaseFunctionSpaceType;
     
      //! maximum polynomial order of functions in this space
      static const int polynomialOrder = Traits::polynomialOrder;
      
      //! type of the base function set(s)
      typedef typename Traits::ShapeFunctionSetType ShapeFunctionSetType;
      // deprecated name 
      typedef ShapeFunctionSetType BaseFunctionSetImp;

      typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
      //! type of the shape function storage 
      typedef Fem::BaseSetLocalKeyStorage< ShapeFunctionSetType > ShapeSetStorageType;

      //! type of base function factory
      typedef LagrangeBaseFunctionFactory< typename BaseFunctionSpaceType::ScalarFunctionSpaceType, dimension, polynomialOrder >
        ScalarFactoryType;
      //! type of singleton base function factory
      typedef BaseFunctionSetSingletonFactory< GeometryType, ShapeFunctionSetType, ScalarFactoryType >
        BaseFunctionSetSingletonFactoryType;
      //! type of singleton list (singleton provider) for base functions
      typedef SingletonList< GeometryType, ShapeFunctionSetType, BaseFunctionSetSingletonFactoryType >
        BaseFunctionSetSingletonProviderType;

      //! type of a Lagrange point set
      typedef LagrangePointSet< GridPartType, polynomialOrder > LagrangePointSetType;

      // type of container for the LagangePointSets 
      typedef CompiledLocalKeyContainer< LagrangePointSetType, 
                  polynomialOrder, polynomialOrder > LagrangePointSetContainerType ;

      // type of local keys for one polynomial order 
      typedef typename LagrangePointSetContainerType :: LocalKeyStorageType  LocalKeyStorageType;

    public:
      //! mapper used to implement mapToGlobal
      typedef typename Traits::MapperType MapperType;

      //! mapper used to for block vector function 
      typedef typename Traits::BlockMapperType BlockMapperType;

      //! size of local blocks
      static const int localBlockSize = Traits::localBlockSize;

      //! type for DoF
      typedef RangeFieldType DofType;
      //! dimension of a value
      static const int dimVal = 1;

      //! mapper singleton key 
      //typedef LagrangeMapperSingletonKey< GridPartType, LagrangePointSetContainerType >
      typedef LagrangeMapperSingletonKey< GridPartType, LocalKeyStorageType >
        MapperSingletonKeyType;

      //! mapper factory 
      typedef LagrangeMapperSingletonFactory< MapperSingletonKeyType, BlockMapperType >
        BlockMapperSingletonFactoryType;

      //! singleton list of mappers 
      typedef SingletonList< MapperSingletonKeyType, BlockMapperType, BlockMapperSingletonFactoryType > BlockMapperProviderType;

    public:
      //! type of identifier for this discrete function space
      typedef int IdentifierType;
      //! identifier of this discrete function space
      static const IdentifierType id = 665;
      
    public:
      using BaseType::gridPart;
      using BaseType::order;

    public:
      //! default communication interface 
      static const InterfaceType defaultInterface = InteriorBorder_InteriorBorder_Interface;

      //! default communication direction 
      static const CommunicationDirection defaultDirection = ForwardCommunication;

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part for the Lagrange space
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       */
      explicit LagrangeDiscreteFunctionSpace
        ( GridPartType &gridPart,
          const InterfaceType commInterface = defaultInterface,
          const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        shapeFunctionSets_(),
        lagrangePointSetContainer_( gridPart ),
        mapper_( 0 ),
        blockMapper_( 0 )
      {
        const IndexSetType &indexSet = gridPart.indexSet();

        AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );
        const std::vector< GeometryType > &geometryTypes = allGeometryTypes.geomTypes( 0 );
        for( unsigned int i = 0; i < geometryTypes.size(); ++i )
        {
          const GeometryType &gt = geometryTypes[ i ];
          // insert shape function set for geometry type 
          shapeFunctionSets_.template insert< BaseFunctionSetSingletonProviderType >( gt );
        }

        MapperSingletonKeyType key( gridPart, lagrangePointSetContainer_.compiledLocalKeys( polynomialOrder ), polynomialOrder );
        blockMapper_ = &BlockMapperProviderType::getObject( key );
        assert( blockMapper_ != 0 );
        mapper_ = new MapperType( *blockMapper_ );
        assert( mapper_ != 0 );
      }

    private:
      // forbid the copy constructor
      LagrangeDiscreteFunctionSpace ( const ThisType& );

    public:
      /** \brief Destructor (freeing base functions and mapper)
          \return 
      **/
      ~LagrangeDiscreteFunctionSpace ()
      {
        delete mapper_;
        BlockMapperProviderType::removeObject( *blockMapper_ );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::contains */
      bool contains ( const int codim ) const
      {
        // forward to mapper since this information is held there 
        return blockMapper().contains( codim );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous () const
      {
        return (polynomialOrder > 0);
      }
      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      inline bool continuous (const IntersectionType &intersection) const
      { 
        return ( (polynomialOrder > 0) && intersection.conforming() );
      }

      /** \brief get the type of this discrete function space
       *  \return DFSpaceIdentifier
       */
      DFSpaceIdentifier type () const
      {
        return LagrangeSpace_id;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const
      {
        return polynomialOrder;
      }

      /** \brief provide access to the base function set for an entity 
       *
       *  \param[in]  entity   entity for which the base function set is requested 
       *
       *  \returns base function set for the specified entity
       */
      template< class EntityType >
      const BaseFunctionSetType baseFunctionSet ( const EntityType &entity ) const
      {
        return baseFunctionSet( entity.type() );
      }

      /** \brief provide access to the base function set for a geometry type
       *
       *  \param[in]  type  type of geometry the base function set is requested for
       *
       *  \returns base function set for the specified geometry
       */
      const BaseFunctionSetType baseFunctionSet ( const GeometryType type ) const
      {
        return BaseFunctionSetType( & shapeFunctionSets_[ type ] );
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
      const LagrangePointSetType &lagrangePointSet ( const GeometryType type ) const
      {
        return lagrangePointSetContainer_.compiledLocalKey( type, polynomialOrder );
      }

      /** \brief get dimension of value
          \return int
      **/
      int dimensionOfValue () const
      {
        return dimVal;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::mapper */
      DUNE_VERSION_DEPRECATED(1,4,remove)
      MapperType &mapper () const
      {
        assert( mapper_ != 0 );
        return *mapper_;
      }

      /** \brief obtain the DoF block mapper of this space
          \return BlockMapperType
      **/
      BlockMapperType &blockMapper () const
      {
        assert( blockMapper_ != 0 );
        return *blockMapper_;
      }

    protected:
      mutable ShapeSetStorageType shapeFunctionSets_;
      LagrangePointSetContainerType lagrangePointSetContainer_;
      MapperType *mapper_;
      BlockMapperType *blockMapper_;
    };
    
  } // namespace Fem 

#if DUNE_FEM_COMPATIBILITY  
// put this in next version 1.4 

  using Fem :: LagrangeDiscreteFunctionSpace ;
#endif // DUNE_FEM_COMPATIBILITY

} // namespace Dune

// include definition of RestrictProlongDefault for Lagrange Space.
#include <dune/fem/space/lagrangespace/adaptmanager.hh>

#endif // #ifndef DUNE_FEM_LAGRANGESPACE_LAGRANGESPACE_HH

#ifndef DUNE_LAGRANGESPACE_LAGRANGESPACE_HH
#define DUNE_LAGRANGESPACE_LAGRANGESPACE_HH

//- system includes 
#include <algorithm>

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>

//- Dune-Fem includes 
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>

//- local includes 
#include "basefunctions.hh"
#include "mapper.hh"
#include "lagrangedatahandle.hh"

namespace Dune
{

  template< class FunctionSpaceImp,
            class GridPartImp,
            int polOrder,
            template< class > class BaseFunctionStorageImp = CachingStorage >
  class LagrangeDiscreteFunctionSpace;



  template< class FunctionSpace, class GridPart, unsigned int polOrder,
            template< class > class BaseFunctionStorage = CachingStorage >
  struct LagrangeDiscreteFunctionSpaceTraits
  {
    CompileTimeChecker< (polOrder > 0) >
      __LagrangeSpace_only_defined_for_polOrder_greater_zero__;
    
    typedef FunctionSpace FunctionSpaceType;
    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType :: ScalarFunctionSpaceType
      ScalarFunctionSpaceType;
    enum { dimRange = FunctionSpaceType :: dimRange };
    
    typedef GridPart GridPartType;
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridPartType :: IndexSetType IndexSetType;
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType
      IteratorType;

    enum { polynomialOrder = polOrder };
    
    typedef LagrangeDiscreteFunctionSpace
      < FunctionSpaceType, GridPartType, polynomialOrder, BaseFunctionStorage >
      DiscreteFunctionSpaceType;
    typedef LagrangeMapper< GridPartType, polynomialOrder, dimRange >
      MapperType;
    
    // mapper for block 
    typedef LagrangeMapper< GridPartType, polynomialOrder, 1 >
      BlockMapperType;
    
    // implementation of basefunction set 
    typedef VectorialBaseFunctionSet< FunctionSpaceType, BaseFunctionStorage >
        BaseFunctionSetImp;

    // exported type 
    typedef SimpleBaseFunctionProxy<BaseFunctionSetImp>  BaseFunctionSetType;

    enum { localBlockSize = dimRange };

    /** \brief defines type of communication data handle for this type of space
     */
    template< class DiscreteFunction,
              class Operation = DFCommunicationOperation :: Add >
    struct CommDataHandle
    {
      //! type of data handle 
      typedef LagrangeCommunicationHandler< DiscreteFunction, Operation > Type;
      //! type of operatation to perform on scatter 
      typedef Operation OperationType;
    };
  };

  //! Key for Mapper singleton list 
  template< class GridPartImp, class LagrangePointSetMapImp >
  class LagrangeMapperSingletonKey 
  {
    const GridPartImp & gridPart_; 
    mutable LagrangePointSetMapImp& pointSet_;
    const int polOrd_;
  public:
    //! constructor taking index set and numDofs 
    LagrangeMapperSingletonKey(const GridPartImp & gridPart, 
                               LagrangePointSetMapImp& pointSet,
                               const int polOrd)
      : gridPart_(gridPart) ,  pointSet_(pointSet) , polOrd_(polOrd) 
    {}
    //! copy constructor 
    LagrangeMapperSingletonKey(const LagrangeMapperSingletonKey &org) 
      : gridPart_(org.gridPart_) , pointSet_(org.pointSet_), polOrd_(org.polOrd_)
    {}
    //! returns true if indexSet pointer and numDofs are equal 
    bool operator == (const LagrangeMapperSingletonKey & otherKey) const 
    {
      return ((&(gridPart_.indexSet()) == &(otherKey.gridPart().indexSet())) 
              && (polOrd_ == otherKey.polOrd_));
    }

    //! return reference to index set 
    const GridPartImp & gridPart() const { return gridPart_; }
    //! return lagrange point map set  
    LagrangePointSetMapImp& pointSet() const { return pointSet_; }

  };


  
  // Factory class for SingletonList to tell how objects are created and
  // how compared.
  template< class Key, class Object >
  struct LagrangeMapperSingletonFactory
  {
    static Object *createObject( const Key &key )
    {
      return new Object( key.gridPart(), key.pointSet() );
    }
    
    static void deleteObject( Object *obj )
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



  /** \class   LagrangeDiscreteFunctionSpace
   *  \ingroup LagrangeDiscreteFunctionSpace
   *  \brief   Lagrange discrete function space
   */
  template< class FunctionSpaceImp,
            class GridPartImp,
            int polOrder,
            template< class > class BaseFunctionStorageImp >
  class LagrangeDiscreteFunctionSpace
  : public DiscreteFunctionSpaceDefault
           < LagrangeDiscreteFunctionSpaceTraits< FunctionSpaceImp,
                                                  GridPartImp,
                                                  polOrder,
                                                  BaseFunctionStorageImp > >
  {
  public:
    //! traits for the discrete function space
    typedef LagrangeDiscreteFunctionSpaceTraits< FunctionSpaceImp,
                                                 GridPartImp,
                                                 polOrder,
                                                 BaseFunctionStorageImp >
      Traits;

    //! type of the discrete function space
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceImp,
                                           GridPartImp,
                                           polOrder,
                                           BaseFunctionStorageImp >
     LagrangeDiscreteFunctionSpaceType;

    typedef typename Traits :: GridPartType GridPartType;
    typedef typename Traits :: GridType GridType;
    typedef typename Traits :: IndexSetType IndexSetType;
    typedef typename Traits :: IteratorType IteratorType;
    //! dimension of the grid (not the world)
    enum { dimension = GridType :: dimension };

    typedef typename Traits :: FunctionSpaceType FunctionSpaceType;
    //! field type for function space's domain
    typedef typename Traits :: DomainFieldType DomainFieldType;
    //! type for function space's domain
    typedef typename Traits :: DomainType DomainType;
    //! field type for function space's range
    typedef typename Traits :: RangeFieldType RangeFieldType;
    //! type for function space's range
    typedef typename Traits :: RangeType RangeType;
    //! dimension of function space's range
    enum { dimRange = FunctionSpaceType :: dimRange };
    //! type of scalar function space
    typedef typename Traits :: ScalarFunctionSpaceType ScalarFunctionSpaceType;
   
    //! maximum polynomial order of functions in this space
    enum { polynomialOrder = Traits :: polynomialOrder };
    
    //! type of the base function set(s)
    typedef typename Traits :: BaseFunctionSetImp BaseFunctionSetImp;
    
    typedef typename Traits :: BaseFunctionSetType BaseFunctionSetType;
    //! type of the base function set map
    typedef std :: map< const GeometryType, const BaseFunctionSetImp* >
      BaseFunctionMapType;
    //! type of base function factory
    typedef LagrangeBaseFunctionFactory
      < ScalarFunctionSpaceType, dimension, polynomialOrder >
      ScalarFactoryType;
    //! type of singleton base function factory
    typedef BaseFunctionSetSingletonFactory
      < GeometryType, BaseFunctionSetImp, ScalarFactoryType >
      BaseFunctionSetSingletonFactoryType;
    //! type of singleton list (singleton provider) for base functions
    typedef SingletonList
      < GeometryType, BaseFunctionSetImp, BaseFunctionSetSingletonFactoryType >
      BaseFunctionSetSingletonProviderType;

    //! type of a Lagrange point set
    typedef LagrangePointSet< GridPartType, polynomialOrder >
      LagrangePointSetType;
    //! type of Lagrange point set map
    typedef std :: map< const GeometryType, const LagrangePointSetType* >
      LagrangePointSetMapType;

    //! mapper used to implement mapToGlobal
    typedef typename Traits :: MapperType MapperType;

    //! mapper used to for block vector function 
    typedef typename Traits :: BlockMapperType BlockMapperType;

    //! size of local blocks
    enum { localBlockSize = Traits :: localBlockSize };

    //! type for DoF
    typedef RangeFieldType DofType;
    //! dimension of a value
    enum { dimVal = 1 };
    //! type of DoF manager
    typedef DofManager< GridType > DofManagerType;
    //! type of DoF manager factory
    typedef DofManagerFactory< DofManagerType > DofManagerFactoryType;

    //! mapper singleton key 
    typedef LagrangeMapperSingletonKey< GridPartType, LagrangePointSetMapType >
      MapperSingletonKeyType;

    //! mapper factory 
    typedef LagrangeMapperSingletonFactory< MapperSingletonKeyType, MapperType >
      MapperSingletonFactoryType;

    //! singleton list of mappers 
    typedef SingletonList
      < MapperSingletonKeyType, MapperType, MapperSingletonFactoryType >
      MapperProviderType;

    //! mapper factory 
    typedef LagrangeMapperSingletonFactory
      < MapperSingletonKeyType, BlockMapperType >
      BlockMapperSingletonFactoryType;

    //! singleton list of mappers 
    typedef SingletonList
      < MapperSingletonKeyType, BlockMapperType, BlockMapperSingletonFactoryType >
      BlockMapperProviderType;

  public:
    //! type of identifier for this discrete function space
    typedef int IdentifierType;
    //! identifier of this discrete function space
    static const IdentifierType id = 665;
    
  private:
    typedef LagrangeDiscreteFunctionSpaceType ThisType;
    typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

  private:
    //! map for the base function sets
    mutable BaseFunctionMapType baseFunctionSet_;

    //! map for the langrage point sets
    mutable LagrangePointSetMapType lagrangePointSet_;
    
    //! corresponding mapper
    MapperType *mapper_;

    //! corresponding mapper
    BlockMapperType *blockMapper_;

  public:
    using BaseType :: gridPart;

  public:
    /** \brief Constructor generating a LagrangeBaseFunctionSet of the requested polynomial order for each element type of the grid 
        \param[in] gridPart 
        \return 
    **/
    inline explicit LagrangeDiscreteFunctionSpace ( GridPartType &gridPart )
    : BaseType( gridPart ),
      baseFunctionSet_(),
      lagrangePointSet_(),
      mapper_( 0 ),
      blockMapper_( 0 )
    {
      const IndexSetType &indexSet = gridPart.indexSet();

      AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );
      const std :: vector< GeometryType >& geometryTypes
        = allGeometryTypes.geomTypes( 0 );
      for( unsigned int i = 0; i < geometryTypes.size(); ++i )
      {
        const GeometryType &geometryType = geometryTypes[ i ];
        
        if( baseFunctionSet_.find( geometryType ) == baseFunctionSet_.end() )
        {
          const BaseFunctionSetImp *baseFunctionSet
            = &(BaseFunctionSetSingletonProviderType
                :: getObject( geometryType ));
          assert( baseFunctionSet != NULL );
          baseFunctionSet_[ geometryType ] = baseFunctionSet;
        }

        if( lagrangePointSet_.find( geometryType ) == lagrangePointSet_.end() )
        {
          const LagrangePointSetType *lagrangePointSet
            = new LagrangePointSetType( geometryType );
          assert( lagrangePointSet != NULL );
          lagrangePointSet_[ geometryType ] = lagrangePointSet;
        }
      }

      MapperSingletonKeyType key( gridPart, lagrangePointSet_, polynomialOrder );
      mapper_ = &MapperProviderType :: getObject(key);
      assert( mapper_ != 0 );
      blockMapper_ = &BlockMapperProviderType :: getObject( key );
      assert( blockMapper_ != 0 );
    }

  private:
    // forbid the copy constructor
    LagrangeDiscreteFunctionSpace ( const ThisType& );

  public:
    /** \brief Destructor (freeing base functions and mapper)
        \return 
    **/
    inline ~LagrangeDiscreteFunctionSpace ()
    {
      if( blockMapper_ ) BlockMapperProviderType::removeObject( *blockMapper_ );
      MapperProviderType::removeObject( *mapper_ );

      typedef typename BaseFunctionMapType :: iterator BFIteratorType;
      BFIteratorType bfend = baseFunctionSet_.end();
      for( BFIteratorType it = baseFunctionSet_.begin(); it != bfend; ++it ) 
      {
        const BaseFunctionSetImp *baseFunctionSet = (*it).second;
        if( baseFunctionSet != NULL )
          BaseFunctionSetSingletonProviderType
          :: removeObject( *baseFunctionSet );
      }

      typedef typename LagrangePointSetMapType :: iterator LPIteratorType;
      LPIteratorType lpend = lagrangePointSet_.end();
      for( LPIteratorType it = lagrangePointSet_.begin(); it != lpend; ++it ) 
      {
        const LagrangePointSetType *lagrangePointSet = (*it).second;
        if( lagrangePointSet != NULL )
          delete lagrangePointSet;
      }
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::contains */
    inline bool contains ( const int codim ) const
    {
      return mapper().contains( codim );
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::continuous */
    inline bool continuous () const
    {
      return (polynomialOrder > 0);
    }

    /** \brief get the type of this discrete function space 
        \return DFSpaceIdentifier
    **/
    inline DFSpaceIdentifier type () const
    {
      return LagrangeSpace_id;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::order */
    inline int order () const
    {
      return polynomialOrder;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::baseFunctionSet(const EntityType &entity) const */
    template< class EntityType >
    inline const BaseFunctionSetType baseFunctionSet ( const EntityType &entity ) const
    {
      return baseFunctionSet( entity.geometry().type() );
    }

    /** \brief provide access to the base function set for a geometry type
     *
     *  \param[in]  type  type of geometry the base function set is requested for
     *
     *  \returns base function set for the specified geometry
     */
    inline const BaseFunctionSetType baseFunctionSet ( const GeometryType type ) const
    {
      assert( baseFunctionSet_.find( type ) != baseFunctionSet_.end() );
      assert( baseFunctionSet_[ type ] != NULL );
      return BaseFunctionSetType( baseFunctionSet_[ type ] );
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
    inline const LagrangePointSetType &lagrangePointSet ( const EntityType &entity ) const
    {
      return this->lagrangePointSet( entity.geometry().type() );
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
    inline const LagrangePointSetType &lagrangePointSet ( const GeometryType type ) const
    {
      assert( lagrangePointSet_.find( type ) != lagrangePointSet_.end() );
      assert( lagrangePointSet_[ type ] != NULL );
      return *lagrangePointSet_[ type ];
    }

    /** \brief get dimension of value
        \return int
    **/
    inline int dimensionOfValue () const
    {
      return dimVal;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::mapper */
    inline MapperType &mapper () const
    {
      assert( mapper_ != 0 );
      return *mapper_;
    }

    /** \brief obtain the DoF block mapper of this space
        \return BlockMapperType
    **/
    inline BlockMapperType &blockMapper () const
    {
      assert( blockMapper_ != 0 );
      return *blockMapper_;
    }
  };
  
} // end Dune namespace  

#endif

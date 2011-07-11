#ifndef DUNE_PADATPTIVELAGRANGESPACE_LAGRANGESPACE_HH
#define DUNE_PADATPTIVELAGRANGESPACE_LAGRANGESPACE_HH

//- system includes 
#include <algorithm>

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>

//- Dune-Fem includes 
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>

//- local includes 
#include <dune/fem/space/lagrangespace/basefunctions.hh>
#include "mapper.hh"

namespace Dune
{

  template< class FunctionSpaceImp,
            class GridPartImp,
            int polOrder,
            template< class > class BaseFunctionStorageImp = CachingStorage >
  class PAdaptiveLagrangeSpace;



  template< class FunctionSpace, class GridPart, unsigned int polOrder,
            template< class > class BaseFunctionStorage = CachingStorage >
  struct PAdaptiveLagrangeSpaceTraits
  {
    dune_static_assert((polOrder > 0), "LagrangeSpace only defined for polOrder > 0" );
    
    typedef FunctionSpace FunctionSpaceType;
    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    enum { dimRange = FunctionSpaceType :: dimRange };
    
    typedef GridPart GridPartType;
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridPartType :: IndexSetType IndexSetType;
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType
      IteratorType;

    // get dimension of local coordinate 
    enum { dimLocal = GridType :: dimension };

    typedef typename ToLocalFunctionSpace< FunctionSpaceType, dimLocal > :: Type 
      BaseFunctionSpaceType;

    enum { polynomialOrder = polOrder };
    
    typedef PAdaptiveLagrangeSpace
      < FunctionSpaceType, GridPartType, polynomialOrder, BaseFunctionStorage >
      DiscreteFunctionSpaceType;

    enum { localBlockSize = dimRange };

    // mapper for block
    typedef LagrangeMapper< GridPartType, polynomialOrder > BlockMapperType;
    typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;
    
    // implementation of basefunction set 
    typedef VectorialBaseFunctionSet< BaseFunctionSpaceType, BaseFunctionStorage >
        BaseFunctionSetImp;

    // exported type 
    typedef SimpleBaseFunctionProxy<BaseFunctionSetImp>  BaseFunctionSetType;

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

  //! Key for Mapper singleton list 
  template< class GridPartImp, class LagrangePointSetMapImp >
  class PAdaptiveMapperSingletonKey 
  {
    const GridPartImp & gridPart_; 
    mutable LagrangePointSetMapImp& pointSet_;
    const int polOrd_;
  public:
    //! constructor taking index set and numDofs 
    PAdaptiveMapperSingletonKey(const GridPartImp & gridPart, 
                               LagrangePointSetMapImp& pointSet,
                               const int polOrd)
      : gridPart_(gridPart) ,  pointSet_(pointSet) , polOrd_(polOrd) 
    {}
    //! copy constructor 
    PAdaptiveMapperSingletonKey(const PAdaptiveMapperSingletonKey &org) 
      : gridPart_(org.gridPart_) , pointSet_(org.pointSet_), polOrd_(org.polOrd_)
    {}
    //! returns true if indexSet pointer and numDofs are equal 
    bool operator == (const PAdaptiveMapperSingletonKey & otherKey) const 
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
  struct PAdaptiveMapperSingletonFactory
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



  /** \addtogroup PAdaptiveLagrangeSpace
   *
   *  Provides access to bse function sets for different element types in
   *  one grid and size of function space and maps from local to global dof
   *  number.
   *
   *  \note This space can only be used with special index sets. If you want
   *  to use the PAdaptiveLagrangeSpace with an index set only
   *  supporting the index set interface you will have to use the
   *  IndexSetWrapper class to provide the required functionality.
   *
   *  \note For adaptive calculations one has to use index sets that are
   *  capable of adaption (i.e. the method adaptive returns true). See also
   *  AdaptiveLeafIndexSet.
   */



  /** \class   PAdaptiveLagrangeSpace
   *  \ingroup PAdaptiveLagrangeSpace
   *  \brief   Lagrange discrete function space
   */
  template< class FunctionSpaceImp,
            class GridPartImp,
            int polOrder,
            template< class > class BaseFunctionStorageImp >
  class PAdaptiveLagrangeSpace
  : public DiscreteFunctionSpaceDefault
           < PAdaptiveLagrangeSpaceTraits< FunctionSpaceImp,
                                                  GridPartImp,
                                                  polOrder,
                                                  BaseFunctionStorageImp > >
  {
  public:
    //! traits for the discrete function space
    typedef PAdaptiveLagrangeSpaceTraits< FunctionSpaceImp,
                                                 GridPartImp,
                                                 polOrder,
                                                 BaseFunctionStorageImp >
      Traits;

    //! type of the discrete function space
    typedef PAdaptiveLagrangeSpace< FunctionSpaceImp,
                                           GridPartImp,
                                           polOrder,
                                           BaseFunctionStorageImp >
     PAdaptiveLagrangeSpaceType;

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
    typedef typename Traits :: BaseFunctionSpaceType BaseFunctionSpaceType;
   
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
      < typename BaseFunctionSpaceType :: ScalarFunctionSpaceType, dimension, polynomialOrder >
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

    //! mapper singleton key 
    typedef PAdaptiveMapperSingletonKey< GridPartType, LagrangePointSetMapType >
      MapperSingletonKeyType;

    //! mapper factory 
    typedef PAdaptiveMapperSingletonFactory< MapperSingletonKeyType, MapperType >
      MapperSingletonFactoryType;

    //! mapper factory 
    typedef PAdaptiveMapperSingletonFactory
      < MapperSingletonKeyType, BlockMapperType >
      BlockMapperSingletonFactoryType;

    //! singleton list of mappers 
    typedef SingletonList
      < MapperSingletonKeyType, BlockMapperType, BlockMapperSingletonFactoryType >
      BlockMapperProviderType;

    typedef std::vector<BaseFunctionMapType> BaseFunctionMapVectorType;
    typedef std::vector<LagrangePointSetMapType> LagrangePointSetMapVectorType;

    template <int pOrd> 
    struct ConstructBaseFunctionSets
    {
      static void apply( BaseFunctionMapVectorType& baseFunctionSetVector, 
                         LagrangePointSetMapVectorType& lagrangePointSetVector,
                         const GeometryType& geometryType ) 
      {
        //std::cout << "Create baseFunctionSet of order " << pOrd << std::endl;
        typedef LagrangeBaseFunctionFactory
          < typename BaseFunctionSpaceType :: ScalarFunctionSpaceType, dimension, pOrd >
          ScalarFactoryType;

        //! type of singleton base function factory
        typedef BaseFunctionSetSingletonFactory
          < GeometryType, BaseFunctionSetImp, ScalarFactoryType >
          BaseFunctionSetSingletonFactoryType;

        //! type of singleton list (singleton provider) for base functions
        typedef SingletonList
          < GeometryType, BaseFunctionSetImp, BaseFunctionSetSingletonFactoryType >
        BaseFunctionSetSingletonProviderType;

        const size_t k = pOrd;
        assert( k < baseFunctionSetVector.size() );
        BaseFunctionMapType& baseFunctionMap = baseFunctionSetVector[ k ];


        if( baseFunctionMap.find( geometryType ) == baseFunctionMap.end() )
        {
          const BaseFunctionSetImp *baseFunctionSet
            = &(BaseFunctionSetSingletonProviderType
               :: getObject( geometryType ));
          assert( baseFunctionSet != NULL );
          baseFunctionMap[ geometryType ] = baseFunctionSet;
        }

        LagrangePointSetMapType& lpSet = lagrangePointSetVector[ k ];
        if( lpSet.find( geometryType ) == lpSet.end() )
        {
          const LagrangePointSetType *lagrangePointSet
            = new LagrangePointSetType( geometryType, k );
          assert( lagrangePointSet != NULL );
          lpSet[ geometryType ] = lagrangePointSet;
        }
      }
    };

  public:
    //! type of identifier for this discrete function space
    typedef int IdentifierType;
    //! identifier of this discrete function space
    static const IdentifierType id = 665;
    
  private:
    typedef PAdaptiveLagrangeSpaceType ThisType;
    typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

  private:
    //! map for the base function sets
    mutable BaseFunctionMapVectorType baseFunctionSet_;

    //! map for the langrage point sets
    mutable LagrangePointSetMapVectorType lagrangePointSet_;
    
    //! corresponding mapper
    MapperType *mapper_;

    //! corresponding mapper
    BlockMapperType *blockMapper_;

  public:
    using BaseType :: gridPart;

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
    explicit PAdaptiveLagrangeSpace
      ( GridPartType &gridPart,
        const InterfaceType commInterface = defaultInterface,
        const CommunicationDirection commDirection = defaultDirection )
    : BaseType( gridPart, commInterface, commDirection ),
      baseFunctionSet_( polynomialOrder+1 ),
      lagrangePointSet_( polynomialOrder+1 ),
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

        ForLoop< ConstructBaseFunctionSets, 1, polynomialOrder > :: 
          apply( baseFunctionSet_, lagrangePointSet_, geometryType );
      }

      MapperSingletonKeyType key( gridPart, lagrangePointSet_[ polynomialOrder ], polynomialOrder );
      blockMapper_ = &BlockMapperProviderType :: getObject( key );
      assert( blockMapper_ != 0 );
      mapper_ = new MapperType( *blockMapper_ );
      assert( mapper_ != 0 );
    }

  private:
    // forbid the copy constructor
    PAdaptiveLagrangeSpace ( const ThisType& );

  public:
    /** \brief Destructor (freeing base functions and mapper)
        \return 
    **/
    ~PAdaptiveLagrangeSpace ()
    {
      delete mapper_;
      BlockMapperProviderType::removeObject( *blockMapper_ );

      for( size_t i = 0; i < baseFunctionSet_.size(); ++ i) 
      {
        typedef typename BaseFunctionMapType :: iterator BFIteratorType;
        BFIteratorType bfend = baseFunctionSet_[ i ].end();
        for( BFIteratorType it = baseFunctionSet_[ i ].begin(); it != bfend; ++it ) 
        {
          const BaseFunctionSetImp *baseFunctionSet = (*it).second;
          if( baseFunctionSet != NULL )
            BaseFunctionSetSingletonProviderType
            :: removeObject( *baseFunctionSet );
        }
      }

      for( size_t i = 0; i < lagrangePointSet_.size(); ++ i) 
      {
        typedef typename LagrangePointSetMapType :: iterator LPIteratorType;
        const LPIteratorType lpend = lagrangePointSet_[ i ].end();
        for( LPIteratorType it = lagrangePointSet_[ i ].begin(); it != lpend; ++it ) 
        {
          const LagrangePointSetType *lagrangePointSet = (*it).second;
          if( lagrangePointSet != NULL )
            delete lagrangePointSet;
        }
      }
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::contains */
    inline bool contains ( const int codim ) const
    {
      // forward to mapper since this information is held there 
      return blockMapper().contains( codim );
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
      return baseFunctionSet( entity.type(), 
                              blockMapper_->polynomOrder( entity ) );
    }

    /** \brief provide access to the base function set for a geometry type
     *
     *  \param[in]  type  type of geometry the base function set is requested for
     *
     *  \returns base function set for the specified geometry
     */
    inline const BaseFunctionSetType baseFunctionSet ( const GeometryType type ) const
    {
      assert( false );
      //assert( baseFunctionSet_.find( type ) != baseFunctionSet_.end() );
      //assert( baseFunctionSet_[ type ] != NULL );
      return BaseFunctionSetType( baseFunctionSet_[ type ], 1 );
    }

    inline const BaseFunctionSetType baseFunctionSet ( const GeometryType type, const int k ) const
    {
      assert( k <= polynomialOrder );
      assert( k > 0 );
      //assert( baseFunctionSet_[ k ].find( type ) != baseFunctionSet_[ k ].end() );
      //assert( baseFunctionSet_[ k ][ type ] != NULL );
      return BaseFunctionSetType( baseFunctionSet_[ k ][ type ] );
    }

    /** \brief provide access to the Lagrange point set for an entity
     *
     *  \note This method is not part of the DiscreteFunctionSpaceInterface. It
     *        is unique to the PAdaptiveLagrangeSpace.
     *
     *  \param[in]  entity  entity the Lagrange point set is requested for
     *  
     *  \returns LagrangePointSet
     */
    template< class EntityType >
    inline const LagrangePointSetType &lagrangePointSet ( const EntityType &entity ) const
    {
      return this->lagrangePointSet( entity.type(),
                                     blockMapper_->polynomOrder( entity ) );
    }

    /** \brief provide access to the Lagrange point set for a geometry type
     *
     *  \note This method is not part of the DiscreteFunctionSpaceInterface. It
     *        is unique to the PAdaptiveLagrangeSpace.
     *
     *  \param[in]  type  type of geometry the Lagrange point set is requested for
     *
     *  \returns LagrangePointSetType
     */
    inline const LagrangePointSetType &lagrangePointSet ( const GeometryType type ) const
    {
      LagrangePointSetMapType& lagrangePointSet = lagrangePointSet_[ polynomialOrder ];
      abort();
      assert( lagrangePointSet.find( type ) != lagrangePointSet.end() );
      assert( lagrangePointSet[ type ] != NULL );
      return *lagrangePointSet[ type ];
    }

    /** \brief provide access to the Lagrange point set for a geometry type
     *
     *  \note This method is not part of the DiscreteFunctionSpaceInterface. It
     *        is unique to the PAdaptiveLagrangeSpace.
     *
     *  \param[in]  type  type of geometry the Lagrange point set is requested for
     *  \param[in]  order polynomial order for given geometry type 
     *
     *  \returns LagrangePointSetType
     */
    inline const LagrangePointSetType &lagrangePointSet ( const GeometryType type, const int order ) const
    {
      LagrangePointSetMapType& lagrangePointSet = lagrangePointSet_[ order ];
      assert( lagrangePointSet.find( type ) != lagrangePointSet.end() );
      assert( lagrangePointSet[ type ] != NULL );
      return *lagrangePointSet[ type ];
    }

    /** \brief get dimension of value
        \return int
    **/
    inline int dimensionOfValue () const
    {
      return dimVal;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::mapper */
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
  };
  
} // end Dune namespace  

// include definition of RestrictProlongDefault for Lagrange Space.
#include "adaptmanager.hh"

#endif // #ifndef DUNE_LAGRANGESPACE_LAGRANGESPACE_HH

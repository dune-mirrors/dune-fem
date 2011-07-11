#ifndef DUNE_PADAPTIVELAGRANGESPACE_MAPPER_HH
#define DUNE_PADAPTIVELAGRANGESPACE_MAPPER_HH

//- Dune includes 
#include <dune/common/geometrytype.hh>
#include <dune/common/exceptions.hh>

//- Dune-Fem includes 
#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/misc/codimmap.hh>
#include <dune/fem/misc/metaprogramming.hh>
#include <dune/fem/space/common/dofmanager.hh>

#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>

//- local includes 
#include <dune/fem/space/lagrangespace/lagrangepoints.hh>

#include <dune/grid/utility/persistentcontainer.hh>

namespace Dune
{
  
  template< class GridPart, int polOrder >
  class PAdaptiveLagrangeMapper;


  
  template< class GridPart, int polOrder >
  struct PAdaptiveLagrangeMapperTraits
  {
    typedef GridPart GridPartType;
    
    static const int polynomialOrder = polOrder;

    typedef typename GridPartType::template Codim< 0 >::IteratorType::Entity EntityType;
    typedef PAdaptiveLagrangeMapper< GridPartType, polynomialOrder > DofMapperType;
    typedef DefaultDofMapIterator< EntityType, DofMapperType > DofMapIteratorType;
  };



  // First Order Lagrange Mapper
  // ---------------------------

  template< class GridPart >
  class PAdaptiveLagrangeMapper< GridPart, 1 >
  : public CodimensionMapper< GridPart, GridPart::GridType::dimension >
  {
    typedef PAdaptiveLagrangeMapper< GridPart, 1 > ThisType;
    typedef CodimensionMapper< GridPart, GridPart::GridType::dimension > BaseType;

  public:
    //! type of the grid part
    typedef typename BaseType::GridPartType GridPartType;

    //! type of entities (codim 0)
    typedef typename BaseType::EntityType EntityType;

    //! type of the underlying grid
    typedef typename GridPartType::GridType GridType;

    //! type of coordinates within the grid
    typedef typename GridType::ctype FieldType;

    //! dimension of the grid
    static const int dimension = GridType::dimension;

    //! order of the Lagrange polynoms
    static const int polynomialOrder = 1; 

    //! type of the Lagrange point set
    typedef LagrangePointSet< GridPartType, polynomialOrder >
      LagrangePointSetType;
    //! type of the map for the Lagrange point sets
    typedef std::map< const GeometryType, const LagrangePointSetType* >
      LagrangePointSetMapType;

    typedef std::vector< LagrangePointSetMapType > LagrangePointSetMapVectorType;

  public:
    //! constructor
    PAdaptiveLagrangeMapper ( const GridPartType &gridPart, LagrangePointSetMapVectorType &lagrangePointSet )
    : BaseType( gridPart )
    {}

    bool fixedDataSize ( const int codim ) const
    {
      return true;
    }
  };



  // Second Order Lagrange Mapper
  // ----------------------------

  template< class GridPart >
  class PAdaptiveLagrangeMapper< GridPart, 2 >
  : public DofMapperDefault< PAdaptiveLagrangeMapperTraits< GridPart, 2 > >
  {
    typedef PAdaptiveLagrangeMapper< GridPart, 2 > ThisType;
    typedef DofMapperDefault< PAdaptiveLagrangeMapperTraits< GridPart, 2 > > BaseType;

  public:
    typedef PAdaptiveLagrangeMapperTraits< GridPart, 2 > Traits;
    
    //! type of the grid part
    typedef typename Traits::GridPartType GridPartType;

    //! type of entities (codim 0)
    typedef typename Traits::EntityType EntityType;

    //! type of DofMapIterator
    typedef typename Traits::DofMapIteratorType DofMapIteratorType;
 
    //! type of the underlying grid
    typedef typename GridPartType::GridType GridType;

    //! type of the index set
    typedef typename GridPartType::IndexSetType IndexSetType;

    //! type of coordinates within the grid
    typedef typename GridType::ctype FieldType;

    //! dimension of the grid
    static const int dimension = GridType::dimension;

    //! order of the Lagrange polynoms
    static const int polynomialOrder = Traits::polynomialOrder;

    //! type of the Lagrange point set
    typedef LagrangePointSet< GridPartType, polynomialOrder >
      LagrangePointSetType;
    //! type of the map for the Lagrange point sets
    typedef std::map< const GeometryType, const LagrangePointSetType* >
      LagrangePointSetMapType;

    typedef std::vector< LagrangePointSetMapType > LagrangePointSetMapVectorType;

    //! type of the DoF manager
    typedef DofManager< GridType > DofManagerType;

    struct EntityDofStorage
    {
      std::vector< std::vector< int > > dofs_;

      EntityDofStorage() : 
        dofs_( polynomialOrder+1 ) 
      {}

      bool exists( const int polOrd ) const 
      {
        return dofs_[ polOrd ].size() > 0 ;
      }

      void insert( const int polOrd, const int numDofs, const int startDof ) 
      {
        assert( ! exists ( polOrd ) );
        {
          dofs_[ polOrd ].resize( numDofs );
          for(int i=0, dof=startDof ; i<numDofs; ++i, ++dof ) 
            dofs_[ polOrd ][ i ] = dof;
        }
      }

      int dof ( const int polOrd, const int dofNumber ) const 
      { 
        return dofs_[ polOrd ][ dofNumber ];
      }

    };

    typedef EntityDofStorage EntityDofStorageType;

    struct PolynomOrderStorage 
    {
      unsigned char k_;
      PolynomOrderStorage() : k_( polynomialOrder ) {}
      PolynomOrderStorage( const int k ) : k_( k ) {}
      int order () const { return k_;}
    };

    typedef PolynomOrderStorage  PolynomOrderStorageType;

    typedef PersistentContainer< GridType, EntityDofStorageType > DofContainerType ;

    typedef PersistentContainer< GridType, PolynomOrderStorageType > PolyOrderContainerType ;
  protected:
    typedef typename LagrangePointSetType::DofInfo DofInfo;

    template < int codim > 
    struct InsertSubEntities
    {
      static void apply( const EntityType& entity, 
                         const LagrangePointSetType& set,
                         const int polOrd,
                         unsigned int&  dofCounter,
                         std::vector< DofContainerType* > dofContainers ) 
      {
        DofContainerType& dofContainer = *dofContainers[ codim ];
        const int count = entity.template count< codim > ();
        for(int i=0; i<count; ++i ) 
        {
          EntityDofStorage& entityDofs = dofContainer( entity, i );
          if( ! entityDofs.exists( polOrd ) ) 
          {
            const int numDofs = set.numDofs( codim, i );
            entityDofs.insert( polOrd, numDofs, dofCounter );
            dofCounter += numDofs;
          }
        }
      }
    };

  public:
    //! constructor
    PAdaptiveLagrangeMapper ( const GridPartType &gridPart,
                     LagrangePointSetMapVectorType &lagrangePointSetVector )
    : gridPart_( gridPart ),
      dm_( DofManagerType :: instance(gridPart.grid()) ),
      indexSet_( gridPart.indexSet() ),
      lagrangePointSet_( lagrangePointSetVector ),
      entityPolynomOrder_( gridPart.grid(), 0 ),
      dofContainer_( dimension+1, (DofContainerType *) 0 ),
      overShoot_( Parameter::getValidValue( "fem.lagrangemapper.overshoot", double( 1.5 ), ValidateNotLess< double >( 1.0 ) ) )
    {
      numDofs_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
        maxDofs_[ codim ] = 0;
      for( int codim = 0; codim <= dimension; ++codim )
        dofContainer_[ codim ] = new DofContainerType( gridPart.grid(), codim );

      for( size_t i=0; i<lagrangePointSet_.size(); ++i ) 
      {
        typedef typename LagrangePointSetMapType :: iterator IteratorType;
        IteratorType end = lagrangePointSet_[ i ].end();
        for( IteratorType it = lagrangePointSet_[ i ].begin(); it != end; ++it )
        {
          const LagrangePointSetType *set = (*it).second;
          if( set == 0 )
            continue;
          
          numDofs_ = std :: max( numDofs_, set->size() );
          for( int codim = 0; codim <= dimension; ++codim )
            maxDofs_[ codim ] = std::max( maxDofs_[ codim ], set->maxDofs( codim ) );
        }
      }

      resize();
      // register to DofManager for adaptation 
      dm_.addIndexSet( *this );
    }

    int polynomOrder( const EntityType& entity ) const 
    {
      return entityPolynomOrder_[ entity ].order();
    }

    DofContainerType& dofContainer( const size_t codim ) const 
    {
      assert( codim < dofContainer_.size() );
      return *(dofContainer_[ codim ]);
    }

    const LagrangePointSetType* 
    lagrangePointSet( const int polOrd, const GeometryType type ) const 
    {
      // add polOrd here 
      return lagrangePointSet_[ polOrd ][ type ];
    }
    
    //! destructor 
    virtual ~PAdaptiveLagrangeMapper ()
    {
      dm_.removeIndexSet( *this );
    }

    //! return overall number of degrees of freedom 
    int size () const
    {
      return size_;
    }

    /** \copydoc Dune::DofMapper::begin(const EntityType &entity) const */
    DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::beginIterator, entity, *this );
    }
    
    /** \copydoc Dune::DofMapper::end(const EntityType &entity) const */
    DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::endIterator, entity, *this );
    }

    /** \copydoc Dune::DofMapper::mapToGlobal */
    int mapToGlobal ( const EntityType &entity, const int localDof ) const
    {
      const int polOrd = polynomOrder( entity );
      const LagrangePointSetType *set = lagrangePointSet( polOrd, entity.type() );

      const DofInfo &dofInfo = set->dofInfo( localDof );

      const unsigned int codim = dofInfo.codim;
      const unsigned int subEntity = dofInfo.subEntity;

      return dofContainer( codim )( entity, subEntity ).dof( polOrd, dofInfo.dofNumber );
    }

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const 
    {
      static const unsigned int codim = Entity::codimension;
      assert( (localDof >= 0) && (localDof < numEntityDofs( entity )) );
      return offset_[ codim ] + indexSet_.index( entity );
    }
    
    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    int maxNumDofs () const
    {
      return numDofs_;
    }

    /** \copydoc Dune::DofMapper::numDofs(const EntityType &entity) const */
    int numDofs ( const EntityType &entity ) const
    {
      const int polOrd = polynomOrder( entity );
      return lagrangePointSet( polOrd, entity.type() )->size();
    }

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      // This implementation only works for nonhybrid grids (simplices or cubes)
      return maxDofs_[ Entity::codimension ];
    }
    
    /** \brief Check, whether any DoFs are associated with a codimension */
    bool contains ( int codim ) const
    {
      return true;
    }

    /** \brief Check, whether the data in a codimension has fixed size */
    bool fixedDataSize ( int codim ) const
    {
      return false;
    }

    /** \copydoc Dune::DofMapper::oldIndex */
    int oldIndex ( const int hole, const int codim ) const
    {
      return offset_[ codim ] + indexSet_.oldIndex( hole, codim );
    }

    /** \copydoc Dune::DofMapper::newIndex */
    int newIndex ( const int hole, const int codim ) const
    {
      return offset_[ codim ] + indexSet_.newIndex( hole, codim );
    }

    /** \copydoc Dune::DofMapper::numberOfHoles */
    int numberOfHoles ( const int codim ) const
    {
      return 0;
    }

    /** \copydoc Dune::DofMapper::numBlocks
     */
    int numBlocks () const 
    {
      return 1;
    }

    /** \copydoc Dune::DofMapper::oldOffset
     */
    int oldOffSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return 0;//oldOffSet_[ block ];
    }

    /** \copydoc Dune::DofMapper::newOffset
     */
    int offSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return 0;//offset_[ block ];
    }

    /** \copydoc Dune::DofMapper::consecutive() const */
    bool consecutive () const
    {
      return BaseType::checkConsecutive( indexSet_ );
    }


    // Adaptation Methods (as for Index Sets)
    void resizeContainers() 
    {
      entityPolynomOrder_.update();
      for( int codim = 0; codim <= dimension; ++codim )
      {
        dofContainer_[ codim ]->update();
      }
    }

    void insertEntity ( const EntityType &entity )
    {
      resizeContainers();
      const int polOrd = polynomOrder( entity );

      const LagrangePointSetType *set = lagrangePointSet( polOrd, entity.type() );
      ForLoop< InsertSubEntities, 0, dimension> :: 
        apply( entity, *set, polOrd, size_, dofContainer_ );
    }

    void removeEntity ( const EntityType &entity )
    {}

    void resize ()
    {
      resizeContainers();
      typedef typename GridPartType :: template Codim< 0 > :: IteratorType IteratorType;
      const IteratorType end = gridPart_.template end<0>();
      for( IteratorType it = gridPart_.template begin<0>(); 
           it != end ; ++it ) 
      {
        insertEntity( *it );
      }
    }

    bool compress ()
    {
      computeOffsets();
      return true;
    }

    void read_xdr ( const char *filename, int timestep )
    {
      computeOffsets();
    }

    void write_xdr ( const char *filename, int timestep )
    {}

  private:
    // prohibit copying and assignment
    PAdaptiveLagrangeMapper ( const ThisType & );
    ThisType &operator=( const ThisType & );

    void computeOffsets ( const double overShoot = 1.0 )
    {
      // store old offset and calculate new offsets
      size_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
      {
        oldOffSet_[ codim ] = offset_[ codim ];
        offset_[ codim ] = size_;

        const unsigned int codimSize = indexSet_.size( codim ) * maxDofs_[ codim ];
        size_ += static_cast< unsigned int >( overShoot * codimSize );
      }
    }

    const GridPartType& gridPart_;
    // reference to dof manager
    DofManagerType& dm_;

    const IndexSetType &indexSet_;
    
    LagrangePointSetMapVectorType &lagrangePointSet_;

    PolyOrderContainerType entityPolynomOrder_; 

    mutable std::vector< DofContainerType* > dofContainer_; 

    // memory overshoot 
    const double overShoot_ ;

    unsigned int maxDofs_[ dimension+1 ];
    mutable unsigned int offset_[ dimension+1 ];
    mutable unsigned int oldOffSet_[ dimension+1 ];
    mutable unsigned int size_;
    unsigned int numDofs_;
  };



  // Higher Order Lagrange Mapper
  // ----------------------------
  //
  // Note: This mapper assumes that the grid is "twist-free".

#ifdef USE_TWISTFREE_MAPPER
  template< class GridPart, int polOrder >
  class PAdaptiveLagrangeMapper
  : public DofMapperDefault< PAdaptiveLagrangeMapperTraits< GridPart, polOrder > >
  {
    typedef PAdaptiveLagrangeMapper< GridPart, polOrder > ThisType;
    typedef DofMapperDefault< PAdaptiveLagrangeMapperTraits< GridPart, polOrder > > BaseType;

  public:
    typedef PAdaptiveLagrangeMapperTraits< GridPart, polOrder > Traits;
    
    //! type of the grid part
    typedef typename Traits::GridPartType GridPartType;

    //! type of entities (codim 0)
    typedef typename Traits::EntityType EntityType;

    //! type of DofMapIterator
    typedef typename Traits::DofMapIteratorType DofMapIteratorType;
 
    //! type of the underlying grid
    typedef typename GridPartType::GridType GridType;

    //! type of coordinates within the grid
    typedef typename GridType::ctype FieldType;

    //! dimension of the grid
    static const int dimension = GridType::dimension;

    //! order of the Lagrange polynoms
    static const int polynomialOrder = Traits::polynomialOrder;

    //! type of the index set
    typedef typename GridPartType::IndexSetType IndexSetType;

    //! type of the Lagrange point set
    typedef LagrangePointSet< GridPartType, polynomialOrder >
      LagrangePointSetType;
    //! type of the map for the Lagrange point sets
    typedef std::map< const GeometryType, const LagrangePointSetType* >
      LagrangePointSetMapType;

    //! type of the DoF manager
    typedef DofManager< GridType > DofManagerType;

  protected:
    typedef typename LagrangePointSetType::DofInfo DofInfo;

  private:
    struct CodimCallInterface;
    
    template< unsigned int codim >
    struct CodimCall;

    typedef CodimMap< dimension+1, CodimCall > CodimCallMapType;

  public:
    //! constructor
    PAdaptiveLagrangeMapper ( const GridPartType &gridPart,
                     LagrangePointSetMapType &lagrangePointSet )
    : dm_( DofManagerType :: instance(gridPart.grid()) ),
      indexSet_( gridPart.indexSet() ),
      lagrangePointSet_( lagrangePointSet ),
      overShoot_( Parameter::getValidValue( "fem.lagrangemapper.overshoot", double( 1.5 ), ValidateNotLess< double >( 1.0 ) ) )
    {
      numDofs_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
        maxDofs_[ codim ] = 0;
      
      typedef typename LagrangePointSetMapType :: iterator IteratorType;
      IteratorType end = lagrangePointSet_.end();
      for( IteratorType it = lagrangePointSet_.begin(); it != end; ++it )
      {
        const LagrangePointSetType *set = (*it).second;
        if( set == 0 )
          continue;

        numDofs_ = std::max( numDofs_, set->size() );
        for( int codim = 0; codim <= dimension; ++codim )
          maxDofs_[ codim ]
            = std::max( maxDofs_[ codim ], set->maxDofs( codim ) );
      }

      computeOffsets();
      dm_.addIndexSet( *this );
    }
    
    //! destructor 
    virtual ~PAdaptiveLagrangeMapper ()
    {
      dm_.removeIndexSet( *this );
    }

    //! return overall number of degrees of freedom 
    int size () const
    {
      return size_;
    }

    /** \copydoc Dune::DofMapper::begin(const EntityType &entity) const */
    DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::beginIterator, entity, *this );
    }
    
    /** \copydoc Dune::DofMapper::end(const EntityType &entity) const */
    DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::endIterator, entity, *this );
    }

    /** \copydoc Dune::DofMapper::mapToGlobal */
    int mapToGlobal ( const EntityType &entity, const int localDof ) const
    {
      // unsigned int codim, subEntity;
      const LagrangePointSetType *set = lagrangePointSet_[ entity.type() ];
      const DofInfo& dofInfo = set->dofInfo( localDof );
      
      const unsigned int codim = dofInfo.codim;
      const int subIndex = indexSet_.subIndex( entity, dofInfo.subEntity, codim );

      return offset_[ codim ] + subIndex * maxDofs_[ codim ] + dofInfo.dofNumber;
    }

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const 
    {
      const unsigned int codim = Entity::codimension;
      assert( localDof < numEntityDofs( entity ) );
      return offset_[ codim ] + indexSet_.index( entity ) * maxDofs_[ codim ] + localDof;
    }
    
    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    int maxNumDofs () const
    {
      return numDofs_;
    }

    /** \copydoc Dune::DofMapper::numDofs(const EntityType &entity) const */
    int numDofs ( const EntityType &entity ) const
    {
      return lagrangePointSet_[ entity.type() ]->size();
    }

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      // This implementation only works for nonhybrid grids (simplices or cubes)
      return maxDofs_[ Entity::codimension ];
    }
    
    /** \brief Check, whether any DoFs are associated with a codimension */
    bool contains ( int codim ) const
    {
      return true;
    }

    /** \brief Check, whether the data in a codimension has fixed size */
    bool fixedDataSize ( int codim ) const
    {
      return false;
    }

    /** \copydoc Dune::DofMapper::oldIndex */
    int oldIndex ( const int hole, const int codim ) const
    {
      const int n = maxDofs_[ codim ];
      return offset_[ codim ] + n*indexSet_.oldIndex( hole / n, codim ) + (hole % n);
    }

    /** \copydoc Dune::DofMapper::newIndex */
    int newIndex ( const int hole, const int codim ) const
    {
      const int n = maxDofs_[ codim ];
      return offset_[ codim ] + n*indexSet_.newIndex( hole / n, codim ) + (hole % n);
    }

    /** \copydoc Dune::DofMapper::numberOfHoles */
    int numberOfHoles ( const int codim ) const
    {
      return maxDofs_[ codim ] * indexSet_.numberOfHoles( codim );
    }

    /** \copydoc Dune::DofMapper::numBlocks
     */
    int numBlocks () const 
    {
      return dimension + 1;
    }

    /** \copydoc Dune::DofMapper::oldOffSet(const int block) const */
    int oldOffSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return oldOffSet_[ block ];
    }

    /** \copydoc Dune::DofMapper::offSet(const int block) const */
    int offSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return offset_[ block ];
    }

    /** \copydoc Dune::DofMapper::consecutive() const */
    bool consecutive () const
    {
      return BaseType::checkConsecutive( indexSet_ );
    }


    // Adaptation Methods (as for Index Sets)

    void insertEntity ( const EntityType &entity )
    {
      // check, whether we need to enlarge any block and compute oversized block size
      unsigned int oldUpperBound = size_;
      for( int codim = dimension ; codim >= 0; --codim )
      {
        const unsigned int codimSize = indexSet_.size( codim ) * maxDofs_[ codim ];
        const unsigned int newUpperBound = offset_[ codim ] + codimSize;
        if( newUpperBound > oldUpperBound )
          return computeOffsets( overShoot_ );
        oldUpperBound = offset_[ codim ];
      }
    }

    void removeEntity ( const EntityType &entity )
    {}

    void resize ()
    {
      computeOffsets();
    }

    bool compress ()
    {
      computeOffsets();
      return true;
    }

    void read_xdr ( const char *filename, int timestep )
    {
      computeOffsets();
    }

    void write_xdr ( const char *filename, int timestep )
    {}

  private:
    // prohibit copying and assignment
    PAdaptiveLagrangeMapper ( const ThisType & );
    ThisType &operator=( const ThisType & );

    void computeOffsets ( const double overShoot = 1.0 )
    {
      // store old offset and calculate new offsets
      size_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
      {
        oldOffSet_[ codim ] = offset_[ codim ];
        offset_[ codim ] = size_;

        const unsigned int codimSize = indexSet_.size( codim ) * maxDofs_[ codim ];
        size_ += static_cast< unsigned int >( overShoot * codimSize );
      }
    }

    // reference to dof manager needed for debug issues 
    DofManagerType& dm_;

    const IndexSetType &indexSet_;
    
    LagrangePointSetMapType &lagrangePointSet_;

    // memory overshoot 
    const double overShoot_ ;

    unsigned int maxDofs_[ dimension+1 ];
    mutable unsigned int offset_[ dimension+1 ];
    mutable unsigned int oldOffSet_[ dimension+1 ];
    mutable unsigned int size_;
    unsigned int numDofs_;
  };
#endif // #ifdef USE_TWISTFREE_MAPPER

} // end namespace Dune 

#endif // #ifndef DUNE_LAGRANGESPACE_MAPPER_HH

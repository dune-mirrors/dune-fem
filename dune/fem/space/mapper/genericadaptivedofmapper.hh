#ifndef DUNE_GENERICPADAPTIVEDOFMAPPER_HH
#define DUNE_GENERICPADAPTIVEDOFMAPPER_HH

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

  // Generic adaptive dof mapper 
  // ----------------------------

  template< class TraitsImp >
  class GenericAdaptiveDofMapper 
  : public DofMapperDefault< TraitsImp > 
  {
  public:
    typedef TraitsImp Traits;

  protected:  
    typedef GenericAdaptiveDofMapper< Traits > ThisType;
    typedef DofMapperDefault< Traits > BaseType;

    using BaseType :: asImp;

  public:
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
      CompiledLocalKeyType;
    //! type of the map for the Lagrange point sets
    typedef std::map< const GeometryType, const CompiledLocalKeyType* >
      CompiledLocalKeyMapType;

    typedef std::vector< CompiledLocalKeyMapType > 
      CompiledLocalKeyVectorType;

    //! type of the DoF manager
    typedef DofManager< GridType > DofManagerType;

    enum { minOrder = 1 };
    enum { maxOrder = polynomialOrder };
    enum { numOrders = maxOrder - minOrder + 1 };

    struct EntityDofStorage
    {
      typedef std::vector< int > DofVectorType;
      std::vector< DofVectorType > dofs_;

      GeometryType type_;
      char used_[ numOrders ];

      EntityDofStorage() : 
        dofs_( numOrders, DofVectorType() ),
        type_()
      {
        // set used to zero 
        for( int i=0; i<numOrders; ++i ) 
          used_[ i ] = 0;
      }

      void assign( const EntityDofStorage& other ) 
      {
        assert( dofs_.size() ==  other.dofs_.size() );
        type_ = other.type_;

        // set used to zero 
        for( int k=0; k<numOrders; ++k ) 
        {
          used_[ k ] = other.used_[ k ];
          DofVectorType& dofs = dofs_[ k ];
          const DofVectorType& otherDofs = other.dofs_[ k ];
          const int dofSize = otherDofs.size();
          dofs.resize( dofSize );
          for( int d = 0; d<dofSize; ++d ) 
            dofs[ d ] = otherDofs[ d ];
        }
      }

      EntityDofStorage( const EntityDofStorage& other ) 
        : dofs_( numOrders )
      {
        assign( other );
      }

      EntityDofStorage& operator= ( const EntityDofStorage& other )
      {
        assign( other );
        return *this;
      }

      bool exists( const int codim, const int polOrd ) const 
      {
        const int entry = determineVectorEntry( codim, polOrd );
        return dofs_[ entry ].size() > 0 ;
      }

      //! returns true if entry has a reference count of 1 
      bool use( const int codim, const int polOrd ) 
      {
        const int entry = determineVectorEntry( codim, polOrd ) ;
        ++used_[ entry ];
        return ( used_[ entry ] == 1 );
      }

      void insert( const GeometryType type, 
                   const int codim,
                   const int polOrd, 
                   const int numDofs, const int startDof ) 
      {
        use( codim, polOrd );
        assert( ! exists ( codim, polOrd ) );
        {
          type_ = type ;
          DofVectorType& dofs = dofs_[ determineVectorEntry( codim, polOrd ) ];

          dofs.resize( numDofs );
          for(int i=0, dof=startDof ; i<numDofs; ++i, ++dof ) 
            dofs[ i ] = dof;
        }
      }

      int determineVectorEntry( const int codim, const int polOrd ) const 
      {
        assert( codim >= 0 );
        assert( codim <= dimension );
        // also for codim == 0 we more then one storage because of different number of
        // dofs per polynmomial degree 
        return (codim < dimension) ? (polOrd-minOrder) : 0;
      }

      const GeometryType& type () const { return type_ ; }

      void remove( const int codim, const int polOrd ) 
      {
        const int entry = determineVectorEntry( codim, polOrd );
        if( used_[ entry ] > 0 ) 
          --used_[ entry ] ;
      }

      void reset() 
      {
        for( int k=0; k<numOrders; ++k ) 
          used_[ k ] = 0;
      }

      int dof ( const int codim, const int polOrd, const size_t dofNumber ) const 
      { 
        const int entry = determineVectorEntry( codim, polOrd );
        assert( dofNumber < dofs_[ entry ].size() );
        return dofs_[ entry ][ dofNumber ];
      }

      int entityDof ( int dofNumber ) const 
      { 
        for( int k = 0; k<numOrders; ++k ) 
        {
          const int dofSize = dofs_[ k ].size();
          if( dofNumber < dofSize ) 
            return dofs_[ k ][ dofNumber ];
          else 
            dofNumber -= dofSize;
        }
        // we should not get here 
        assert( false );
        abort();
        return -1;
      }

      int entityDofs () const 
      {
        int dofSize = 0;
        for( int k = 0; k<numOrders; ++k ) 
        {
          dofSize += dofs_[ k ].size();
        }
        return dofSize;
      }

      template <class VectorType> 
      void detectUnusedDofs( VectorType& isHole, 
                             const int actSize ) 
      {
        for( int k=0; k<numOrders; ++k )
        {
          DofVectorType& dofs = dofs_[ k ];
          const int dofSize = dofs.size();

          if( dofSize > 0 ) 
          {
            if( used_[  k ] )
            {
              for( int d = 0; d<dofSize; ++d ) 
              {
                const int dof = dofs[ d ] ;
                if( dof < actSize )
                {
                  assert( dof < (int)isHole.size() );
                  isHole[ dof ] = false ;
                }
              }
            }
            else 
            {
              dofs.resize( 0 );
            }
          }
        }
      }

      void printDofs() const 
      {
        for( int k = 0; k<numOrders; ++k )
        {
          const DofVectorType& dofs = dofs_[ k ];
          const int dofSize = dofs.size();
          for( int d = 0; d<dofSize; ++d ) 
            std::cout << dofs[ d ] << " dofs " << std::endl;
        }
      }

      template <class VectorType> 
      bool removeHoles( VectorType& oldIdx, VectorType& newIdx, 
                        VectorType& holesVec, int& currentHole,  
                        const int usedSize, int& holes ) 
      {
        bool haveToCopy = false ;
        for( int k=0; k<numOrders; ++k )
        {
          DofVectorType& dofs = dofs_[ k ];
          const int dofSize = dofs.size();
          for( int dof = 0; dof<dofSize; ++dof ) 
          {
            assert( used_[ k ] );
            // get global DoF number 
            int& currDof = dofs[ dof ] ;

            // if dof >= usedSize it has to be migrated to a hole 
            if( currDof >= usedSize ) 
            {
              // serach next hole that is smaler than actual size 
              --currentHole;

              // if currentHole < 0 then error, because we have index larger then
              // actual size 
              assert(currentHole >= 0);
              assert( holesVec[currentHole] < usedSize );

              // remember old and new index 
              oldIdx[ holes ] = currDof;
              currDof = holesVec[ currentHole ];
              newIdx[ holes ] = currDof ;

              // increase number of holes 
              ++holes;

              haveToCopy = true;
            }
          }
        }
        return haveToCopy;
      }
    };

    typedef EntityDofStorage EntityDofStorageType;

    struct PolynomOrderStorage 
    {
      unsigned char k_;
      unsigned char active_ ;
      PolynomOrderStorage() : k_( maxOrder ), active_( 0 ) {}
      PolynomOrderStorage( const int k ) : k_( k ), active_( 0 ) {}
      int order () const { return k_;}
      void set ( const int k ) { k_ = k; active_ = 1 ; }
      void activate() { active_ = 1; }
      bool active () const { return active_; }
      bool deactivate ( int& k ) 
      { 
        k = k_;
        if( active_ ) 
        {
          active_ = 0 ;
          return true ;
        }
        return false;
      }
    };

    typedef PolynomOrderStorage  PolynomOrderStorageType;

    typedef PersistentContainer< GridType, EntityDofStorageType > DofContainerType ;

    typedef PersistentContainer< GridType, PolynomOrderStorageType > PolyOrderContainerType ;
  protected:
    typedef typename CompiledLocalKeyType::DofInfo DofInfo;

    template < int codim > 
    struct InsertSubEntities
    {
      static void insertDofs( const EntityType& entity, 
                              const CompiledLocalKeyType& clk,
                              const int polOrd,
                              const int subEntity,
                              unsigned int& globalSize,
                              unsigned int& notAlreadyCounted,
                              EntityDofStorage& entityDofs )
      {
        const int numDofs = clk.numDofs( codim, subEntity );
        // only if dofs exists on this entity do something
        if( numDofs > 0 ) 
        {
          bool notCountedYet = false ;
          if( ! entityDofs.exists( codim, polOrd ) ) 
          {
            entityDofs.insert( entity.type(), codim, polOrd, numDofs, globalSize );
            globalSize += numDofs;
            notCountedYet = true ;
          }
          else 
          {
            // if refcount is only 1 then true is returned 
            notCountedYet = entityDofs.use( codim, polOrd );
          }

          // if not counted yet, count !
          if( notCountedYet ) 
          {
            notAlreadyCounted += numDofs;
          }
        }
      }

      static void apply( const EntityType& entity, 
                         const CompiledLocalKeyType& clk,
                         const int polOrd,
                         unsigned int& globalSize,
                         unsigned int& notAlreadyCounted,
                         std::vector< DofContainerType* > dofContainers ) 
      {
        DofContainerType& dofContainer = *dofContainers[ codim ];
        if( codim == 0 ) 
        {
          insertDofs( entity, clk, polOrd, 0, globalSize, 
                      notAlreadyCounted, dofContainer[ entity ] );
        }
        else 
        {
          const int count = entity.template count< codim > ();
          for(int i=0; i<count; ++i ) 
          {
            insertDofs( entity, clk, polOrd, i, globalSize, 
                        notAlreadyCounted, dofContainer( entity, i ) );
          }
        }
      }
    };

    template < int codim > 
    struct RemoveSubEntities
    {
      static void apply( const EntityType& entity, 
                         const int polOrd,
                         std::vector< DofContainerType* > dofContainers ) 
      {
        DofContainerType& dofContainer = *dofContainers[ codim ];
        const int count = entity.template count< codim > ();
        for(int i=0; i<count; ++i ) 
        {
          EntityDofStorage& entityDofs = dofContainer( entity, i );
          entityDofs.remove( codim, polOrd );
        }
      }
    };

  public:
    //! constructor
    GenericAdaptiveDofMapper ( const GridPartType &gridPart,
                              CompiledLocalKeyVectorType &compiledLocalKeyVector )
    : gridPart_( gridPart ),
      dm_( DofManagerType :: instance(gridPart.grid()) ),
      compiledLocalKeys_( compiledLocalKeyVector ),
      entityPolynomOrder_( gridPart.grid(), 0 ),
      dofContainer_( dimension+1, (DofContainerType *) 0 ),
      numberOfHoles_( 0 ),
      oldIndex_(),
      newIndex_(),
      size_(0),
      maxNumDofs_( 0 ),
      sequence_( dm_.sequence() )
    {
      /*
      PolynomOrderStorage p;
      std::cout << sizeof( p ) << " size of polStorage" << std::endl;
      EntityDofStorage en;
      std::cout << sizeof( en ) << " size of enStorage" << std::endl;
      GeometryType type;
      std::cout << sizeof( type) << " size of GeomType " << std::endl;
      */

      for( int codim = 0; codim <= dimension; ++codim )
        dofContainer_[ codim ] = new DofContainerType( gridPart.grid(), codim );

      for( size_t i=0; i<compiledLocalKeys_.size(); ++i ) 
      {
        typedef typename CompiledLocalKeyMapType :: iterator IteratorType;
        IteratorType end = compiledLocalKeys_[ i ].end();
        for( IteratorType it = compiledLocalKeys_[ i ].begin(); it != end; ++it )
        {
          const CompiledLocalKeyType *clk = (*it).second;
          if( clk == 0 ) continue;
          
          maxNumDofs_ = std :: max( maxNumDofs_, clk->size() );
        }
      }

      resize();
      // register to DofManager for adaptation 
      dm_.addIndexSet( asImp() );
    }

    //! sort of copy constructor
    GenericAdaptiveDofMapper ( const GenericAdaptiveDofMapper& other,
                               CompiledLocalKeyVectorType &compiledLocalKeyVector )
    : gridPart_( other.gridPart_ ),
      dm_( other.dm_ ),
      compiledLocalKeys_( compiledLocalKeyVector ),
      entityPolynomOrder_( other.entityPolynomOrder_ ),
      dofContainer_( dimension+1, (DofContainerType *) 0 ),
      numberOfHoles_( other.numberOfHoles_ ),
      oldIndex_( other.oldIndex_ ),
      newIndex_( other.newIndex_ ),
      size_( other.size_ ),
      maxNumDofs_( other.maxNumDofs_ ),
      sequence_( other.sequence_ )
    {
      for( int codim = 0; codim <= dimension; ++codim )
        dofContainer_[ codim ] = new DofContainerType( *(other.dofContainer_[ codim ]) );

      dm_.addIndexSet( asImp() );
    }

    int polynomOrder( const EntityType& entity ) const 
    {
      return entityPolynomOrder_[ entity ].order();
    }

    void setPolynomOrder( const EntityType& entity, const int polOrd ) 
    {
      if( polOrd < 1 || polOrd > polynomialOrder ) 
        return ;

      entityPolynomOrder_[ entity ].set( polOrd );
    }

    DofContainerType& dofContainer( const size_t codim ) const 
    {
      assert( codim < dofContainer_.size() );
      return *(dofContainer_[ codim ]);
    }

    const CompiledLocalKeyType* 
    compiledLocalKey( const int polOrd, const GeometryType type ) const 
    {
      // add polOrd here 
      return compiledLocalKeys_[ polOrd ][ type ];
    }
    
    //! destructor 
    virtual ~GenericAdaptiveDofMapper ()
    {
      dm_.removeIndexSet( asImp() );
    }

    //! return overall number of degrees of freedom 
    int size () const
    {
      return size_;
    }

    /** \copydoc Dune::DofMapper::begin(const EntityType &entity) const */
    DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::beginIterator, entity, asImp() );
    }
    
    /** \copydoc Dune::DofMapper::end(const EntityType &entity) const */
    DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::endIterator, entity, asImp() );
    }

    /** \copydoc Dune::DofMapper::mapToGlobal */
    int mapToGlobal ( const EntityType &entity, const int localDof ) const
    {
      const int polOrd = polynomOrder( entity );
      const DofInfo &dofInfo = compiledLocalKey( polOrd, entity.type() )->dofInfo( localDof );

      const unsigned int codim = dofInfo.codim;
      const unsigned int subEntity = dofInfo.subEntity;

      return dofContainer( codim )( entity, subEntity ).dof( codim, polOrd, dofInfo.dofNumber );
    }

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const 
    {
      return dofContainer( Entity :: codimension )[ entity ].entityDof( localDof );
    }
    
    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    int maxNumDofs () const
    {
      return maxNumDofs_;
    }

    /** \copydoc Dune::DofMapper::numDofs(const EntityType &entity) const */
    int numDofs ( const EntityType &entity ) const
    {
      const int polOrd = polynomOrder( entity );
      return compiledLocalKey( polOrd, entity.type() )->size();
    }

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      // This implementation only works for nonhybrid grids (simplices or cubes)
      return dofContainer( Entity :: codimension )[ entity ].entityDofs();
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
    int oldIndex ( const int hole, const int block ) const
    {
      assert( oldIndex_[ hole ] >= 0 );
      return oldIndex_[ hole ];
    }

    /** \copydoc Dune::DofMapper::newIndex */
    int newIndex ( const int hole, const int block ) const
    {
      assert( newIndex_[ hole ] >= 0 );
      return newIndex_[ hole ];
    }

    /** \copydoc Dune::DofMapper::numberOfHoles */
    int numberOfHoles ( const int block ) const
    {
      return numberOfHoles_;
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
      return 0;
    }

    /** \copydoc Dune::DofMapper::newOffset
     */
    int offSet ( const int block ) const
    {
      return 0;
    }

    /** \copydoc Dune::DofMapper::consecutive() const */
    bool consecutive () const
    {
      return true;
    }

    // Adaptation Methods (as for Index Sets)
    void resizeContainers() 
    {
      entityPolynomOrder_.update();
      for( int codim = 0; codim <= dimension; ++codim )
      {
        dofContainer( codim ).update();
      }
    }

    void insertEntity ( const EntityType &entity )
    {
      resizeContainers();
      insertEntityDofs( entity );
    }

    // return number of local dofs that were not visited yet 
    unsigned int insertEntityDofs( const EntityType &entity )
    {
      PolynomOrderStorageType& polyStorage = entityPolynomOrder_[ entity ];
      if( ! polyStorage.active() )
      {
        unsigned int notAlreadyCounted = 0;

        const int polOrd = polyStorage.order();
        // get lagrange point set 
        const CompiledLocalKeyType *clk = compiledLocalKey( polOrd, entity.type() );

        //std::cout << "Insert Entity " << gridPart_.grid().localIdSet().id( entity ) << std::endl;
        
        polyStorage.activate();
        ForLoop< InsertSubEntities, 0, dimension> :: 
          apply( entity, *clk, polOrd, size_, notAlreadyCounted, dofContainer_ );

        //printEntityDofs( entity );
        return notAlreadyCounted ;
      }

      return 0;
    }

    void removeEntity ( const EntityType &entity )
    {
      int polOrd; 
      if( entityPolynomOrder_[ entity ].deactivate( polOrd ) )
      {
        ForLoop< RemoveSubEntities, 0, dimension> :: 
          apply( entity, polOrd, dofContainer_ );
      }
    }

    void resize ()
    {
      resizeContainers();
      insertAllUsed();
    }

    //! adjust mapper to newly set polynomial orders 
    void adapt()
    {
      resizeContainers();
      sequence_ = -1;
      compress();
    }

    //! return number of DoFs currently used for space
    size_t insertAllUsed() 
    {
      // reset all dof entries 
      setUnused();
      
      // count current size 
      size_t usedSize = 0;

      typedef typename GridPartType :: template Codim< 0 > :: IteratorType IteratorType;
      const IteratorType end = gridPart_.template end<0>();
      for( IteratorType it = gridPart_.template begin<0>(); 
           it != end ; ++it ) 
      {
        // number of new dofs (not already counted) is returned 
        usedSize += insertEntityDofs( *it );
      }

      //std::cout << "Insert Size = " << size_ << std::endl;
      //printDofs();
      //std::cout << "Insert done! " << std::endl;
      return usedSize ;
    }

    void printDofs() const 
    {
      //std::cout << "Size " << size_ << std::endl;
      typedef typename GridPartType :: template Codim< 0 > :: IteratorType IteratorType;
      const IteratorType end = gridPart_.template end<0>();
      for( IteratorType it = gridPart_.template begin<0>(); 
           it != end ; ++it ) 
      {
        printEntityDofs( *it );
      }
    }

    void printEntityDofs( const EntityType& entity ) const 
    {
      std::cout << "Print entity " << gridPart_.grid().localIdSet().id( entity ) << " with " << std::endl;
      for( int i = 0; i<numDofs( entity ); ++i ) 
      {
        std::cout << "en[ " << i << " ] = " << mapToGlobal( entity, i ) << std::endl;
      }
    }

    //! reset all used flags of all DoF entries 
    void setUnused() 
    {
      {
        typedef typename PolyOrderContainerType :: Iterator Iterator;
        const Iterator endit = entityPolynomOrder_.end();
        for( Iterator it = entityPolynomOrder_.begin(); it != endit; ++it )
        {
          PolynomOrderStorageType& p = *it;
          int pOrd; 
          p.deactivate( pOrd );
        }
      }

      for( int codim = 0; codim <= dimension; ++codim ) 
      {
        DofContainerType& codimContainer = dofContainer( codim );
        typedef typename DofContainerType :: Iterator Iterator;
        const Iterator endit = codimContainer.end();
        for( Iterator it = codimContainer.begin(); it != endit; ++it )
        {
          it->reset();
        }
      }
    }

    // --compress 
    bool compress ()
    {
      if( sequence_ == dm_.sequence() ) 
      {
        numberOfHoles_ = 0;
        return false ;
      }

      // make all dofs that are currently used 
      // (corresponding to GridPart's iterator)
      const size_t usedSize = insertAllUsed();

      // reset number of holes 
      numberOfHoles_ = 0;

      // true if a least one dof must be copied 
      bool haveToCopy = false;

      std::vector<int> validHoles;

      {
        // default is true which means entry is hole
        std::vector< bool > holeMarker( size_, true );

        // mark holes 
        for( int codim = 0; codim <= dimension; ++codim ) 
        {
          DofContainerType& codimContainer = dofContainer( codim );
          typedef typename DofContainerType :: Iterator Iterator;
          const Iterator endit = codimContainer.end();
          for( Iterator it = codimContainer.begin(); it != endit; ++it )
          {
            EntityDofStorageType& dof = *it;
            // store holes if unused dofs exits, otherwise increase actSize 
            dof.detectUnusedDofs( holeMarker, usedSize );
          }
        }

        // check for real holes 
        validHoles.reserve( usedSize );
        validHoles.resize( 0 );

        // check the vector of hole flags and store numbers
        for( int i=0; i<usedSize; ++i ) 
        {
          // it's only a valid hole, if it was not marked otherwise 
          if( holeMarker[ i ] )  
          {
            //std::cout << "Dof number " << i << " is not used" << std::endl;
            validHoles.push_back( i );
          }
        }
      }

      if( validHoles.size() > 0 ) 
      {
        // counter for current hole to use 
        int currentHole = validHoles.size();

        // resize old-new index storage to correct size 
        oldIndex_.resize( currentHole, -1) ;
        newIndex_.resize( currentHole, -1) ;

        for( int codim = 0; codim <= dimension; ++codim ) 
        {
          DofContainerType& codimContainer = dofContainer( codim );
          typedef typename DofContainerType :: Iterator Iterator;
          const Iterator endit = codimContainer.end();
          for( Iterator it = codimContainer.begin(); it != endit; ++it )
          {
            // get dof storage 
            EntityDofStorageType& dof = *it;

            // replace DoF larger than size with DoF from holes 
            // set haveToCopy to true if true is returned 
            haveToCopy |= dof.removeHoles( oldIndex_, newIndex_, 
                                           validHoles, currentHole, 
                                           usedSize, numberOfHoles_ );
          }
        }
      }

      // store new size 
      size_ = usedSize ;

      //std::cout << "Size " << size_ << " holes " << numberOfHoles_ << std::endl;
      //printDofs();

      sequence_ = dm_.sequence();

      return haveToCopy;
    }

    void read_xdr ( const char *filename, int timestep )
    {
    }

    void write_xdr ( const char *filename, int timestep )
    {}

  private:
    // prohibit copying and assignment
    GenericAdaptiveDofMapper ( const ThisType & );
    ThisType &operator=( const ThisType & );

    const GridPartType& gridPart_;
    // reference to dof manager
    DofManagerType& dm_;

    CompiledLocalKeyVectorType &compiledLocalKeys_;

    PolyOrderContainerType entityPolynomOrder_; 

    mutable std::vector< DofContainerType* > dofContainer_; 

    int numberOfHoles_ ;
    std::vector< int > oldIndex_ ;
    std::vector< int > newIndex_ ;

    mutable unsigned int size_;
    unsigned int maxNumDofs_;
    int sequence_ ;
  };

} // end namespace Dune 

#endif // #ifndef DUNE_LAGRANGESPACE_MAPPER_HH
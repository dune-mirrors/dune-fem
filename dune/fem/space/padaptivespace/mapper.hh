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

    int polynomOrder( const EntityType& entity ) const 
    {
      return 1;
    }

    void setPolynomOrder( const EntityType& entity, const int polOrd ) 
    {
    }
  };



  // Second Order Lagrange Mapper
  // ----------------------------

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
        dofs_( numOrders ),
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

      void use( const int codim, const int polOrd ) 
      {
        ++used_[ determineVectorEntry( codim, polOrd ) ];
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
        return (codim > 0 && codim < dimension) ? (polOrd-minOrder) : 0;
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
      void detectUnusedDofs( VectorType& holes, int& actHoles, int& actSize ) 
      {
        for( int k=0; k<numOrders; ++k )
        {
          DofVectorType& dofs = dofs_[ k ];
          const int dofSize = dofs.size();

          if( used_[  k ] ) 
          {
            //for( int d = 0; d<dofSize; ++d ) 
            //  std::cout << dofs[ d ] << " used dofs " << std::endl;
            actSize += dofSize ;
          }
          else 
          {

          //if( used_ [ k ] <= 0 && dofSize > 0 )
          //{
            for( int d = 0; d<dofSize; ++d ) 
            {
              holes[ actHoles++ ] = dofs[ d ];
              //std::cout << "add dof " << dofs[ d ] << " to list" << std::endl;
            }
            dofs.resize( 0 );
          }
          //else 
          //  actSize += dofSize ;
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
                        VectorType& holesVec, int& actHole, int& holes,
                        const int actSize ) 
      {
        bool haveToCopy = false ;
        for( int k=0; k<numOrders; ++k )
        {
          DofVectorType& dofs = dofs_[ k ];
          const int dofSize = dofs.size();
          for( int dof = 0; dof<dofSize; ++dof ) 
          {
            assert( used_[ k ] );
            int& currDof = dofs[ dof ] ;
            if( currDof >= actSize ) 
            {
              // serach next hole that is smaler than actual size 
              --actHole;
              // if actHole < 0 then error, because we have index larger then
              // actual size 
              assert(actHole >= 0);
              while ( holesVec[actHole] >= actSize )
              {
                --actHole;
                if(actHole < 0) break;
              }

              assert(actHole >= 0);

              // remember old and new index 
              oldIdx[holes] = currDof;
              currDof = holesVec[actHole];
              newIdx[holes] = currDof ;
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
    typedef typename LagrangePointSetType::DofInfo DofInfo;

    template < int codim > 
    struct InsertSubEntities
    {
      static void insertDofs( const EntityType& entity, 
                              const LagrangePointSetType& set,
                              const int polOrd,
                              const int subEntity,
                              unsigned int&  dofCounter,
                              EntityDofStorage& entityDofs )
      {
        if( ! entityDofs.exists( codim, polOrd ) ) 
        {
          const int numDofs = set.numDofs( codim, subEntity );
          entityDofs.insert( entity.type(), codim, polOrd, numDofs, dofCounter );
          dofCounter += numDofs;
        }
        else 
          entityDofs.use( codim, polOrd );
      }
      static void apply( const EntityType& entity, 
                         const LagrangePointSetType& set,
                         const int polOrd,
                         unsigned int&  dofCounter,
                         std::vector< DofContainerType* > dofContainers ) 
      {
        DofContainerType& dofContainer = *dofContainers[ codim ];
        if( codim == 0 ) 
        {
          insertDofs( entity, set, polOrd, 0, dofCounter, 
                      dofContainer[ entity ] );
        }
        else 
        {
          const int count = entity.template count< codim > ();
          for(int i=0; i<count; ++i ) 
          {
            insertDofs( entity, set, polOrd, i, dofCounter, 
                        dofContainer( entity, i ) );
          }
        }
      }
    };

    template < int codim > 
    struct RemoveSubEntities
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
          entityDofs.remove( codim, polOrd );
        }
      }
    };

  public:
    //! constructor
    PAdaptiveLagrangeMapper ( const GridPartType &gridPart,
                     LagrangePointSetMapVectorType &lagrangePointSetVector )
    : gridPart_( gridPart ),
      dm_( DofManagerType :: instance(gridPart.grid()) ),
      lagrangePointSet_( lagrangePointSetVector ),
      entityPolynomOrder_( gridPart.grid(), 0 ),
      dofContainer_( dimension+1, (DofContainerType *) 0 ),
      numberOfHoles_( 0 ),
      oldIndex_(),
      newIndex_(),
      size_(0),
      sequence_( dm_.sequence() )
    {
      PolynomOrderStorage p;
      std::cout << sizeof( p ) << " size of polStorage" << std::endl;
      EntityDofStorage en;
      std::cout << sizeof( en ) << " size of enStorage" << std::endl;
      GeometryType type;
      std::cout << sizeof( type) << " size of GeomType " << std::endl;

      maxNumDofs_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
        dofContainer_[ codim ] = new DofContainerType( gridPart.grid(), codim );

      for( size_t i=0; i<lagrangePointSet_.size(); ++i ) 
      {
        typedef typename LagrangePointSetMapType :: iterator IteratorType;
        IteratorType end = lagrangePointSet_[ i ].end();
        for( IteratorType it = lagrangePointSet_[ i ].begin(); it != end; ++it )
        {
          const LagrangePointSetType *set = (*it).second;
          if( set == 0 ) continue;
          
          maxNumDofs_ = std :: max( maxNumDofs_, set->size() );
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
      const DofInfo &dofInfo = lagrangePointSet( polOrd, entity.type() )->dofInfo( localDof );

      const unsigned int codim = dofInfo.codim;
      const unsigned int subEntity = dofInfo.subEntity;

      return dofContainer( codim )( entity, subEntity ).dof( codim, polOrd, dofInfo.dofNumber );
    }

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const 
    {
      return dofContainer( Entity :: codimension )( entity ).entityDof( localDof );
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
      return lagrangePointSet( polOrd, entity.type() )->size();
    }

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      // This implementation only works for nonhybrid grids (simplices or cubes)
      return dofContainer( Entity :: codimension )( entity ).entityDofs();
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
      return oldIndex_[ hole ];
    }

    /** \copydoc Dune::DofMapper::newIndex */
    int newIndex ( const int hole, const int block ) const
    {
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

    void insertEntityDofs( const EntityType &entity )
    {
      PolynomOrderStorageType& polyStorage = entityPolynomOrder_[ entity ];
      if( ! polyStorage.active() )
      {
        const int polOrd = polyStorage.order();

        //std::cout << "Insert Entity " << gridPart_.grid().localIdSet().id( entity ) << std::endl;
        
        polyStorage.activate();
        const LagrangePointSetType *set = lagrangePointSet( polOrd, entity.type() );
        ForLoop< InsertSubEntities, 0, dimension> :: 
          apply( entity, *set, polOrd, size_, dofContainer_ );

        //printEntityDofs( entity );
      }
    }

    void removeEntity ( const EntityType &entity )
    {
      int polOrd; 
      if( entityPolynomOrder_[ entity ].deactivate( polOrd ) )
      {
        //std::cout << "Remove Entity " << gridPart_.grid().localIdSet().id( entity ) << " with " << polOrd
        //  << std::endl;

        const LagrangePointSetType *set = lagrangePointSet( polOrd, entity.type() );
        ForLoop< RemoveSubEntities, 0, dimension> :: 
          apply( entity, *set, polOrd, size_, dofContainer_ );
      }
    }

    void resize ()
    {
      resizeContainers();
      insertAllUsed();
    }

    void insertAllUsed() 
    {
      setUnused();

      typedef typename GridPartType :: template Codim< 0 > :: IteratorType IteratorType;
      const IteratorType end = gridPart_.template end<0>();
      for( IteratorType it = gridPart_.template begin<0>(); 
           it != end ; ++it ) 
      {
        insertEntityDofs( *it );
      }

      //std::cout << "Size " << size_ << std::endl;
      //printDofs();
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
      insertAllUsed();

      numberOfHoles_ = 0;
      const int sizeOfVecs = size_;

      // true if a least one dof must be copied 
      bool haveToCopy = false;

      int countValidHoles = 0;
      int actHole = 0;
      int newActSize = 0;
      std::vector<int> validHoles; 
      {
        std::vector< int > allHoles( sizeOfVecs, -1 ) ;
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
            dof.detectUnusedDofs( allHoles, actHole, newActSize );
          }
        }

        // check for real holes 
        validHoles.resize( actHole, -1 );

        // remove invalid holes from list 
        for( int i=0; i<actHole; ++i ) 
        {
          // it's only a valid hole, if it is smaller then the current size 
          if( allHoles[ i ] < newActSize ) 
          {
            validHoles[ countValidHoles ] = allHoles[ i ];
            ++countValidHoles;
          }
        }
      }

      if( countValidHoles > 0 ) 
      {
        oldIndex_.resize( countValidHoles, -1) ;
        newIndex_.resize( countValidHoles, -1) ;
        for( int codim = 0; codim <= dimension; ++codim ) 
        {
          DofContainerType& codimContainer = dofContainer( codim );
          typedef typename DofContainerType :: Iterator Iterator;
          const Iterator endit = codimContainer.end();
          for( Iterator it = codimContainer.begin(); it != endit; ++it )
          {
            EntityDofStorageType& dof = *it;
            bool haveTo = dof.removeHoles( oldIndex_, newIndex_, validHoles, countValidHoles, numberOfHoles_, newActSize );
            if( haveTo ) haveToCopy = true ;
          }
        }
      }

      size_ = newActSize;
      /*
      for( int i=0; i<numberOfHoles_; ++i ) 
      {
        std::cout << "old[ " << i << " ] = " << oldIndex_[ i ];
        std::cout << " new[ " << i << " ] = " << newIndex_[ i ] <<std::endl;
      }
      */

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
    PAdaptiveLagrangeMapper ( const ThisType & );
    ThisType &operator=( const ThisType & );

    const GridPartType& gridPart_;
    // reference to dof manager
    DofManagerType& dm_;

    LagrangePointSetMapVectorType &lagrangePointSet_;

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

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
      GeometryType type_;
      std::vector< std::vector< int > > dofs_;
      std::vector< int > used_;

      EntityDofStorage() : 
        type_(),
        dofs_( polynomialOrder+1 ),
        used_(  polynomialOrder+1, int(0) )
      {}

      bool exists( const int polOrd ) const 
      {
        return dofs_[ polOrd ].size() > 0 ;
      }

      bool used ( const int polOrd ) const 
      {
        return used_[ polOrd ] > 0;
      }

      void insert( const GeometryType type, const int polOrd, const int numDofs, const int startDof ) 
      {
        assert( ! exists ( polOrd ) );
        {
          type_ = type ;
          ++used_[ polOrd ] ;
          dofs_[ polOrd ].resize( numDofs );
          for(int i=0, dof=startDof ; i<numDofs; ++i, ++dof ) 
            dofs_[ polOrd ][ i ] = dof;
        }
      }

      const GeometryType& type () const { return type_ ; }

      void remove( const int polOrd ) 
      {
        if( used_[ polOrd ]  > 0 ) 
          --used_[ polOrd ] ;
      }

      int dof ( const int polOrd, const int dofNumber ) const 
      { 
        return dofs_[ polOrd ][ dofNumber ];
      }

      int entityDof ( int dofNumber ) const 
      { 
        for( int k = 1; k<=polynomialOrder; ++k ) 
        {
          const int dofSize = dofs_[ k ].size();
          if( dofNumber < dofSize ) 
            return dofs_[ k ][ dofNumber ];
          else 
            dofNumber -= dofSize;
        }
        assert( false );
        abort();
        return -1;
      }

      int entityDofs () const 
      {
        int dofSize = 0;
        for( int k = 1; k<=polynomialOrder; ++k ) 
        {
          dofSize += dofs_[ k ].size();
        }
        return dofSize;
      }

      template <class VectorType> 
      void fillHoles( VectorType& holes, int& actHoles, int& actSize ) 
      {
        for( int k = 1; k<=polynomialOrder; ++k )
        {
          const int dofSize = dofs_[ k ].size();
          if( used_ [ k ] <= 0 && dofSize > 0 )
          {
            for( int d = 0; d<dofSize; ++d, ++actHoles ) 
            {
              holes[ actHoles ] = dofs_[ k ][ d ];
            }
            dofs_[ k ].resize( 0 );
          }
          else 
            actSize += dofSize ;
        }
      }

      template <class VectorType> 
      bool checkRanges( VectorType& oldIdx, VectorType& newIdx, 
                        VectorType& holesVec, int& actHole, int& holes,
                        const int actSize ) 
      {
        bool haveToCopy = false ;
        for( int k=1; k<=polynomialOrder; ++k )
        {
          const int dofSize = dofs_[ k ].size();
          for( int dof = 0; dof<dofSize; ++dof ) 
          {
            int& currDof = dofs_[ k ][ dof ] ;
            if( currDof >= actSize)
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
      PolynomOrderStorage() : k_( polynomialOrder ), active_( 1 ) {}
      PolynomOrderStorage( const int k ) : k_( k ), active_( 1 ) {}
      int order () const { return k_;}
      bool deactivate ( int& polOrd ) 
      { 
        polOrd = k_;
        return true ;
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
            entityDofs.insert( entity.type(), polOrd, numDofs, dofCounter );
            dofCounter += numDofs;
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
          entityDofs.remove( polOrd );
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
      overShoot_( Parameter::getValidValue( "fem.lagrangemapper.overshoot", double( 1.5 ), ValidateNotLess< double >( 1.0 ) ) ),
      size_(0),
      numberOfHoles_( 0 )
    {
      numDofs_ = 0;
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
      const DofInfo &dofInfo = lagrangePointSet( polOrd, entity.type() )->dofInfo( localDof );

      const unsigned int codim = dofInfo.codim;
      const unsigned int subEntity = dofInfo.subEntity;

      return dofContainer( codim )( entity, subEntity ).dof( polOrd, dofInfo.dofNumber );
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
      const int polOrd = polynomOrder( entity );

      //std::cout << "Insert Entity " << indexSet_.index( entity ) << std::endl;

      const LagrangePointSetType *set = lagrangePointSet( polOrd, entity.type() );
      ForLoop< InsertSubEntities, 0, dimension> :: 
        apply( entity, *set, polOrd, size_, dofContainer_ );
    }

    void removeEntity ( const EntityType &entity )
    {
      int polOrd; 
      if( entityPolynomOrder_[ entity ].deactivate( polOrd ) )
      {
        std::cout << "Remove Entity " << indexSet_.index( entity ) << " with " << polOrd
          << std::endl;

        const LagrangePointSetType *set = lagrangePointSet( polOrd, entity.type() );
        ForLoop< RemoveSubEntities, 0, dimension> :: 
          apply( entity, *set, polOrd, size_, dofContainer_ );
      }
    }

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

      std::cout << "Size " << size_ << std::endl;
    }

    // --compress 
    bool compress ()
    {
      numberOfHoles_ = 0;
      const int sizeOfVecs = size_;
      std::vector< int > holes_( sizeOfVecs, -1 ) ;

      // true if a least one dof must be copied 
      bool haveToCopy = false;

      // mark holes 
      int actHole = 0;
      int newActSize = 0;
      for( int codim = 0; codim <= dimension; ++codim ) 
      {
        DofContainerType& codimContainer = dofContainer( codim );
        typedef typename DofContainerType :: Iterator Iterator;
        const Iterator endit = codimContainer.end();
        for( Iterator it = codimContainer.begin(); it != endit; ++it )
        {
          EntityDofStorageType& dof = *it;
          // store holes if unused dofs exits, otherwise increase actSize 
          dof.fillHoles( holes_, actHole, newActSize );
        }
      }

      oldIndex_.resize( actHole, -1) ;
      newIndex_.resize( actHole, -1) ;

      std::cout << "Current real size " << newActSize << std::endl;
      //for( int i=0; i<actHole; ++i ) 
      //  std::cout << "Hole[ " << i << " ] = " << holes_[ i ] << std::endl;

      if( actHole > 0 ) 
      {
        for( int codim = 0; codim <= dimension; ++codim ) 
        {
          DofContainerType& codimContainer = dofContainer( codim );
          typedef typename DofContainerType :: Iterator Iterator;
          const Iterator endit = codimContainer.end();
          for( Iterator it = codimContainer.begin(); it != endit; ++it )
          {
            EntityDofStorageType& dof = *it;
            bool haveTo = dof.checkRanges( oldIndex_, newIndex_, holes_, actHole, numberOfHoles_, newActSize );
            if( haveTo ) haveToCopy = true ;
          }
        }
      }

      oldIndex_.resize( numberOfHoles_ );
      newIndex_.resize( numberOfHoles_ );

      size_ = newActSize;

      for( int i=0; i<numberOfHoles_; ++i ) 
      {
        std::cout << "old[ " << i << " ] = " << oldIndex_[ i ];
        std::cout << " new[ " << i << " ] = " << newIndex_[ i ] <<std::endl;
      }

      std::cout << "Size " << size_ << " holes " << numberOfHoles_ << std::endl;
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

    const IndexSetType &indexSet_;
    
    LagrangePointSetMapVectorType &lagrangePointSet_;

    PolyOrderContainerType entityPolynomOrder_; 

    mutable std::vector< DofContainerType* > dofContainer_; 

    int numberOfHoles_ ;
    std::vector< int > oldIndex_ ;
    std::vector< int > newIndex_ ;

    // memory overshoot 
    const double overShoot_ ;

    mutable unsigned int size_;
    unsigned int numDofs_;
  };

} // end namespace Dune 

#endif // #ifndef DUNE_LAGRANGESPACE_MAPPER_HH

#ifndef DUNE_FEM_GENERICPADAPTIVEDOFMAPPER_HH
#define DUNE_FEM_GENERICPADAPTIVEDOFMAPPER_HH

#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/utility/persistentcontainer.hh>

#include <dune/fem/common/forloop.hh>
#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/misc/metaprogramming.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/localkey.hh>
#include <dune/fem/space/lagrange/lagrangepoints.hh>
#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>

namespace Dune
{

  namespace Fem
  {

    // GenericAdaptiveDofMapper
    // ------------------------

    template< class TraitsImp >
    class GenericAdaptiveDofMapper
    : public AdaptiveDofMapper< TraitsImp >
    {
      typedef GenericAdaptiveDofMapper< TraitsImp > ThisType;
      typedef AdaptiveDofMapper< TraitsImp > BaseType;

    public:
      typedef TraitsImp Traits;
      typedef std::size_t SizeType;

      // true if all dofs are associated with entities of codim 0
      static const bool discontinuousMapper = Traits :: discontinuousMapper ;

    protected:
      using BaseType::asImp;

    public:
      //! type of the grid part
      typedef typename Traits::GridPartType GridPartType;

      //! type of entities (codim 0)
      typedef typename Traits::ElementType ElementType;

      //! type of the underlying grid
      typedef typename GridPartType::GridType GridType;

      //! type of the index set
      typedef typename GridPartType::IndexSetType IndexSetType;

      //! type of global key
      typedef typename Traits :: GlobalKeyType  GlobalKeyType ;

      //! type of coordinates within the grid
      typedef typename GridType::ctype FieldType;

      //! dimension of the grid
      static const int dimension = GridType::dimension;

      //! highest codimension used to attach dofs
      static const int highestDimension = ( discontinuousMapper ) ? 0 : dimension ;

      //! order of the Lagrange polynoms
      static const int polynomialOrder = Traits::polynomialOrder;

      //! default partition iterator type for index setup
      static const PartitionIteratorType pitype = GridPartType :: indexSetPartitionType ;

      //! type of vector containing compiled local keys
      typedef typename Traits :: CompiledLocalKeyVectorType  CompiledLocalKeyVectorType;

      //! compiled local key type
      typedef typename CompiledLocalKeyVectorType :: value_type :: value_type  CompiledLocalKeyType;

      //! type of the DoF manager
      typedef DofManager< GridType > DofManagerType;

      enum { minOrder = 1, maxOrder = polynomialOrder, numOrders = maxOrder - minOrder + 1 };

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
          assert( codim <= highestDimension );
          // also for codim == 0 we have more then
          // one storage because of different number of
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
          const unsigned int entry = determineVectorEntry( codim, polOrd );
          assert( entry < dofs_.size() );
          assert( type_ != GeometryType() );
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

      struct PolynomialOrderStorage
      {
      private:
        // PolynomialOrderStorage() : k_( minOrder ), active_( -std::abs(k_) ) {}
        signed char k_; // stores current polynomial order
        signed char active_ ; // stores active/non-active and suggested pol order

      public:
        // PolynomialOrderStorage() : k_( minOrder ), active_( -std::abs(k_) ) {}
        PolynomialOrderStorage( const int k ) : k_( k ), active_( -std::abs(k_) ) {}
        int order () const { return k_; }
        void suggest ( const int k )
        {
          if( active() )
            active_ = std::abs( k );
          else
            active_ = -std::abs( k );
        }
        void set ( const int k ) { k_ = k; active_ = std::abs( k ) ; }
        void activate() { active_ = std::abs( active_ ); }
        bool active () const { return active_ > 0; }
        bool deactivate ( int& k )
        {
          k = k_;
          if( active() )
          {
            active_ = -active_;
            return true ;
          }
          return false;
        }

        int suggested () const { return std::abs( active_ ); }
        void update() { set( suggested() ); }
      };

      typedef PolynomialOrderStorage  PolynomialOrderStorageType;

      typedef PersistentContainer< GridType, EntityDofStorageType > DofContainerType ;

      typedef PersistentContainer< GridType, PolynomialOrderStorageType > PolyOrderContainerType ;
    protected:
      template <int codim, bool dg>
      struct NumDofs
      {
        static int numDofs( const ElementType& entity,
                            const CompiledLocalKeyType& clk,
                            const int subEntity )
        {
          return clk.numDofs( codim, subEntity );
        }
      };

      template <int codim>
      struct NumDofs<codim, true>
      {
        static int numDofs( const ElementType& entity,
                            const CompiledLocalKeyType& clk,
                            const int subEntity )
        {
          if( codim == 0 )
            return clk.size();
          else
            return 0;
        }
      };

      template < int codim >
      struct InsertSubEntities
      {
        static void insertDofs( const ElementType& entity,
                                const CompiledLocalKeyType& clk,
                                const int polOrd,
                                const int subEntity,
                                unsigned int& globalSize,
                                unsigned int& notAlreadyCounted,
                                EntityDofStorage& entityDofs )
        {
          const int numDofs = NumDofs<codim, discontinuousMapper>::numDofs( entity, clk, subEntity );

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

        static void apply( const ElementType& entity,
                           const CompiledLocalKeyType& clk,
                           const int polOrd,
                           unsigned int& globalSize,
                           unsigned int& notAlreadyCounted,
                           std::vector< DofContainerType* > dofContainers )
        {
          DofContainerType &dofContainer = *dofContainers[ codim ];
          if( codim == 0 )
          {
            insertDofs( entity, clk, polOrd, 0, globalSize,
                        notAlreadyCounted, dofContainer[ entity ] );
          }
          else
          {
            const int count = entity.subEntities( codim );
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
        static void apply( const ElementType& entity,
                           const int polOrd,
                           std::vector< DofContainerType* > dofContainers )
        {
          DofContainerType &dofContainer = *dofContainers[ codim ];
          const int count = entity.subEntities( codim );
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
                                 const int order,
                                 CompiledLocalKeyVectorType &compiledLocalKeyVector )
      : gridPart_( gridPart ),
        dm_( DofManagerType :: instance(gridPart.grid()) ),
        compiledLocalKeys_( compiledLocalKeyVector ),
        order_( ( order >= minOrder && order <= maxOrder ) ? order : maxOrder ),
        entityPolynomOrder_( gridPart.grid(), 0, PolynomialOrderStorageType( order_ ) ),
        dofContainer_( dimension+1, nullptr ),
        numberOfHoles_( 0 ),
        oldIndex_(),
        newIndex_(),
        size_(0),
        maxNumDofs_( 0 ),
        sequence_( dm_.sequence() )
      {
        for( int codim = 0; codim <= highestDimension; ++codim )
          dofContainer_[ codim ] = new DofContainerType( gridPart.grid(), codim );

        for( size_t i=0; i<compiledLocalKeys_.size(); ++i )
        {
          maxNumDofs_ = std :: max( maxNumDofs_, compiledLocalKeys_[ i ].maxSize() );
        }

        resize();
        // register to DofManager for adaptation
        dm_.addIndexSet( asImp() );
      }

      //! sort of copy constructor
      GenericAdaptiveDofMapper ( const GenericAdaptiveDofMapper& other,
                                 const int order,
                                 CompiledLocalKeyVectorType &compiledLocalKeyVector )
      : gridPart_( other.gridPart_ ),
        dm_( other.dm_ ),
        compiledLocalKeys_( compiledLocalKeyVector ),
        order_( ( order >= minOrder && order <= maxOrder ) ? order : maxOrder ),
        entityPolynomOrder_( other.entityPolynomOrder_ ),
        dofContainer_( dimension+1, nullptr ),
        numberOfHoles_( other.numberOfHoles_ ),
        oldIndex_( other.oldIndex_ ),
        newIndex_( other.newIndex_ ),
        size_( other.size_ ),
        maxNumDofs_( other.maxNumDofs_ ),
        sequence_( other.sequence_ )
      {
        for( int codim = 0; codim <= highestDimension; ++codim )
          dofContainer_[ codim ] = new DofContainerType( *(other.dofContainer_[ codim ]) );

        dm_.addIndexSet( asImp() );
      }

      int polynomOrder( const ElementType& entity ) const
      {
        return entityPolynomOrder_[ entity ].order();
      }

      int suggestedOrder( const ElementType& entity ) const
      {
        return entityPolynomOrder_[ entity ].suggestedOrder();
      }

      void suggestPolynomOrder( const ElementType& entity, const int polOrd )
      {
        // minOrder is static but order_ is dynamically set
        if( polOrd < minOrder || polOrd > order_ )
          return ;

        entityPolynomOrder_[ entity ].suggest( polOrd );
      }

      DofContainerType &dofContainer ( const std::size_t codim ) const
      {
        assert( codim < dofContainer_.size() );
        assert( dofContainer_[ codim ] );
        return *(dofContainer_[ codim ]);
      }

      const CompiledLocalKeyType&
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

      template< class Functor >
      void mapEach ( const ElementType &element, Functor f ) const
      {
        const int n = numDofs( element );
        for( int i = 0; i < n; ++i )
          f( i, mapToGlobal( element, i ) );
      }

      /** \copydoc Dune::Fem::DofMapper::mapToGlobal */
      int mapToGlobal ( const ElementType &entity, const int localDof ) const
      {
        if( discontinuousMapper )
        {
          // get polynomial order
          const int polOrd = polynomOrder( entity );

          // return global dof
          return dofContainer( 0 )[ entity ].dof( 0, polOrd, localDof );
        }
        else
        {
          // the continuous case
          const int polOrd = polynomOrder( entity );

          const CompiledLocalKeyType& compLocalKey = compiledLocalKey( polOrd, entity.type() );
          // get dof info for entity and local dof
          const Fem::LocalKey &dofInfo = compLocalKey.localKey( localDof );

          const unsigned int codim = dofInfo.codim();
          const unsigned int subEntity = dofInfo.subEntity();

          unsigned int index = dofInfo.index() ;

          // account for possible twists in the grid (only 2d)
          if( dimension == 2 && codim == 1 )
          {
            auto refElem = referenceElement< FieldType, dimension >( entity.type() );

#ifndef NDEBUG
            const int vxSize = refElem.size( subEntity, codim, dimension );
            // two vertices per edge in 2d
            assert( vxSize == 2 );
#endif
            const int vx[ 2 ] = { refElem.subEntity ( subEntity, codim, 0, dimension ),
                                  refElem.subEntity ( subEntity, codim, 1, dimension ) };

            // flip index if face is twisted
            if( gridPart_.grid().localIdSet().subId( entity, vx[ 0 ], dimension ) >
                gridPart_.grid().localIdSet().subId( entity, vx[ 1 ], dimension ) )
            {
              const unsigned int numDofsSubEntity = compLocalKey.numDofs( codim, subEntity );
              index = numDofsSubEntity - index - 1;
            }
          }

          assert( index < compLocalKey.numDofs( codim, subEntity ) );
          return dofContainer( codim )( entity, subEntity ).dof( codim, polOrd, index );
        }
      }

      /** \copydoc Dune::Fem::DofMapper::mapEachEntityDof */
      template< class Entity, class Functor >
      void mapEachEntityDof ( const Entity &entity, Functor f ) const
      {
        const int n = numEntityDofs( entity );
        for( int i = 0; i < n; ++i )
          f( i, dofContainer( Entity::codimension )[ entity ].entityDof( i ) );
      }

      void onSubEntity ( const ElementType &element, int i, int c, std::vector< bool > &indices ) const
      {
        indices.resize( numDofs(element) );
        if( discontinuousMapper )
        {
          if (c == 0)
            std::fill(indices.begin(),indices.end(),true);
          else
            std::fill(indices.begin(),indices.end(),false);
        }
        else
        {
          DUNE_THROW( NotImplemented, "Method onSubEntity(...) not yet implemented for TupleMapper" );
        }
      }

      void map ( const ElementType &element, std::vector< SizeType > &indices ) const
      {
        indices.resize( numDofs( element ) );
        mapEach( element, AssignFunctor< std::vector< SizeType > >( indices ) );
      }

      template< class Entity >
      void mapEntityDofs ( const Entity &entity, std::vector< SizeType > &indices ) const
      {
        indices.resize( numEntityDofs( entity ) );
        mapEachEntityDof( entity, AssignFunctor< std::vector< SizeType > >( indices ) );
      }

      /** \copydoc Dune::Fem::DofMapper::maxNumDofs() const */
      int maxNumDofs () const
      {
        return maxNumDofs_;
      }

      /** \copydoc Dune::Fem::DofMapper::numDofs(const ElementType &entity) const */
      int numDofs ( const ElementType &entity ) const
      {
        const int polOrd = polynomOrder( entity );
        return compiledLocalKey( polOrd, entity.type() ).size();
      }

      /** \copydoc Dune::Fem::DofMapper::numEntityDofs(const Entity &entity) const */
      template< class Entity >
      int numEntityDofs ( const Entity &entity ) const
      {
        if( discontinuousMapper )
        {
          if( Entity :: codimension == 0 )
            return dofContainer( 0 )[ entity ].entityDofs();
          else
            return 0;
        }
        else
        {
          // !!! to be revised
          // This implementation only works for nonhybrid grids (simplices or cubes)
          return dofContainer( Entity :: codimension )[ entity ].entityDofs();
        }
      }

      /** \brief Check, whether any DoFs are associated with a codimension */
      bool contains ( int codim ) const
      {
        return (discontinuousMapper) ? (codim == 0) : true;
      }

      /** \brief Check, whether the data in a codimension has fixed size */
      bool fixedDataSize ( int codim ) const
      {
        return false;
      }

      /** \copydoc Dune::Fem::DofMapper::oldIndex */
      int oldIndex ( const int hole, const int block ) const
      {
        assert( oldIndex_[ hole ] >= 0 );
        return oldIndex_[ hole ];
      }

      /** \copydoc Dune::Fem::DofMapper::newIndex */
      int newIndex ( const int hole, const int block ) const
      {
        assert( newIndex_[ hole ] >= 0 );
        return newIndex_[ hole ];
      }

      /** \copydoc Dune::Fem::DofMapper::numberOfHoles */
      int numberOfHoles ( const int block ) const
      {
        return numberOfHoles_;
      }

      /** \copydoc Dune::Fem::DofMapper::numBlocks
       */
      int numBlocks () const
      {
        return 1;
      }

      /** \copydoc Dune::Fem::DofMapper::oldOffset
       */
      int oldOffSet ( const int block ) const
      {
        return 0;
      }

      /** \copydoc Dune::Fem::DofMapper::newOffset
       */
      int offSet ( const int block ) const
      {
        return 0;
      }

      /** \copydoc Dune::Fem::DofMapper::consecutive() const */
      bool consecutive () const
      {
        return true;
      }

      // Adaptation Methods (as for Index Sets)
      void resizeContainers()
      {
        entityPolynomOrder_.resize( PolynomialOrderStorageType( order_ ) );
        entityPolynomOrder_.shrinkToFit();
        for( int codim = 0; codim <= highestDimension; ++codim )
        {
          dofContainer( codim ).resize();
          dofContainer( codim ).shrinkToFit();
        }
      }

      void insertEntity ( const ElementType &entity )
      {
        resizeContainers();
        insertEntityDofs( entity );
      }

      // return number of local dofs that were not visited yet
      unsigned int insertEntityDofs( const ElementType &entity )
      {
        PolynomialOrderStorageType& polyStorage = entityPolynomOrder_[ entity ];
        if( ! polyStorage.active() )
        {
          unsigned int notAlreadyCounted = 0;

          const int polOrd = polyStorage.order();
          // get lagrange point set
          const CompiledLocalKeyType& clk = compiledLocalKey( polOrd, entity.type() );

          //std::cout << "Insert Entity " << gridPart_.grid().localIdSet().id( entity ) << std::endl;

          // activate dofs for this entity
          polyStorage.activate();

          // insert for all sub entities
          Fem::ForLoop< InsertSubEntities, 0, highestDimension>::
            apply( entity, clk, polOrd, size_, notAlreadyCounted, dofContainer_ );

          //printEntityDofs( entity );
          return notAlreadyCounted ;
        }

        return 0;
      }

      void removeEntity ( const ElementType &entity )
      {
        int polOrd;
        // polOrd ist set on call of deactivate
        if( entityPolynomOrder_[ entity ].deactivate( polOrd ) )
        {
          Fem::ForLoop< RemoveSubEntities, 0, highestDimension>::
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
        // set new polynomial orders for entities
        for( auto& pol : entityPolynomOrder_ )
          pol.update();

        sequence_ = -1;
        compress();
      }

      // insert father element when conforming refinement is enabled
      // (due to refinement of more than one level)
      unsigned int insertFather( const ElementType &entity )
      {
        if( entity.hasFather() )
        {
          // if father is a new element, insert it
          ElementType dad = entity.father();
          if( dad.isNew() )
          {
            unsigned int usedSize = insertEntityDofs( dad );
            // also insert dad's fathers
            usedSize += insertFather( dad );
            return usedSize;
          }
        }
        return 0;
      }

      //! return true if elements can be refined more than once during adaptation
      bool considerHierarchy () const
      {
        return DGFGridInfo< GridType > :: refineStepsForHalf() > 1 ||
               ! Dune::Capabilities::hasSingleGeometryType<GridType> :: v  ;
      }

      //! return number of DoFs currently used for space
      size_t insertAllUsed()
      {
        // reset all dof entries
        setUnused();

        // count current size
        size_t usedSize = 0;

        const bool considerHierarchyOfElements = considerHierarchy();

        // pitype is determined by the GridPart, see above
        const auto end = gridPart_.template end<0, pitype>();
        for( auto it = gridPart_.template begin<0, pitype>();
             it != end ; ++it )
        {
          const ElementType &element = *it;
          if( considerHierarchyOfElements )
          {
            // insert father elements (conforming and hybrid grids only)
            usedSize += insertFather( element );
          }

          // number of new dofs (not already counted) is returned
          usedSize += insertEntityDofs( element );
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

      void printEntityDofs( const ElementType& entity ) const
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
        // deactivate all entries in the polynomial order container
        {
          typedef typename PolyOrderContainerType :: Iterator Iterator;
          const Iterator endit = entityPolynomOrder_.end();
          for( Iterator it = entityPolynomOrder_.begin(); it != endit; ++it )
          {
            PolynomialOrderStorageType& p = *it;
            int pOrd;
            p.deactivate( pOrd );
          }
        }

        for( int codim = 0; codim <= highestDimension; ++codim )
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

        // adjust size of containers
        resizeContainers();

        // make all dofs that are currently used
        // (corresponding to GridPart's iterator)
        const size_t usedSize = insertAllUsed();

        //std::cout << "Size " << size_ << " holes " << numberOfHoles_ << std::endl;
        //printDofs();

        // reset number of holes
        numberOfHoles_ = 0;

        // true if a least one dof must be copied
        bool haveToCopy = false;

        std::vector<int> validHoles;

        {
          // default is true which means entry is hole
          std::vector< bool > holeMarker( size_, true );

          // mark holes
          for( int codim = 0; codim <= highestDimension; ++codim )
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
          for( size_t i=0; i<usedSize; ++i )
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

          for( int codim = 0; codim <= highestDimension; ++codim )
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

      void backup () const
      {}

      void restore ()
      {}

      template< class InStream >
      void read ( InStream &in )
      {}

      template< class OutStream >
      void write ( OutStream &out )
      {}

      GenericAdaptiveDofMapper ( const ThisType& ) = delete;
      ThisType& operator=( const ThisType& ) = delete;

    private:

      const GridPartType& gridPart_;
      // reference to dof manager
      DofManagerType& dm_;

      CompiledLocalKeyVectorType &compiledLocalKeys_;

      // default polynomial order to be set to newly inserted elements
      const int order_;

      PolyOrderContainerType entityPolynomOrder_;

      mutable std::vector< DofContainerType* > dofContainer_;

      int numberOfHoles_ ;
      std::vector< int > oldIndex_ ;
      std::vector< int > newIndex_ ;

      mutable unsigned int size_;
      unsigned int maxNumDofs_;
      int sequence_ ;
    }; // class GenericAdaptiveDofMapper


    namespace Capabilities
    {
      // isConsecutiveIndexSet
      // ---------------------

      template< class Traits >
      struct isConsecutiveIndexSet< GenericAdaptiveDofMapper< Traits > >
      {
        static const bool v = true;
      };

      template< class Traits >
      struct isAdaptiveDofMapper< GenericAdaptiveDofMapper< Traits > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GENERICPADAPTIVEDOFMAPPER_HH

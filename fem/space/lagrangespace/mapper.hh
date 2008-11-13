#ifndef DUNE_LAGRANGESPACE_MAPPER_HH
#define DUNE_LAGRANGESPACE_MAPPER_HH

//- Dune includes 
#include <dune/common/geometrytype.hh>
#include <dune/common/exceptions.hh>

//- Dune-Fem includes 
#include <dune/fem/misc/codimmap.hh>
#include <dune/fem/misc/gridhelper.hh>
#include <dune/fem/misc/metaprogramming.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/dofmapper.hh>

//- local includes 
#include "lagrangepoints.hh"

namespace Dune
{
  
  template< class GridPart, unsigned int polOrder, unsigned int dimR >
  class LagrangeMapper;


  
  template< class GridPart, unsigned int polOrder, unsigned int dimR >
  struct LagrangeMapperTraits
  {
    typedef GridPart GridPartType;
    
    enum { polynomialOrder = polOrder };

    //! dimension of the discrete function space's range
    enum { dimRange = dimR };

    typedef typename GridPartType :: template Codim< 0 > :: IteratorType :: Entity
      EntityType;

    typedef LagrangeMapper< GridPartType, polynomialOrder, dimRange >
      DofMapperType;

    typedef DefaultDofMapIterator< EntityType, DofMapperType >
      DofMapIteratorType;
  };


  
  template< class GridPart, unsigned int dimR >
  class LagrangeMapper< GridPart, 1, dimR >
  : public DofMapperDefault< LagrangeMapperTraits< GridPart, 1, dimR > >
  {
  public:
    typedef LagrangeMapperTraits< GridPart, 1, dimR > Traits;
    
    //! type of the grid part
    typedef typename Traits :: GridPartType GridPartType;

    //! type of entities (codim 0)
    typedef typename Traits :: EntityType EntityType;

    //! type of DofMapIterator
    typedef typename Traits :: DofMapIteratorType DofMapIteratorType;
    
    //! type of the underlying grid
    typedef typename GridPartType :: GridType GridType;

    //! type of coordinates within the grid
    typedef typename GridType :: ctype FieldType;

    //! dimension of the grid
    enum { dimension = GridType :: dimension };

    //! order of the Lagrange polynoms
    enum { polynomialOrder = Traits :: polynomialOrder };

    //! dimension of the discrete function space's range
    enum { dimRange = Traits :: dimRange };

  private:
    typedef LagrangeMapper< GridPartType, polynomialOrder, dimRange > ThisType;
    typedef DofMapperDefault< Traits > BaseType;

  public:
    //! type of the index set
    typedef typename GridPartType :: IndexSetType IndexSetType;

    //! type of the Lagrange point set
    typedef LagrangePointSet< GridPartType, polynomialOrder >
      LagrangePointSetType;
    //! type of the map for the Lagrange point sets
    typedef std :: map< const GeometryType, const LagrangePointSetType* >
      LagrangePointSetMapType;

  private:
    const IndexSetType &indexSet_;
    unsigned int maxDofs_;

  public:
    //! constructor
    LagrangeMapper ( const GridPartType &gridPart,
                     LagrangePointSetMapType &lagrangePointSet )
    : indexSet_( gridPart.indexSet() )
    , maxDofs_ ( 0 )
    {
      typedef typename LagrangePointSetMapType :: iterator IteratorType;
      IteratorType end = lagrangePointSet.end();
      for( IteratorType it = lagrangePointSet.begin(); it != end; ++it )
      {
        const LagrangePointSetType *set = (*it).second;
        if( set == NULL )
          continue;
        
        const unsigned int setDofs = set->numDofs( dimension );
        maxDofs_ = (maxDofs_ >= setDofs) ? maxDofs_ : setDofs;
      }
    }
   
    //! destructor
    virtual ~LagrangeMapper ()
    {
    }

    /** \copydoc Dune::DofMapper::size() const */
    int size () const
    {
      return dimRange * indexSet_.size( dimension );
    }

    /** \copydoc Dune::DofMapper::begin(const EntityType &entity) const */
    inline DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: beginIterator, entity, *this );
    }
    
    /** \copydoc Dune::DofMapper::end(const EntityType &entity) const */
    inline DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: endIterator, entity, *this );
    }

    /** \copydoc Dune::DofMapper::mapToGlobal */
    inline int mapToGlobal ( const EntityType &entity,
                             const int localDof ) const
    {
      const int coordinate = localDof % dimRange;
      const int localPoint = localDof / dimRange;
      const int globalPoint
        = indexSet_.template subIndex< dimension >( entity, localPoint );
      return dimRange * globalPoint + coordinate;
    }

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal */
    template< class Entity >
    inline int mapEntityDofToGlobal ( const Entity &entity,
                                      const int localDof ) const
    {
      if( Entity :: codimension != (unsigned int) dimension )
        DUNE_THROW( RangeError, "No such local DoF." );

      assert( (localDof >= 0) && (localDof < dimRange) );
      return dimRange * indexSet_.index( entity ) + localDof;
    }
    
    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    int maxNumDofs () const
    {
      return dimRange * maxDofs_;
    }

    using BaseType :: numDofs;
    
    /** \copydoc Dune::DofMapper::numDofs(const EntityType &entity) const */
    int numDofs ( const EntityType &entity ) const
    {
      return dimRange * entity.template count< dimension >();
    }

    template< class Entity >
    inline int numEntityDofs ( const Entity &entity ) const
    {
      return (Entity :: codimension == (unsigned int) dimension ? dimRange : 0);
    }

    /** \brief Check, whether any DoFs are associated with a codimension */
    inline bool contains ( unsigned int codim ) const
    {
      return (codim == dimension);
    }

    /** \brief Check, whether the data in a codimension has fixed size */
    inline bool fixedDataSize ( unsigned int codim ) const
    {
      return true;
    }
   
    /** \copydoc Dune::DofMapper::oldIndex */
    int oldIndex ( int hole, int ) const
    {
      const int coordinate = hole % dimRange;
      const int setHole = hole / dimRange;
      const int setIndex = indexSet_.oldIndex( setHole, dimension );
      return setIndex * dimRange + coordinate;
    }

    /** \copydoc Dune::DofMapper::newIndex */
    int newIndex ( int hole , int ) const
    {
      const int coordinate = hole % dimRange;
      const int setHole = hole / dimRange;
      const int setIndex = indexSet_.newIndex( setHole, dimension );
      return setIndex * dimRange + coordinate;
    }

    /** \copydoc Dune::DofMapper::numberOfHoles
     */
    int numberOfHoles ( int ) const
    {
      return dimRange * indexSet_.numberOfHoles( dimension );
    }

    /** \copydoc Dune::DofMapper::newSize
     */
    int newSize () const
    {
      return this->size();
    }

    /** \copydoc Dune::DofMapper::consecutive() const */
    bool consecutive () const
    {
      return BaseType :: checkConsecutive( indexSet_ );
    }
  };



  template< class GridPart, unsigned int dimR >
  class LagrangeMapper< GridPart, 2, dimR >
  : public DofMapperDefault< LagrangeMapperTraits< GridPart, 2, dimR > >
  {
  public:
    typedef LagrangeMapperTraits< GridPart, 2, dimR > Traits;
    
    //! type of the grid part
    typedef typename Traits :: GridPartType GridPartType;

    //! type of entities (codim 0)
    typedef typename Traits :: EntityType EntityType;

    //! type of DofMapIterator
    typedef typename Traits :: DofMapIteratorType DofMapIteratorType;
 
    //! type of the underlying grid
    typedef typename GridPartType :: GridType GridType;

    //! type of coordinates within the grid
    typedef typename GridType :: ctype FieldType;

    //! dimension of the grid
    enum { dimension = GridType :: dimension };

    //! order of the Lagrange polynoms
    enum { polynomialOrder = Traits :: polynomialOrder };

    //! dimension of the discrete function space's range
    enum { dimRange = Traits :: dimRange };

  private:
    typedef LagrangeMapper< GridPartType, polynomialOrder, dimRange > ThisType;
    typedef DofMapperDefault< Traits > BaseType;

  public:
    //! type of the index set
    typedef typename GridPartType :: IndexSetType IndexSetType;

    //! type of the Lagrange point set
    typedef LagrangePointSet< GridPartType, polynomialOrder >
      LagrangePointSetType;
    //! type of the map for the Lagrange point sets
    typedef std :: map< const GeometryType, const LagrangePointSetType* >
      LagrangePointSetMapType;

    //! type of the DoF manager
    typedef DofManager< GridType > DofManagerType;
    //! type of the DoF manager factory
    typedef DofManagerFactory< DofManagerType > DMFactoryType;

  protected:
    typedef typename LagrangePointSetType :: DofInfo DofInfo;

  private:
    struct CodimCallInterface;
    
    template< unsigned int codim >
    struct CodimCall;

    typedef CodimMap< dimension+1, CodimCall > CodimCallMapType;

  private:
    // reference to dof manager needed for debug issues 
    const DofManagerType& dm_;

    const IndexSetType &indexSet_;
    const CodimCallMapType codimCall_;
    
    LagrangePointSetMapType &lagrangePointSet_;

    unsigned int maxDofs_[ dimension+1 ];
    mutable unsigned int offset_[ dimension+1 ];
    mutable unsigned int oldOffSet_[ dimension+1 ];
    mutable unsigned int size_;
    unsigned int numDofs_;

    // for debugging only 
    mutable int sequence_;

  public:
    //! constructor
    LagrangeMapper ( const GridPartType &gridPart,
                     LagrangePointSetMapType &lagrangePointSet )
    : dm_( DMFactoryType :: getDofManager(gridPart.grid()) ),
      indexSet_( gridPart.indexSet() ),
      lagrangePointSet_( lagrangePointSet ),
      sequence_( dm_.sequence() )
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
        
        numDofs_ = std :: max( numDofs_, set->size() );
        for( int codim = 0; codim <= dimension; ++codim )
          maxDofs_[ codim ]
            = std :: max( maxDofs_[ codim ], set->maxDofs( codim ) );
      }

      size_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
      {
        offset_[ codim ] = size_;
        oldOffSet_[ codim ] = size_;
        size_ += indexSet_.size( codim ) * maxDofs_[ codim ];
      }
    }
    
    //! destructor 
    virtual ~LagrangeMapper ()
    {}

    //! return overall number of degrees of freedom 
    int size () const
    {
      return dimRange * size_;
    }

    /** \copydoc Dune::DofMapperInterface::begin(const EntityType &entity) const */
    inline DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: beginIterator, entity, *this );
    }
    
    /** \copydoc Dune::DofMapperInterface::end(const EntityType &entity) const */
    inline DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: endIterator, entity, *this );
    }

    /** \copydoc Dune::DofMapperInterface::mapToGlobal */
    int mapToGlobal ( const EntityType &entity, const int local ) const
    {
      const int coordinate = local % dimRange;
      const int localDof = local / dimRange;
      
      // unsigned int codim, subEntity;
      const LagrangePointSetType *set
        = lagrangePointSet_[ entity.geometry().type() ];
      const DofInfo& dofInfo = set->dofInfo( localDof );

      const unsigned int codim = dofInfo.codim;
      const int subIndex
        = codimCall_[ codim ].subIndex( *this, entity, dofInfo.subEntity );

      const int globalDof
        = dimRange * (offset_[ codim ] + subIndex) + coordinate;

      assert( codimCall_[ codim ].checkMapEntityDofToGlobal
                ( *this, entity, dofInfo.subEntity, coordinate, globalDof ) );
      return globalDof;
    }

    /** \copydoc Dune::DofMapperInterface::mapEntityDofToGlobal */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const 
    {
      assert( localDof < numEntityDofs( entity ) );
      const int globalDofPt
        = offset_[ Entity :: codimension ] + indexSet_.index( entity );
      return dimRange * globalDofPt + localDof;
    }
    
    /** \copydoc Dune::DofMapperInterface::maxNumDofs() const */
    inline int maxNumDofs () const
    {
      return dimRange * numDofs_;
    }

    using BaseType :: numDofs;

    /** \copydoc Dune::DofMapperInterface::numDofs(const EntityType &entity) const */
    inline int numDofs ( const EntityType &entity ) const
    {
      return dimRange * lagrangePointSet_[ entity.geometry().type() ]->size();
    }

    /** \copydoc Dune::DofMapperInterface::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    inline int numEntityDofs ( const Entity &entity ) const
    {
      // This implementation only works for nonhybrid grids (simplices or cubes)
      return dimRange * maxDofs_[ Entity :: codimension ];
    }
    
    /** \brief Check, whether any DoFs are associated with a codimension */
    inline bool contains ( unsigned int codim ) const
    {
      return true;
    }

    /** \brief Check, whether the data in a codimension has fixed size */
    inline bool fixedDataSize ( unsigned int codim ) const
    {
      return false;
    }

    /** \copydoc Dune::DofMapperInterface::oldIndex */
    int oldIndex ( const int num, const int codim ) const
    {
      // corresponding number of set is newn 
      const int newn  = static_cast<int> (num / dimRange);
      // local number of dof is local 
      const int local = (num % dimRange);
      // codim to be revised 
      return dimRange * (oldOffSet_[codim] + indexSet_.oldIndex(newn,codim)) + local;
    }

    /** \copydoc Dune::DofMapperInterface::newIndex */
    int newIndex ( const int num , const int codim) const
    {
      // corresponding number of set is newn 
      const int newn  = static_cast<int> (num / dimRange);
      // local number of dof is local 
      const int local = (num % dimRange);
      // codim to be revised 
      return dimRange * (offset_[codim] + indexSet_.newIndex(newn,codim)) + local;
    }

    /** \copydoc Dune::DofMapperInterface::numberOfHoles */
    int numberOfHoles ( const int codim ) const
    {
      return (maxDofs_[ codim ] > 0) ? 
        (dimRange * indexSet_.numberOfHoles( codim )) : 0;
    }

    /** \copydoc Dune::DofMapperInterface::update */
    void update()
    {
      // assure that update is only called once per 
      // dof manager resize or compress 
      if( sequence_ != dm_.sequence() )
      {
        // calculate new size 
        size_ = 0;
        for( int codim = 0; codim <= dimension; ++codim )
        {
          oldOffSet_[ codim ] = offset_[ codim ];
          offset_[ codim ] = size_;
          size_ += indexSet_.size( codim ) * maxDofs_[ codim ];
        }
        sequence_ = dm_.sequence();
      }
    }

    /** \copydoc Dune::DofMapperInterface::numBlocks
     */
    int numBlocks () const 
    {
      return dimension + 1;
    }

    /** \copydoc Dune::DofMapperInterface::oldOffset
     */
    int oldOffSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return dimRange * oldOffSet_[ block ];
    }

    /** \copydoc Dune::DofMapperInterface::newOffset
     */
    int offSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return dimRange * offset_[ block ];
    }

    /** \copydoc Dune::DofMapperInterface::newSize
     */ 
    int newSize () const
    {
      int newSize = 0;
      for( int codim = 0; codim <= dimension; ++codim )
        newSize += indexSet_.size( codim ) * maxDofs_[ codim ];
      return dimRange * newSize;
    }

    /** \copydoc Dune::DofMapper::consecutive() const */
    bool consecutive () const
    {
      return BaseType :: checkConsecutive( indexSet_ );
    }
  };


/** \cond */
  template< class GridPart, unsigned int dimR >
  struct LagrangeMapper< GridPart, 2, dimR > :: CodimCallInterface
  {
    typedef LagrangeMapper< GridPart, 2, dimR > MapperType;

    virtual ~CodimCallInterface ()
    {}

    virtual int subIndex ( const MapperType &mapper,
                           const EntityType &entity,
                           int i ) const = 0;

    virtual bool
    checkMapEntityDofToGlobal ( const MapperType &mapper,
                                const EntityType &entity,
                                const int subEntity,
                                const int localDof,
                                const int globalDof ) const = 0;
  };
/** \endcond */



/** \cond */
  template< class GridPart, unsigned int dimR >
  template< unsigned int codim >
  struct LagrangeMapper< GridPart, 2, dimR > :: CodimCall
  : public CodimCallInterface
  {
    typedef LagrangeMapper< GridPart, 2, dimR > MapperType;

    typedef CodimCallInterface BaseType;

    virtual int subIndex ( const MapperType &mapper,
                           const EntityType &entity,
                           int i ) const
    {
      return mapper.indexSet_.template subIndex< codim >( entity, i );
    }

    virtual bool checkMapEntityDofToGlobal ( const MapperType &mapper,
                                             const EntityType &entity,
                                             const int subEntity,
                                             const int localDof,
                                             const int globalDof ) const
    {
      typedef Capabilities :: hasEntity< GridType, codim > HasEntity;

      return checkMapEntityDofToGlobal
        ( mapper, entity, subEntity, localDof, globalDof,
          MetaBool< HasEntity :: v >() );
    }

    inline bool
    checkMapEntityDofToGlobal ( const MapperType &mapper,
                                const EntityType &entity,
                                const int subEntity,
                                const int localDof,
                                const int globalDof,
                                const MetaBool< true > hasEntity ) const
    {
      typedef typename EntityType :: template Codim< codim > :: EntityPointer
        SubEntityPtrType;

      const SubEntityPtrType subEntityPtr
        = entity.template entity< codim >( subEntity );
      const int globalEntityDof
        = mapper.mapEntityDofToGlobal( *subEntityPtr, localDof );
      return (globalEntityDof == globalDof);
    }

    inline bool
    checkMapEntityDofToGlobal ( const MapperType &mapper,
                                const EntityType &entity,
                                const int subEntity,
                                const int localDof,
                                const int globalDof,
                                const MetaBool< false > hasEntity ) const
    {
      return true;
    }
  };
  /** \endcond */



  // Higher Order Lagrange Mapper
  // ----------------------------
  //
  // Note: This mapper assumes that the grid is "twist-free".

#ifdef USE_TWISTFREE_MAPPER
  template< class GridPart, unsigned int polOrder, unsigned int dimR >
  class LagrangeMapper
  : public DofMapperDefault< LagrangeMapperTraits< GridPart, polOrder, dimR > >
  {
  public:
    typedef LagrangeMapperTraits< GridPart, polOrder, dimR > Traits;
    
    //! type of the grid part
    typedef typename Traits :: GridPartType GridPartType;

    //! type of entities (codim 0)
    typedef typename Traits :: EntityType EntityType;

    //! type of DofMapIterator
    typedef typename Traits :: DofMapIteratorType DofMapIteratorType;
 
    //! type of the underlying grid
    typedef typename GridPartType :: GridType GridType;

    //! type of coordinates within the grid
    typedef typename GridType :: ctype FieldType;

    //! dimension of the grid
    enum { dimension = GridType :: dimension };

    //! order of the Lagrange polynoms
    enum { polynomialOrder = Traits :: polynomialOrder };

    //! dimension of the discrete function space's range
    enum { dimRange = Traits :: dimRange };

  private:
    typedef LagrangeMapper< GridPartType, polynomialOrder, dimRange > ThisType;
    typedef DofMapperDefault< Traits > BaseType;

  public:
    //! type of the index set
    typedef typename GridPartType :: IndexSetType IndexSetType;

    //! type of the Lagrange point set
    typedef LagrangePointSet< GridPartType, polynomialOrder >
      LagrangePointSetType;
    //! type of the map for the Lagrange point sets
    typedef std :: map< const GeometryType, const LagrangePointSetType* >
      LagrangePointSetMapType;

    //! type of the DoF manager
    typedef DofManager< GridType > DofManagerType;
    //! type of the DoF manager factory
    typedef DofManagerFactory< DofManagerType > DMFactoryType;

  protected:
    typedef typename LagrangePointSetType :: DofInfo DofInfo;

  private:
    struct CodimCallInterface;
    
    template< unsigned int codim >
    struct CodimCall;

    typedef CodimMap< dimension+1, CodimCall > CodimCallMapType;

  private:
    // reference to dof manager needed for debug issues 
    const DofManagerType& dm_;

    const IndexSetType &indexSet_;
    const CodimCallMapType codimCall_;
    
    LagrangePointSetMapType &lagrangePointSet_;

    unsigned int maxDofs_[ dimension+1 ];
    mutable unsigned int offset_[ dimension+1 ];
    mutable unsigned int oldOffSet_[ dimension+1 ];
    mutable unsigned int size_;
    unsigned int numDofs_;

    // for debugging only 
    mutable int sequence_;

  public:
    //! constructor
    LagrangeMapper ( const GridPartType &gridPart,
                     LagrangePointSetMapType &lagrangePointSet )
    : dm_( DMFactoryType :: getDofManager(gridPart.grid()) ),
      indexSet_( gridPart.indexSet() ),
      lagrangePointSet_( lagrangePointSet ),
      sequence_( dm_.sequence() )
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

        numDofs_ = std :: max( numDofs_, set->size() );
        for( int codim = 0; codim <= dimension; ++codim )
          maxDofs_[ codim ]
            = std :: max( maxDofs_[ codim ], set->maxDofs( codim ) );
      }

      size_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
      {
        offset_[ codim ] = size_;
        oldOffSet_[ codim ] = size_;
        size_ += indexSet_.size( codim ) * maxDofs_[ codim ];
      }
    }
    
    //! destructor 
    virtual ~LagrangeMapper ()
    {}

    //! return overall number of degrees of freedom 
    int size () const
    {
      return dimRange * size_;
    }

    /** \copydoc Dune::DofMapper::begin(const EntityType &entity) const */
    inline DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: beginIterator, entity, *this );
    }
    
    /** \copydoc Dune::DofMapper::end(const EntityType &entity) const */
    inline DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: endIterator, entity, *this );
    }

    /** \copydoc Dune::DofMapper::mapToGlobal */
    int mapToGlobal ( const EntityType &entity, const int local ) const
    {
      const int coordinate = local % dimRange;
      const int localDof = local / dimRange;
      
      // unsigned int codim, subEntity;
      const LagrangePointSetType *set
        = lagrangePointSet_[ entity.geometry().type() ];
      const DofInfo& dofInfo = set->dofInfo( localDof );
      
      const int entityDof = dofInfo.dofNumber * dimRange + coordinate;

      const unsigned int codim = dofInfo.codim;
      const int subIndex
        = codimCall_[ codim ].subIndex( *this, entity, dofInfo.subEntity );

      const int globalDof
        = dimRange * (offset_[ codim ] + subIndex * maxDofs_[ codim ])
          + entityDof;

      assert( codimCall_[ codim ].checkMapEntityDofToGlobal
                ( *this, entity, dofInfo.subEntity, entityDof, globalDof ) );
      return globalDof;
    }

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const 
    {
      const unsigned int codim = Entity :: codimension;

      assert( localDof < numEntityDofs( entity ) );
      const int offset = offset_[ codim ]
                         + indexSet_.index( entity ) * maxDofs_[ codim ];
      return dimRange * offset + localDof;
    }
    
    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    inline int maxNumDofs () const
    {
      return dimRange * numDofs_;
    }

    using BaseType :: numDofs;

    /** \copydoc Dune::DofMapper::numDofs(const EntityType &entity) const */
    inline int numDofs ( const EntityType &entity ) const
    {
      return dimRange * lagrangePointSet_[ entity.geometry().type() ]->size();
    }

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    inline int numEntityDofs ( const Entity &entity ) const
    {
      // This implementation only works for nonhybrid grids (simplices or cubes)
      return dimRange * maxDofs_[ Entity :: codimension ];
    }
    
    /** \brief Check, whether any DoFs are associated with a codimension */
    inline bool contains ( unsigned int codim ) const
    {
      return true;
    }

    /** \brief Check, whether the data in a codimension has fixed size */
    inline bool fixedDataSize ( unsigned int codim ) const
    {
      return false;
    }

    /** \copydoc Dune::DofMapper::oldIndex */
    int oldIndex ( const int num, const int codim ) const
    {
      // corresponding number of set is newn 
      const int newn  = static_cast<int> (num / dimRange);
      // local number of dof is local 
      const int local = (num % dimRange);
      // codim to be revised 
      return dimRange * (oldOffSet_[codim] + indexSet_.oldIndex(newn,codim)) + local;
    }

    /** \copydoc Dune::DofMapper::newIndex */
    int newIndex ( const int num , const int codim) const
    {
      // corresponding number of set is newn 
      const int newn  = static_cast<int> (num / dimRange);
      // local number of dof is local 
      const int local = (num % dimRange);
      // codim to be revised 
      return dimRange * (offset_[codim] + indexSet_.newIndex(newn,codim)) + local;
    }

    /** \copydoc Dune::DofMapper::numberOfHoles */
    int numberOfHoles ( const int codim ) const
    {
      return (maxDofs_[ codim ] > 0) ? 
        (dimRange * indexSet_.numberOfHoles( codim )) : 0;
    }

    /** \copydoc Dune::DofMapper::update */
    void update()
    {
      // assure that update is only called once per 
      // dof manager resize or compress 
      if( sequence_ != dm_.sequence() )
      {
        // calculate new size 
        size_ = 0;
        for( int codim = 0; codim <= dimension; ++codim )
        {
          oldOffSet_[ codim ] = offset_[ codim ];
          offset_[ codim ] = size_;
          size_ += indexSet_.size( codim ) * maxDofs_[ codim ];
        }
        sequence_ = dm_.sequence();
      }
    }

    /** \copydoc Dune::DofMapper::numBlocks
     */
    int numBlocks () const 
    {
      return dimension + 1;
    }

    /** \copydoc Dune::DofMapper::oldOffset
     */
    int oldOffSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return dimRange * oldOffSet_[ block ];
    }

    /** \copydoc Dune::DofMapper::newOffset
     */
    int offSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return dimRange * offset_[ block ];
    }

    /** \copydoc Dune::DofMapper::newSize
     */ 
    int newSize () const
    {
      int newSize = 0;
      for( int codim = 0; codim <= dimension; ++codim )
        newSize += indexSet_.size( codim ) * maxDofs_[ codim ];
      return dimRange * newSize;
    }

    /** \copydoc Dune::DofMapper::consecutive() const */
    bool consecutive () const
    {
      return BaseType :: checkConsecutive( indexSet_ );
    }
  };



  template< class GridPart, unsigned int polOrder, unsigned int dimR >
  struct LagrangeMapper< GridPart, polOrder, dimR > :: CodimCallInterface
  {
    typedef LagrangeMapper< GridPart, polOrder, dimR > MapperType;

    virtual ~CodimCallInterface ()
    {}

    virtual int subIndex ( const MapperType &mapper,
                           const EntityType &entity,
                           int i ) const = 0;

    virtual bool
    checkMapEntityDofToGlobal ( const MapperType &mapper,
                                const EntityType &entity,
                                const int subEntity,
                                const int localDof,
                                const int globalDof ) const = 0;
  };



  template< class GridPart, unsigned int polOrder, unsigned int dimR >
  template< unsigned int codim >
  struct LagrangeMapper< GridPart, polOrder, dimR > :: CodimCall
  : public CodimCallInterface
  {
    typedef LagrangeMapper< GridPart, polOrder, dimR > MapperType;

    typedef CodimCallInterface BaseType;

    virtual int subIndex ( const MapperType &mapper,
                           const EntityType &entity,
                           int i ) const
    {
      return mapper.indexSet_.template subIndex< codim >( entity, i );
    }

    virtual bool checkMapEntityDofToGlobal ( const MapperType &mapper,
                                             const EntityType &entity,
                                             const int subEntity,
                                             const int localDof,
                                             const int globalDof ) const
    {
      typedef Capabilities :: hasEntity< GridType, codim > HasEntity;

      return checkMapEntityDofToGlobal
        ( mapper, entity, subEntity, localDof, globalDof,
          MetaBool< HasEntity :: v >() );
    }

    inline bool
    checkMapEntityDofToGlobal ( const MapperType &mapper,
                                const EntityType &entity,
                                const int subEntity,
                                const int localDof,
                                const int globalDof,
                                const MetaBool< true > hasEntity ) const
    {
      typedef typename EntityType :: template Codim< codim > :: EntityPointer
        SubEntityPtrType;

      const SubEntityPtrType subEntityPtr
        = entity.template entity< codim >( subEntity );
      const int globalEntityDof
        = mapper.mapEntityDofToGlobal( *subEntityPtr, localDof );
      return (globalEntityDof == globalDof);
    }

    inline bool
    checkMapEntityDofToGlobal ( const MapperType &mapper,
                                const EntityType &entity,
                                const int subEntity,
                                const int localDof,
                                const int globalDof,
                                const MetaBool< false > hasEntity ) const
    {
      return true;
    }
  };
#endif

} // end namespace Dune 

#endif

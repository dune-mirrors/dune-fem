#ifndef DUNE_LAGRANGESPACE_MAPPER_HH
#define DUNE_LAGRANGESPACE_MAPPER_HH

//- Dune includes 
#include <dune/common/geometrytype.hh>
#include <dune/common/exceptions.hh>

//- Dune-Fem includes 
#include <dune/fem/misc/codimmap.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/dofmapperinterface.hh>

//- local includes 
#include "lagrangepoints.hh"

namespace Dune
{
  
  template< class GridPart, unsigned int polOrder, unsigned int dimrange >
  class LagrangeMapper;


  
  template< class GridPart, unsigned int polOrder, unsigned int dimrange >
  struct LagrangeMapperTraits
  {
    typedef GridPart GridPartType;
    
    enum { polynomialOrder = polOrder };

    //! dimension of the discrete function space's range
    enum { dimRange = dimrange };

    typedef typename GridPartType :: template Codim< 0 > :: IteratorType :: Entity
      EntityType;

    typedef LagrangeMapper< GridPartType, polynomialOrder, dimRange >
      DofMapperType;

    typedef DefaultDofMapIterator< EntityType, DofMapperType >
      DofMapIteratorType;
  };


  
  template< class GridPart, unsigned int dimrange >
  class LagrangeMapper< GridPart, 1, dimrange >
  : public DofMapperDefault< LagrangeMapperTraits< GridPart, 1, dimrange > >
  {
  public:
    typedef LagrangeMapperTraits< GridPart, 1, dimrange > Traits;
    
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
    enum { DimRange = Traits :: dimRange };

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

    /** \copydoc Dune::DofMapperInterface::size() const */
    int size () const
    {
      return dimRange * indexSet_.size( dimension );
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
    inline int mapToGlobal ( const EntityType &entity,
                             const int localDof ) const
    {
      const int coordinate = localDof % dimRange;
      const int localPoint = localDof / dimRange;
      const int globalPoint
        = indexSet_.template subIndex< dimension >( entity, localPoint );
      return dimRange * globalPoint + coordinate;
    }

    /** \copydoc Dune::DofMapperInterface::mapEntityDofToGlobal */
    template< class Entity >
    inline int mapEntityDofToGlobal ( const Entity &entity,
                                      const int localDof ) const
    {
      if( Entity :: codimension != (unsigned int) dimension )
        DUNE_THROW( RangeError, "No such local DoF." );

      assert( (localDof >= 0) && (localDof < dimRange) );
      return dimRange * indexSet_.index( entity ) + localDof;
    }
    
    /** \copydoc Dune::DofMapperInterface::maxNumDofs() const */
    int maxNumDofs () const
    {
      return dimRange * maxDofs_;
    }

    using BaseType :: numDofs;
    
    /** \copydoc Dune::DofMapperInterface::numDofs(const EntityType &entity) const */
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
   
    /** \copydoc Dune::DofMapperInterface::oldIndex */
    int oldIndex ( int hole, int ) const
    {
      const int coordinate = hole % dimRange;
      const int setHole = hole / dimRange;
      const int setIndex = indexSet_.oldIndex( setHole, dimension );
      return setIndex * dimRange + coordinate;
    }

    /** \copydoc Dune::DofMapperInterface::newIndex */
    int newIndex ( int hole , int ) const
    {
      const int coordinate = hole % dimRange;
      const int setHole = hole / dimRange;
      const int setIndex = indexSet_.newIndex( setHole, dimension );
      return setIndex * dimRange + coordinate;
    }

    /** \copydoc Dune::DofMapperInterface::numberOfHoles
     */
    int numberOfHoles ( int ) const
    {
      return dimRange * indexSet_.numberOfHoles( dimension );
    }

    /** \copydoc Dune::DofMapperInterface::newSize
     */
    int newSize () const
    {
      return this->size();
    }

    /** \copydoc Dune::DofMapperInterface::needsCompress
     */
    bool needsCompress () const
    {
      return indexSet_.needsCompress();
    }
  };



  template< class GridPart, unsigned int dimrange >
  class LagrangeMapper< GridPart, 2, dimrange >
  : public DofMapperDefault< LagrangeMapperTraits< GridPart, 2, dimrange > >
  {
  public:
    typedef LagrangeMapperTraits< GridPart, 2, dimrange > Traits;
    
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
    enum { DimRange = Traits :: dimRange };

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
    : dm_( DMFactoryType :: getDofManager(gridPart.grid()) )
    , indexSet_( gridPart.indexSet() )
    , lagrangePointSet_( lagrangePointSet )
    , sequence_( dm_.sequence() )
    {
      for( int codim = 0; codim <= dimension; ++codim )
        maxDofs_[ codim ] = 0;
      
      typedef typename LagrangePointSetMapType :: iterator IteratorType;
      IteratorType end = lagrangePointSet_.end();
      for( IteratorType it = lagrangePointSet_.begin(); it != end; ++it )
      {
        const LagrangePointSetType *set = (*it).second;
        if( set == NULL )
          continue;
        
        for( int codim = 0; codim <= dimension; ++codim )
        {
          const unsigned int setDofs = set->numDofs( codim );
          unsigned int &maxDofs = maxDofs_[ codim ];
          maxDofs = (maxDofs >= setDofs) ? maxDofs : setDofs;
        }
      }

      size_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
      {
        offset_[ codim ] = size_;
        oldOffSet_[ codim ] = size_;
        size_ += indexSet_.size( codim ) * maxDofs_[ codim ];
      }

      numDofs_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
        numDofs_ += maxDofs_[ codim ];
    }
    
    //! destructor 
    virtual ~LagrangeMapper ()
    {
    }

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
      
      const int subIndex
        = codimCall_[ dofInfo.codim ].subIndex( *this, entity, dofInfo.subEntity );
      
      const int globalDof
        = dimRange * (offset_[ dofInfo.codim ] + subIndex) + coordinate;
      assert( globalDof == codimCall_[ dofInfo.codim ].mapEntityDofToGlobal
                             ( *this, entity, dofInfo.subEntity, coordinate ) );
      return globalDof;
    }

    /** \copydoc Dune::DofMapperInterface::mapEntityDofToGlobal */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const 
    {
      assert( localDof < dimRange );
      const int globalDofPt
        = offset_[ Entity :: codimension ] + indexSet_.index( entity );
      return dimRange * globalDofPt + localDof;
    }
    
    /** \copydoc Dune::DofMapperInterface::maxNumDofs() const */
    inline int maxNumDofs () const
    {
      return numDofs_;
    }

    using BaseType :: numDofs;

    /** \copydoc Dune::DofMapperInterface::numDofs(const EntityType &entity) const */
    inline int numDofs ( const EntityType &entity ) const
    {
      return lagrangePointSet_[ entity.geometry().type() ]->size();
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
      {
        newSize += indexSet_.size( codim ) * maxDofs_[ codim ];
      }
      return dimRange * newSize;
    }

    /** \copydoc Dune::DofMapperInterface::needsCompress
     */
    bool needsCompress () const
    {
      return indexSet_.needsCompress();
    }
  };


  
  template< class GridPart, unsigned int dimrange >
  struct LagrangeMapper< GridPart, 2, dimrange > :: CodimCallInterface
  {
    typedef LagrangeMapper< GridPart, 2, dimrange > MapperType;

    virtual ~CodimCallInterface ()
    {}

    virtual int subIndex ( const MapperType &mapper,
                           const EntityType &entity,
                           int i ) const = 0;

    virtual int mapEntityDofToGlobal ( const MapperType &mapper,
                                       const EntityType &entity,
                                       int subEntity,
                                       int localDof ) const = 0;
  };



  template< class GridPart, unsigned int dimrange >
  template< unsigned int codim >
  struct LagrangeMapper< GridPart, 2, dimrange > :: CodimCall
  : public CodimCallInterface
  {
    typedef LagrangeMapper< GridPart, 2, dimrange > MapperType;

    typedef CodimCallInterface BaseType;

    virtual int subIndex ( const MapperType &mapper,
                           const EntityType &entity,
                           int i ) const
    {
      return mapper.indexSet_.template subIndex< codim >( entity, i );
    }

    virtual int mapEntityDofToGlobal ( const MapperType &mapper,
                                       const EntityType &entity,
                                       int subEntity,
                                       int localDof ) const
    {
      typedef typename EntityType :: template Codim< codim > :: EntityPointer
        SubEntityPtrType;

      const SubEntityPtrType subEntityPtr
        = entity.template entity< codim >( subEntity );
      return mapper.mapEntityDofToGlobal( *subEntityPtr, localDof );
    }
  };

 
} // end namespace Dune 

#endif

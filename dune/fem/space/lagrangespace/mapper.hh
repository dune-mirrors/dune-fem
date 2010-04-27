#ifndef DUNE_LAGRANGESPACE_MAPPER_HH
#define DUNE_LAGRANGESPACE_MAPPER_HH

//- Dune includes 
#include <dune/common/geometrytype.hh>
#include <dune/common/exceptions.hh>

//- Dune-Fem includes 
#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/misc/codimmap.hh>
#include <dune/fem/misc/metaprogramming.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/dofmapper.hh>

//- local includes 
#include "lagrangepoints.hh"

namespace Dune
{
  
  template< class GridPart, int polOrder, int dimR >
  class LagrangeMapper;


  
  template< class GridPart, int polOrder, int dimR >
  struct LagrangeMapperTraits
  {
    typedef GridPart GridPartType;
    
    static const int polynomialOrder = polOrder;

    //! dimension of the discrete function space's range
    static const int dimRange = dimR;

    typedef typename GridPartType::template Codim< 0 >::IteratorType::Entity EntityType;
    typedef LagrangeMapper< GridPartType, polynomialOrder, dimRange > DofMapperType;
    typedef DefaultDofMapIterator< EntityType, DofMapperType > DofMapIteratorType;
  };


  
  template< class GridPart, int dimR >
  class LagrangeMapper< GridPart, 1, dimR >
  : public DofMapperDefault< LagrangeMapperTraits< GridPart, 1, dimR > >
  {
    typedef LagrangeMapper< GridPart, 1, dimR > ThisType;
    typedef DofMapperDefault< LagrangeMapperTraits< GridPart, 1, dimR > > BaseType;

  public:
    typedef LagrangeMapperTraits< GridPart, 1, dimR > Traits;
    
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

    //! dimension of the discrete function space's range
    static const int dimRange = Traits::dimRange;

    //! type of the Lagrange point set
    typedef LagrangePointSet< GridPartType, polynomialOrder >
      LagrangePointSetType;
    //! type of the map for the Lagrange point sets
    typedef std::map< const GeometryType, const LagrangePointSetType* >
      LagrangePointSetMapType;

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
    {}

    /** \copydoc Dune::DofMapper::size() const */
    int size () const
    {
      return dimRange * indexSet_.size( dimension );
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
      const int coordinate = localDof % dimRange;
      const int localPoint = localDof / dimRange;
      const int globalPoint = indexSet_.subIndex( entity, localPoint, dimension );
      return dimRange * globalPoint + coordinate;
    }

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const
    {
      if( Entity::codimension != dimension )
        DUNE_THROW( RangeError, "No such local DoF." );

      assert( (localDof >= 0) && (localDof < dimRange) );
      return dimRange * indexSet_.index( entity ) + localDof;
    }
    
    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    int maxNumDofs () const
    {
      return dimRange * maxDofs_;
    }

    using BaseType::numDofs;
    
    /** \copydoc Dune::DofMapper::numDofs(const EntityType &entity) const */
    int numDofs ( const EntityType &entity ) const
    {
      return dimRange * entity.template count< dimension >();
    }

    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      return (Entity::codimension == dimension ? dimRange : 0);
    }

    /** \brief Check, whether any DoFs are associated with a codimension */
    bool contains ( int codim ) const
    {
      return (codim == dimension);
    }

    /** \brief Check, whether the data in a codimension has fixed size */
    bool fixedDataSize ( int codim ) const
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

    /** \copydoc Dune::DofMapper::consecutive() const */
    bool consecutive () const
    {
      return BaseType::checkConsecutive( indexSet_ );
    }

  private:
    const IndexSetType &indexSet_;
    unsigned int maxDofs_;
  };



  template< class GridPart, int dimR >
  class LagrangeMapper< GridPart, 2, dimR >
  : public DofMapperDefault< LagrangeMapperTraits< GridPart, 2, dimR > >
  {
    typedef LagrangeMapper< GridPart, 2, dimR > ThisType;
    typedef DofMapperDefault< LagrangeMapperTraits< GridPart, 2, dimR > > BaseType;

  public:
    typedef LagrangeMapperTraits< GridPart, 2, dimR > Traits;
    
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

    //! dimension of the discrete function space's range
    static const int dimRange = Traits::dimRange;

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

  public:
    //! constructor
    LagrangeMapper ( const GridPartType &gridPart,
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
        
        numDofs_ = std :: max( numDofs_, set->size() );
        for( int codim = 0; codim <= dimension; ++codim )
          maxDofs_[ codim ] = std::max( maxDofs_[ codim ], set->maxDofs( codim ) );
      }

      computeOffsets();
      dm_.addIndexSet( *this );
    }
    
    //! destructor 
    virtual ~LagrangeMapper ()
    {
      dm_.removeIndexSet( *this );
    }

    //! return overall number of degrees of freedom 
    int size () const
    {
      return dimRange * size_;
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
    int mapToGlobal ( const EntityType &entity, const int local ) const
    {
      const int coordinate = local % dimRange;
      const int localDof = local / dimRange;
      
      // unsigned int codim, subEntity;
      const LagrangePointSetType *set = lagrangePointSet_[ entity.type() ];
      const DofInfo &dofInfo = set->dofInfo( localDof );

      const unsigned int codim = dofInfo.codim;
      const int subIndex = indexSet_.subIndex( entity, dofInfo.subEntity, codim );

      const int globalDof = dimRange * (offset_[ codim ] + subIndex) + coordinate;

      return globalDof;
    }

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const 
    {
      static const unsigned int codim = Entity::codimension;
      assert( localDof < numEntityDofs( entity ) );
      const int globalDofPt = offset_[ codim ] + indexSet_.index( entity );
      return dimRange * globalDofPt + localDof;
    }
    
    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    int maxNumDofs () const
    {
      return dimRange * numDofs_;
    }

    using BaseType::numDofs;

    /** \copydoc Dune::DofMapper::numDofs(const EntityType &entity) const */
    int numDofs ( const EntityType &entity ) const
    {
      return dimRange * lagrangePointSet_[ entity.type() ]->size();
    }

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      // This implementation only works for nonhybrid grids (simplices or cubes)
      return dimRange * maxDofs_[ Entity::codimension ];
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
    int oldIndex ( const int num, const int codim ) const
    {
      // corresponding number of set is newn 
      const int newn  = static_cast<int> (num / dimRange);
      // local number of dof is local 
      const int local = (num % dimRange);
      // codim to be revised 
      return dimRange * (offset_[codim] + indexSet_.oldIndex(newn,codim)) + local;
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
    LagrangeMapper ( const ThisType & );
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

    // reference to dof manager
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



  // Higher Order Lagrange Mapper
  // ----------------------------
  //
  // Note: This mapper assumes that the grid is "twist-free".

#ifdef USE_TWISTFREE_MAPPER
  template< class GridPart, int polOrder, int dimR >
  class LagrangeMapper
  : public DofMapperDefault< LagrangeMapperTraits< GridPart, polOrder, dimR > >
  {
    typedef LagrangeMapper< GridPart, polOrder, dimR > ThisType;
    typedef DofMapperDefault< LagrangeMapperTraits< GridPart, polOrder, dimR > > BaseType;

  public:
    typedef LagrangeMapperTraits< GridPart, polOrder, dimR > Traits;
    
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

    //! dimension of the discrete function space's range
    static const int dimRange = Traits::dimRange;

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
    LagrangeMapper ( const GridPartType &gridPart,
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
    virtual ~LagrangeMapper ()
    {
      dm_.removeIndexSet( *this );
    }

    //! return overall number of degrees of freedom 
    int size () const
    {
      return dimRange * size_;
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
      const int coordinate = localDof % dimRange;
      const int scalarDof = localDof / dimRange;
      
      // unsigned int codim, subEntity;
      const LagrangePointSetType *set = lagrangePointSet_[ entity.type() ];
      const DofInfo& dofInfo = set->dofInfo( scalarDof );
      
      const int entityDof = dofInfo.dofNumber * dimRange + coordinate;

      const unsigned int codim = dofInfo.codim;
      const int subIndex = indexSet_.subIndex( entity, dofInfo.subEntity, codim );

      const int globalDof = dimRange * (offset_[ codim ] + subIndex * maxDofs_[ codim ]) + entityDof;

      return globalDof;
    }

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const 
    {
      const unsigned int codim = Entity::codimension;

      assert( localDof < numEntityDofs( entity ) );
      const int offset = offset_[ codim ] + indexSet_.index( entity ) * maxDofs_[ codim ];
      return dimRange * offset + localDof;
    }
    
    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    int maxNumDofs () const
    {
      return dimRange * numDofs_;
    }

    using BaseType::numDofs;

    /** \copydoc Dune::DofMapper::numDofs(const EntityType &entity) const */
    int numDofs ( const EntityType &entity ) const
    {
      return dimRange * lagrangePointSet_[ entity.type() ]->size();
    }

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      // This implementation only works for nonhybrid grids (simplices or cubes)
      return dimRange * maxDofs_[ Entity::codimension ];
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
    int oldIndex ( const int num, const int codim ) const
    {
      // corresponding number of set is newn 
      const int newn  = static_cast<int> (num / dimRange);
      // local number of dof is local 
      const int local = (num % dimRange);
      // codim to be revised 
      return dimRange * (offset_[codim] + indexSet_.oldIndex(newn,codim)) + local;
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
      return dimRange * oldOffSet_[ block ];
    }

    /** \copydoc Dune::DofMapper::offSet(const int block) const */
    int offSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return dimRange * offset_[ block ];
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
    LagrangeMapper ( const ThisType & );
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

#ifndef DUNE_FEM_CODIMENSIONMAPPER_HH
#define DUNE_FEM_CODIMENSIONMAPPER_HH

#if HAVE_DUNE_GEOMETRY
#include <dune/geometry/referenceelements.hh>
#else
#include <dune/grid/common/genericreferenceelements.hh>
#endif

#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/mapper/dofmapper.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class GridPart, int codim >
  class CodimensionMapper;

  template< class Element, class DofMapper >
  class CodimensionDofMapIterator;



  // CodimensionMapperTraits
  // -----------------------

  template< class GridPart, int codim >
  struct CodimensionMapperTraits
  {
    typedef GridPart GridPartType;

    // we still need entities of codimension 0 here
    typedef typename GridPartType::template Codim< 0 >::EntityType ElementType;

    typedef typename GridPartType::IndexSetType IndexSetType;
    
    typedef CodimensionMapper< GridPartType, codim > DofMapperType;

    typedef DefaultDofMapIterator< ElementType, DofMapperType > DofMapIteratorType;
  };


  template< class GridPart >
  struct CodimensionMapperTraits< GridPart, 0 >
  {
    typedef GridPart GridPartType;

    typedef typename GridPartType::template Codim< 0 >::EntityType ElementType;

    typedef typename GridPartType::IndexSetType IndexSetType;
    
    typedef CodimensionMapper< GridPartType, 0 > DofMapperType;

    typedef CodimensionDofMapIterator< ElementType, DofMapperType > DofMapIteratorType;
  };



  /** \class CodimensionMapper
   *  \brief mapper allocating one DoF per subentity of a given codimension
   *
   *  \tparam  GridPart  grid part, the mapper shall be used on
   *  \tparam  cdim      codimension
   */
  template< class GridPart, int cdim >
  class CodimensionMapper
  : public DofMapperDefault< CodimensionMapperTraits< GridPart, cdim > >
  {
    typedef CodimensionMapper< GridPart, cdim > ThisType;
    typedef DofMapperDefault< CodimensionMapperTraits< GridPart, cdim > > BaseType;

  public:
    typedef typename BaseType::Traits Traits;

    typedef typename Traits::ElementType ElementType;

    typedef typename Traits::GridPartType GridPartType;

    typedef typename Traits::IndexSetType IndexSetType;

    typedef typename Traits::DofMapIteratorType DofMapIteratorType;

    //! dimension of the grid
    static const int dimension = GridPartType::dimension;

    //! codimension that is mapped 
    static const int codimension = cdim;

    //! Constructor
    explicit CodimensionMapper( const GridPartType &gridPart );

    /** \copydoc DofMapper::contains */
    bool contains ( const int codim ) const 
    {
      return (codim == codimension);
    }

    /** \copydoc DofMapper::size */
    int size () const
    {
      // return number of dofs for codimension 
      return indexSet_.size( codimension );
    }

    /** \copydoc Dune::DofMapper::begin(const ElementType &entity) const */
    DofMapIteratorType begin ( const ElementType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::beginIterator, entity, *this );
    }
    
    /** \copydoc Dune::DofMapper::end(const ElementType &entity) const */
    DofMapIteratorType end ( const ElementType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::endIterator, entity, *this );
    }

    /** \copydoc DofMapper::mapToGlobal */
    int mapToGlobal ( const ElementType &entity, const int localDof ) const
    {
      // we only have one local dof 
      assert( localDof < maxNumDofs() );
      if ( codimension == 0 ) 
        return indexSet_.index( entity );
      else 
        return indexSet_.subIndex( entity, localDof, codimension );
    }

  protected:
    template <int codim1, int codim>
    struct IndexExtractor 
    {
      template <class Entity >
      static inline int index(const IndexSetType& indexSet, const Entity& entity )
      {
        DUNE_THROW(InvalidStateException,"Wrong codimension selected"); 
        return -1;
      }
    };

    //! return index of entity (specialized because of Grids like Yasp which are not
    //! supporting every codimension entity)
    template <int codim>
    struct IndexExtractor< codim, codim >
    {
      template< class Entity >
      static inline int index(const IndexSetType& indexSet, const Entity& entity )
      {
        return indexSet.index( entity );
      }
    };

  public:  
    /** \copydoc DofMapper::mapEntityDofToGlobal */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const
    {
      assert( codimension == Entity::codimension );
      return IndexExtractor< codimension,  Entity::codimension > :: index( indexSet_, entity ); 
    }

    /** \copydoc DofMapper::maxNumDofs() const */
    int maxNumDofs () const
    {
      return maxNumberOfDofs_;
    }

    /** \copydoc Dune::DofMapper::numDofs(const ElementType &entity) const */
    int numDofs ( const ElementType &entity ) const
    {
      return entity.template count< codimension >();
    }

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      return (contains( Entity::codimension ) ? 1 : 0);
    }
   
    /** \copydoc DofMapper::oldIndex */
    int oldIndex ( const int hole, int ) const
    {
      // forward to index set 
      return indexSet_.oldIndex( hole, codimension ) ;
    }

    /** \copydoc DofMapper::newIndex */
    int newIndex ( const int hole, int ) const
    {
      // forward to index set 
      return indexSet_.newIndex( hole, codimension );
    }

    /** \copydoc DofMapper::numberOfHoles */
    int numberOfHoles ( const int ) const
    {
      return indexSet_.numberOfHoles( codimension );
    }

    /** \copydoc DofMapper::consecutive */
    bool consecutive () const 
    {
      return BaseType::checkConsecutive( indexSet_ );
    }

  protected:
    // index set for the grid
    const IndexSetType &indexSet_;

    // maximal number of local dofs
    int maxNumberOfDofs_;
  };


  template< class GridPart, int cdim >
  inline CodimensionMapper< GridPart, cdim >
    ::CodimensionMapper( const GridPartType &gridPart )
  : indexSet_( gridPart.indexSet() ),
    maxNumberOfDofs_( 0 )
  {
    typedef typename GridPartType::ctype ctype;
    typedef GenericReferenceElements< ctype, dimension > RefElements;

    AllGeomTypes< IndexSetType, typename GridPartType::GridType > allTypes( indexSet_ );
    const std::vector< GeometryType > &types = allTypes.geomTypes( 0 );
    const unsigned int numTypes = types.size();
    for( unsigned int i = 0; i < numTypes; ++i )
    {
      const GeometryType &type = types[ i ];
      
      const int numSubEntities = RefElements::general( type ).size( codimension );
      maxNumberOfDofs_ = std::max( maxNumberOfDofs_, numSubEntities );
    }
    assert( maxNumberOfDofs_ > 0 );
  }



  // CodimensionDofMapIterator
  // -------------------------

  template< class Element, class DofMapper >
  class CodimensionDofMapIterator
  {
    typedef CodimensionDofMapIterator< Element, DofMapper > ThisType;

  public:
    typedef Element ElementType;
    typedef DofMapper DofMapperType;

    enum IteratorType { beginIterator, endIterator };

    CodimensionDofMapIterator ( const IteratorType type,
                                const ElementType &entity,
                                const DofMapperType &mapper )
    : baseIndex_( (type == endIterator) ? -1 : mapper.mapToGlobal( entity, 0 ) ),
      dof_( (type == endIterator) ? mapper.numDofs( entity ) : 0 )
    {}

    CodimensionDofMapIterator ( const ThisType &other )
    : baseIndex_( other.baseIndex_ ),
      dof_( other.dof_ )
    {}

    ThisType &operator++ ()
    {
      ++dof_;
      return *this;
    }

    bool operator== ( const ThisType &other ) const
    {
      return dof_ == other.dof_;
    }

    bool operator!= ( const ThisType &other ) const
    {
      return dof_ != other.dof_;
    }

    int local () const
    {
      return dof_;
    }

    int global () const
    {
      return baseIndex_ + dof_;
    }

  protected:
    const int baseIndex_;
    int dof_;
  };



  // CodimensionMapperSingletonFactory
  // ---------------------------------

  template< class GridPart, int cdim >
  struct CodimensionMapperSingletonFactory
  {
    typedef CodimensionMapper< GridPart, cdim > Object;

    struct Key
    {
      Key ( const GridPart &gp )
      : gridPart( gp )
      {}

      bool operator== ( const Key &other )
      {
        return (&gridPart.indexSet() == &other.gridPart.indexSet() );
      }

      const GridPart &gridPart;
    };

    //! create new mapper  
    static Object *createObject ( const Key &key )
    {
      return new Object ( key.gridPart );
    }

    //! delete mapper object 
    static void deleteObject ( Object *object )
    {
      delete object;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_CODIMENSIONMAPPER_HH

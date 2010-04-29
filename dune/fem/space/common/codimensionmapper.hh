#ifndef DUNE_CODIMENSIONMAPPER_HH
#define DUNE_CODIMENSIONMAPPER_HH

#include <dune/fem/space/common/dofmapper.hh>

namespace Dune
{

  template< class GridPartImp, int polOrder, int codim >
  class CodimensionMapper;

  class CodimensionDofMapIterator;



  template< class GridPartImp, int polOrder, int codim >
  struct CodimensionMapperTraits
  {
    typedef GridPartImp GridPartType;

    typedef typename GridPartType :: template Codim< codim > :: IteratorType :: Entity
      EntityType;

    typedef typename GridPartType :: IndexSetType IndexSetType;
    
    typedef CodimensionDofMapIterator DofMapIteratorType;

    typedef CodimensionMapper< GridPartType, polOrder, codim >
      DofMapperType;
  };



  //***************************************************************************
  //
  //!  The CodimensionMapper for mapping of local dof numbers to global dof numbers, 
  //!  i.e. the entry in the vector of unknowns
  //
  //***************************************************************************
  template< class GridPartImp, int polOrder, int codim >
  class CodimensionMapper
  : public DofMapperDefault< CodimensionMapperTraits< GridPartImp, polOrder, codimension> >
  {
  public:
    typedef CodimensionMapperTraits< GridPartImp, polOrder, codimension> Traits;

    typedef typename Traits :: EntityType EntityType;

    typedef typename Traits :: GridPartType GridPartType;

    typedef typename Traits :: IndexSetType IndexSetType;

    typedef typename Traits :: DofMapIteratorType DofMapIteratorType;

    //! codimension that is mapped 
    enum { codimension = codim };

  private:
    typedef CodimensionMapper< GridPartType, polOrder, codimension > ThisType;
    typedef DofMapperDefault< Traits > BaseType;

  protected:
    // index set of grid, i.e. the element indices 
    const IndexSetType &indexSet_;

    // number of dofs on element 
    const int numberOfDofs_;

  public:
    //! Constructor 
    inline CodimensionMapper( const GridPartType &gridPart,
                              const int numDofs )
    : indexSet_( gridPart.indexSet() ),
      numberOfDofs_( numDofs )
    {}

    //! return size of function space 
    //! see dofmanager.hh for definition of IndexSet, which 
    //! is a wrapper for en.index 
    /** \copydoc DofMapper::size */
    int size () const
    {
      // return number of dofs * number of elements 
      return (numberOfDofs_ * indexSet_.size( codimension ));
    }

    /** \copydoc Dune::DofMapper::begin(const EntityType &entity) const */
    inline DofMapIteratorType begin ( const EntityType &entity ) const
    {
      const int baseIndex = indexSet_.index( entity ) * numberOfDofs_;
      typename DofMapIteratorType :: BeginIterator type;
      return DofMapIteratorType( type, baseIndex );
    }
    
    /** \copydoc Dune::DofMapper::end(const EntityType &entity) const */
    inline DofMapIteratorType end ( const EntityType &entity ) const
    {
      typename DofMapIteratorType :: EndIterator type;
      return DofMapIteratorType( type, numberOfDofs_ );
    }

    /** \copydoc DofMapper::mapToGlobal */
    int mapToGlobal ( const EntityType &entity, const int localDof ) const
    {
      const int baseIndex = indexSet_.index( entity ) * numberOfDofs_;
      return baseIndex + localDof;
    }

    /** \copydoc DofMapper::maxNumDofs() const */
    int maxNumDofs () const
    {
      return numberOfDofs_;
    }
   
    /** \copydoc DofMapper::oldIndex */
    int oldIndex (const int hole, int ) const
    {
      // corresponding number of set is newn 
      const int newn  = static_cast<int> (hole / numberOfDofs_);
      // local number of dof is local 
      const int local = (hole % numberOfDofs_);
      return (numberOfDofs_ * indexSet_.oldIndex( newn, codimension )) + local;
    }

    /** \copydoc DofMapper::newIndex */
    int newIndex (const int hole, int ) const
    {
      // corresponding number of set is newn 
      const int newn = static_cast<int> (hole / numberOfDofs_);
      // local number of dof is local 
      const int local = (hole % numberOfDofs_);
      return (numberOfDofs_ * indexSet_.newIndex( newn, codimension )) + local;
    }

    /** \copydoc DofMapper::numberOfHoles */
    int numberOfHoles ( const int ) const
    {
      // this index set works only for codim = 0 at the moment
      return numberOfDofs_ * indexSet_.numberOfHoles( codimension );
    }

    /** \copydoc DofMapper::consecutive */
    bool consecutive() const 
    {
      return BaseType::checkConsecutive(indexSet_);
    }
  };



  class CodimensionDofMapIterator
  {
  public:
    struct BeginIterator {};
    struct EndIterator {};

  private:
    typedef CodimensionDofMapIterator ThisType;
    
  protected:
    const int baseIndex_;
    int dof_;
    
  public:
    inline CodimensionDofMapIterator ( const BeginIterator type,
                              const int baseIndex )
    : baseIndex_( baseIndex ),
      dof_( 0 )
    {}

    inline CodimensionDofMapIterator ( const EndIterator type,
                              const int numDofs )
    : baseIndex_( -1 ),
      dof_( numDofs )
    {}

    inline CodimensionDofMapIterator ( const ThisType &other )
    : baseIndex_( other.baseIndex_ ),
      dof_( other.dof_ )
    {}

    inline ThisType &operator++ ()
    {
      ++dof_;
      return *this;
    }

    inline bool operator== ( const ThisType &other ) const
    {
      return dof_ == other.dof_;
    }

    inline bool operator!= ( const ThisType &other ) const
    {
      return dof_ != other.dof_;
    }

    inline int local () const
    {
      return dof_;
    }

    inline int global () const
    {
      return baseIndex_ + dof_;
    }
  };


} // end namespace Dune

#endif

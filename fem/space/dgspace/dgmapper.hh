#ifndef DUNE_DGMAPPER_HH
#define DUNE_DGMAPPER_HH

#include <dune/fem/space/common/dofmapper.hh>

namespace Dune
{

  template< class GridPartImp, int polOrder, int dimRange >
  class DGMapper;

  class DGDofMapIterator;



  template< class GridPartImp, int polOrder, int dimRange >
  struct DGMapperTraits
  {
    typedef GridPartImp GridPartType;

    typedef typename GridPartType :: template Codim< 0 > :: IteratorType :: Entity
      EntityType;

    typedef typename GridPartType :: IndexSetType IndexSetType;
    
    typedef DGDofMapIterator DofMapIteratorType;

    typedef DGMapper< GridPartType, polOrder, dimRange >
      DofMapperType;
  };



  //***************************************************************************
  //
  //!  DG Mapper for mapping of local dof numbers to global dof numbers, 
  //!  i.e. the entry in the vector of unknowns
  //
  //***************************************************************************
  template< class GridPartImp, int polOrder, int dimRange >
  class DGMapper
  : public DofMapperDefault< DGMapperTraits< GridPartImp, polOrder, dimRange > >
  {
  public:
    typedef DGMapperTraits< GridPartImp, polOrder, dimRange > Traits;

    typedef typename Traits :: EntityType EntityType;

    typedef typename Traits :: GridPartType GridPartType;

    typedef typename Traits :: IndexSetType IndexSetType;

    typedef typename Traits :: DofMapIteratorType DofMapIteratorType;

  private:
    typedef DGMapper< GridPartType, polOrder, dimRange > ThisType;
    typedef DofMapperDefault< Traits > BaseType;

  protected:
    // index set of grid, i.e. the element indices 
    const IndexSetType &indexSet_;

    // number of dofs on element 
    const int numberOfDofs_;

  public:
    //! Constructor 
    inline DGMapper( const GridPartType &gridPart,
                     int numDofs )
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
      return (numberOfDofs_ * indexSet_.size( 0 ));
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
   
    /** \copydoc DofMapper::newSize */
    int newSize() const
    {
      return this->size();
    }

    /** \copydoc DofMapper::oldIndex */
    int oldIndex (const int hole, int ) const
    {
      // corresponding number of set is newn 
      const int newn  = static_cast<int> (hole / numberOfDofs_);
      // local number of dof is local 
      const int local = (hole % numberOfDofs_);
      return (numberOfDofs_ * indexSet_.oldIndex(newn,0)) + local;
    }

    /** \copydoc DofMapper::newIndex */
    int newIndex (const int hole, int ) const
    {
      // corresponding number of set is newn 
      const int newn = static_cast<int> (hole / numberOfDofs_);
      // local number of dof is local 
      const int local = (hole % numberOfDofs_);
      return (numberOfDofs_ * indexSet_.newIndex(newn,0)) + local;
    }

    /** \copydoc DofMapper::numberOfHoles */
    int numberOfHoles ( int ) const
    {
      // this index set works only for codim = 0 at the moment
      return numberOfDofs_ * indexSet_.numberOfHoles(0);
    }

    /** \copydoc DofMapper::consecutive */
    bool consecutive() const 
    {
      return BaseType::checkConsecutive(indexSet_);
    }
  };



  class DGDofMapIterator
  {
  public:
    struct BeginIterator {};
    struct EndIterator {};

  private:
    typedef DGDofMapIterator ThisType;
    
  protected:
    const int baseIndex_;
    int dof_;
    
  public:
    inline DGDofMapIterator ( const BeginIterator type,
                              const int baseIndex )
    : baseIndex_( baseIndex ),
      dof_( 0 )
    {}

    inline DGDofMapIterator ( const EndIterator type,
                              const int numDofs )
    : baseIndex_( -1 ),
      dof_( numDofs )
    {}

    inline DGDofMapIterator ( const ThisType &other )
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

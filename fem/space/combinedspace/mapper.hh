#ifndef DUNE_FEM_COMBINEDSPACE_MAPPER_HH
#define DUNE_FEM_COMBINEDSPACE_MAPPER_HH

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include <dune/fem/space/common/dofmapper.hh>
#include <dune/fem/space/combinedspace/combineddofstorage.hh>

namespace Dune
{
  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  class CombinedMapper;

  template< class CombinedMapperTraits >
  class CombinedDofMapIterator;



  template< class ContainedMapper, int N, DofStoragePolicy policy >
  struct CombinedMapperTraits
  {
    typedef ContainedMapper ContainedMapperType;

    enum { numComponents = N };

    typedef CombinedMapperTraits< ContainedMapperType, numComponents, policy >
      Traits;

    typedef typename ContainedMapperType :: EntityType EntityType;

    typedef CombinedDofConversionUtility< ContainedMapperType, policy >
      GlobalDofConversionUtilityType;
    
    typedef CombinedDofMapIterator< Traits > DofMapIteratorType;

    typedef CombinedMapper< ContainedMapperType, numComponents, policy >
      DofMapperType;
  };

  

  template< class CombinedMapperTraits >
  class CombinedDofMapIterator
  {
  public:
    typedef CombinedMapperTraits Traits;

  private:
    typedef CombinedDofMapIterator< Traits > ThisType;

  public:
    typedef typename Traits :: ContainedMapperType ContainedMapperType;

    typedef typename ContainedMapperType :: DofMapIteratorType
      ContainedDofMapIteratorType;

    typedef typename Traits :: GlobalDofConversionUtilityType
      GlobalDofConversionUtilityType;

    enum { numComponents = Traits :: numComponents };

  protected:
    const GlobalDofConversionUtilityType dofUtil_;

    ContainedDofMapIteratorType containedIterator_;
    int component_;
    
  public:
    inline
    CombinedDofMapIterator ( const ContainedDofMapIteratorType &containedIterator,
                             const GlobalDofConversionUtilityType &dofUtil )
    : dofUtil_( dofUtil ),
      containedIterator_( containedIterator ),
      component_( 0 )
    {}

    inline CombinedDofMapIterator ( const ThisType &other )
    : dofUtil_( other.dofUtil_ ),
      containedIterator_( other.containedIterator_ ),
      component_( other.component_ )
    {}

    inline ThisType &operator++ ()
    {
      ++component_;
      if( component_ >= numComponents )
      {
        ++containedIterator_;
        component_ = 0;
      }
      return *this;
    }

    inline bool operator== ( const ThisType &other ) const
    {
      return (containedIterator_ == other.containedIterator_)
             && (component_ == other.component_);
    }

    inline bool operator!= ( const ThisType &other ) const
    {
      return !(*this == other);
    }

    inline int local () const
    {
      return containedIterator_.local() * numComponents + component_;
    }

    inline int global () const
    {
      return dofUtil_.combinedDof( containedIterator_.global(), component_ );
    }

    inline int component () const
    {
      return component_;
    }

    inline int localScalar () const
    {
      return containedIterator_.local();
    }

    inline int globalScalar () const
    {
      return containedIterator_.global();
    }
  };



  /** \class CombinedMapper
   *  \ingroup CombinedSpace
   *  \brief DofMapper for the CombinedSpace
   */
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  class CombinedMapper
  : public DofMapperDefault< CombinedMapperTraits< ContainedMapper, N, policy > >
  {
  public:
    typedef CombinedMapperTraits< ContainedMapper, N, policy > Traits;

    enum { numComponents = Traits :: numComponents };

    typedef typename Traits :: ContainedMapperType ContainedMapperType;

  private:
    typedef CombinedMapper< ContainedMapperType, numComponents, policy >
      ThisType;
    typedef DofMapperDefault< Traits > BaseType;

    template< class ContainedSpace, int NC, DofStoragePolicy p >
    friend class CombinedSpace;

  public:
    typedef typename Traits :: EntityType EntityType;

    typedef typename Traits :: DofMapIteratorType DofMapIteratorType;

    typedef typename Traits :: GlobalDofConversionUtilityType
      GlobalDofConversionUtilityType;
    
    typedef CombinedDofConversionUtility< ContainedMapperType, PointBased >
      LocalDofConversionUtilityType;

  protected:
    //- Data members
    ContainedMapperType *containedMapper_;

    const LocalDofConversionUtilityType utilLocal_;
    GlobalDofConversionUtilityType utilGlobal_;
    //int oldSize_, size_;

  public:
    //! Constructor
    inline explicit CombinedMapper ( ContainedMapperType &mapper );

  private:
    // prohibit copying
    inline CombinedMapper ( const ThisType & );

  public:
    //! Total number of degrees of freedom
    inline int size () const;
    
    /** \copydoc Dune::DofMapper::begin(const EntityType &entity) const */
    inline DofMapIteratorType begin ( const EntityType &entity ) const;
    
    /** \copydoc Dune::DofMapper::end(const EntityType &entity) const */
    inline DofMapIteratorType end ( const EntityType &entity ) const;

    /** \copydoc Dune::DofMapper::mapToGlobal(const EntityType &entity,const int localDof) const */
    inline int mapToGlobal( const EntityType &entity,
                            const int localDof ) const;

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal(const Entity &entity,const int localDof) const */
    template< class Entity >
    inline int mapEntityDofToGlobal( const Entity &entity,
                                     const int localDof ) const;

    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    inline int maxNumDofs () const;

    /** \copydoc Dune::DofMapper::numDofs(const EntityType &entity) const */
    inline int numDofs ( const EntityType &entity ) const;

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    inline int numEntityDofs ( const Entity &entity ) const;

    //- Method inherited from mapper interface
    //! if grid has changed determine new size 
    //! (to be called once per timestep, therefore virtual )
    inline int newSize() const;
  
    //! return new index in dof array 
    inline int newIndex ( const int hole, const int block ) const;

    //! return old index in dof array of given index ( for dof compress ) 
    inline int oldIndex ( const int hole, const int block ) const;

    //! return number of holes in the data 
    inline int numberOfHoles ( const int block ) const;
  
    //! returnn number of mem blocks 
    inline int numBlocks () const; 

    //! update offset information
    inline void update (); 
    
    //! return current old offset of block 
    inline int oldOffSet ( const int block ) const;

    //! return current offset of block 
    inline int offSet ( const int block ) const;

    //! return true if compress will affect data  
    inline bool consecutive () const;

  protected:
    inline ContainedMapperType &containedMapper () const;

    inline static int chooseSize ( int pointBased,
                                   int variableBased,
                                   Int2Type<PointBased> );

    inline static int chooseSize ( int pointBased,
                                   int variableBased,
                                   Int2Type<VariableBased> );

  }; // end class CombinedMapper
  
  /** @} **/  
  
} // end namespace Dune

#include "mapper_inline.hh"

#endif

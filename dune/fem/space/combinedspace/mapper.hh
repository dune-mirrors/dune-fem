#ifndef DUNE_FEM_COMBINEDSPACE_MAPPER_HH
#define DUNE_FEM_COMBINEDSPACE_MAPPER_HH

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include <dune/fem/space/mapper/dofmapper.hh>
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

    typedef typename ContainedMapperType::ElementType ElementType;

    typedef CombinedDofConversionUtility< ContainedMapperType, numComponents, policy >
      GlobalDofConversionUtilityType;
    
    typedef CombinedDofMapIterator< Traits > DofMapIteratorType;

    typedef CombinedMapper< ContainedMapperType, numComponents, policy >
      DofMapperType;
  };

  

  template< class CombinedMapperTraits >
  class CombinedDofMapIterator
  {
    typedef CombinedDofMapIterator< CombinedMapperTraits > ThisType;

  public:
    typedef CombinedMapperTraits Traits;

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

    CombinedDofMapIterator ( const ThisType &other )
    : dofUtil_( other.dofUtil_ ),
      containedIterator_( other.containedIterator_ ),
      component_( other.component_ )
    {}

    ThisType &operator++ ()
    {
      ++component_;
      if( component_ >= numComponents )
      {
        ++containedIterator_;
        component_ = 0;
      }
      return *this;
    }

    bool operator== ( const ThisType &other ) const
    {
      return (containedIterator_ == other.containedIterator_)
             && (component_ == other.component_);
    }

    bool operator!= ( const ThisType &other ) const
    {
      return !(*this == other);
    }

    int local () const
    {
      return containedIterator_.local() * numComponents + component_;
    }

    int global () const
    {
      return dofUtil_.combinedDof( containedIterator_.global(), component_ );
    }

    int component () const
    {
      return component_;
    }

    int localScalar () const
    {
      return containedIterator_.local();
    }

    int globalScalar () const
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
    typedef CombinedMapper< ContainedMapper, N, policy > ThisType;
    typedef DofMapperDefault< CombinedMapperTraits< ContainedMapper, N, policy > > BaseType;

    template< class, int, DofStoragePolicy >
    friend class CombinedSpace;

    template< class Functor >
    struct FunctorWrapper;

  public:
    typedef CombinedMapperTraits< ContainedMapper, N, policy > Traits;

    enum { numComponents = Traits :: numComponents };

    typedef typename Traits :: ContainedMapperType ContainedMapperType;

    typedef typename Traits::ElementType ElementType;

    typedef typename Traits :: DofMapIteratorType DofMapIteratorType;

    typedef typename Traits :: GlobalDofConversionUtilityType
      GlobalDofConversionUtilityType;
    
    typedef CombinedDofConversionUtility< ContainedMapperType, numComponents, PointBased >
      LocalDofConversionUtilityType;

  protected:
    //- Data members
    ContainedMapperType *containedMapper_;

    const LocalDofConversionUtilityType utilLocal_;
    GlobalDofConversionUtilityType utilGlobal_;
    //int oldSize_, size_;

  public:
    //! Constructor
    explicit CombinedMapper ( ContainedMapperType &mapper );

  private:
    // prohibit copying
    CombinedMapper ( const ThisType & );

  public:
    //! Total number of degrees of freedom
    int size () const;
    
    /** \copydoc Dune::DofMapper::begin(const ElementType &entity) const */
    DofMapIteratorType begin ( const ElementType &entity ) const;
    
    /** \copydoc Dune::DofMapper::end(const ElementType &entity) const */
    DofMapIteratorType end ( const ElementType &entity ) const;

    /** \copydoc Dune::DofMapper::mapToGlobal(const ElementType &entity,const int localDof) const */
    int mapToGlobal( const ElementType &entity, const int localDof ) const;

    /** \copydoc DofMapper::mapEach */
    template< class Functor >
    void mapEach ( const ElementType &element, Functor f ) const;

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal(const Entity &entity,const int localDof) const */
    template< class Entity >
    int mapEntityDofToGlobal( const Entity &entity, const int localDof ) const;

    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    int maxNumDofs () const;

    /** \copydoc Dune::DofMapper::numDofs(const ElementType &entity) const */
    int numDofs ( const ElementType &entity ) const;

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const;

    //! return new index in dof array 
    int newIndex ( const int hole, const int block ) const;

    //! return old index in dof array of given index ( for dof compress ) 
    int oldIndex ( const int hole, const int block ) const;

    //! return number of holes in the data 
    int numberOfHoles ( const int block ) const;
  
    //! returnn number of mem blocks 
    int numBlocks () const; 

    //! return current old offset of block 
    int oldOffSet ( const int block ) const;

    //! return current offset of block 
    int offSet ( const int block ) const;

    //! return true if compress will affect data  
    bool consecutive () const;

  protected:
    ContainedMapperType &containedMapper () const;

    static int chooseSize ( int pointBased, int variableBased,
                            integral_constant<int,PointBased> );

    static int chooseSize ( int pointBased, int variableBased,
                            integral_constant<int,VariableBased> );

  }; // end class CombinedMapper
  
  /** @} **/  



  // CombinedMapper::FunctorWrapper
  // ------------------------------

  template< class ContainedMapper, int N, DofStoragePolicy policy >
  template< class Functor >
  struct CombinedMapper< ContainedMapper, N, policy >::FunctorWrapper
  {
    FunctorWrapper ( const GlobalDofConversionUtilityType &dofUtil, Functor functor )
    : dofUtil_( dofUtil ),
      functor_( functor )
    {}

    void operator() ( int localBlock, int globalBlock ) 
    {
      int localDof = localBlock*numComponents;
      for( int component = 0; component < numComponents; ++component, ++localDof )
        functor_( localDof, dofUtil_.combinedDof( globalBlock, component ) );
    }

  private:
    GlobalDofConversionUtilityType dofUtil_;
    Functor functor_;
  };



  // Implementation of CombinedMapper
  // --------------------------------
  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline CombinedMapper< ContainedMapper, N, policy >
    :: CombinedMapper ( ContainedMapperType &mapper )
  : containedMapper_( &mapper ),
    utilLocal_( containedMapper(), numComponents ),
    utilGlobal_( containedMapper(), numComponents )
    //oldSize_( spc_.size() ),
    //size_( spc_.size() )
  {}


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy > :: size() const
  {
    return numComponents * containedMapper().size();
  }


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline typename CombinedMapper< ContainedMapper, N, policy > ::  DofMapIteratorType
  CombinedMapper< ContainedMapper, N, policy >
    :: begin ( const ElementType &entity ) const
  {
    return DofMapIteratorType( containedMapper().begin( entity ), utilGlobal_ );
  }


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline typename CombinedMapper< ContainedMapper, N, policy > ::  DofMapIteratorType
  CombinedMapper< ContainedMapper, N, policy >
    :: end ( const ElementType &entity ) const
  {
    return DofMapIteratorType( containedMapper().end( entity ), utilGlobal_ );
  }


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  template< class Functor >
  inline void CombinedMapper< ContainedMapper, N, policy >
    ::mapEach ( const ElementType &element, Functor f ) const
  {
    containedMapper().mapEach( element, FunctorWrapper< Functor >( utilGlobal_, f ) );
  }

  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: mapToGlobal ( const ElementType &entity, const int localDof ) const 
  {
    const int component = utilLocal_.component( localDof );
    const int containedLocal = utilLocal_.containedDof( localDof );
 
    const int containedGlobal
      = containedMapper().mapToGlobal( entity, containedLocal );
    
    return utilGlobal_.combinedDof( containedGlobal, component );
  }


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  template< class Entity >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: mapEntityDofToGlobal ( const Entity &entity, int localDof ) const 
  {
    const int component = utilLocal_.component( localDof );
    const int containedLocal = utilLocal_.containedDof( localDof );
 
    const int containedGlobal = 
      containedMapper().mapEntityDofToGlobal( entity, containedLocal );
    
    return utilGlobal_.combinedDof( containedGlobal, component );
  }


  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: maxNumDofs () const
  {
    return numComponents * containedMapper().maxNumDofs();
  }


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: numDofs ( const ElementType &entity ) const
  {
    return numComponents * containedMapper().numDofs( entity );
  }


  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  template< class Entity >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: numEntityDofs ( const Entity &entity ) const
  {
    return numComponents * containedMapper().numEntityDofs( entity );
  }

  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: newIndex ( const int hole, const int block ) const
  {
    if( policy == PointBased )
    {
      const int component = utilGlobal_.component( hole );
      const int containedHole = utilGlobal_.containedDof( hole );
      const int containedIndex = containedMapper().newIndex( containedHole, block );
      return utilGlobal_.combinedDof( containedIndex, component );
    }
    else
    {
      const int numContainedBlocks = containedMapper().numBlocks();
      const int containedBlock = block % numContainedBlocks;
      const int component = block / numContainedBlocks;
      
      const int containedOffset = containedMapper().newIndex( hole, containedBlock );
      return containedOffset + component * containedMapper().size();
    }
  }


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: oldIndex ( const int hole, const int block ) const
  {
    if( policy == PointBased )
    {
      const int component = utilGlobal_.component( hole );
      const int containedHole = utilGlobal_.containedDof( hole );
      const int containedIndex = containedMapper().oldIndex( containedHole, block );
      return utilGlobal_.combinedDof( containedIndex, component );
    }
    else
    {
      const int numContainedBlocks = containedMapper().numBlocks();
      const int containedBlock = block % numContainedBlocks;
      const int component = block / numContainedBlocks;
      
      const int containedOffset = containedMapper().oldIndex( hole, containedBlock );
      return containedOffset + component * containedMapper().size();
    }

#if 0
    assert(policy != PointBased || block==0);
    return (policy == PointBased) ?
     (mapper_.oldIndex(hole/N,0)*N+hole%N) :
     (mapper_.oldIndex(hole,block%mapper_.numBlocks())+
      oldSize_*block/mapper_.numBlocks());

    DofConversionUtility<policy> 
      tmpUtilGlobal(chooseSize(N, mapper_.size(), integral_constant<int,policy>()));

    const int component = tmpUtilGlobal.component(hole);
    const int contained = tmpUtilGlobal.containedDof(hole);

    // const int containedNew = mapper_.oldIndex(contained,block);
    const int containedNew = mapper_.oldIndex(contained,0);

    return tmpUtilGlobal.combinedDof(containedNew, component);
#endif
  }


  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: numberOfHoles ( const int block ) const
  {
    if( policy == PointBased )
      return numComponents * containedMapper().numberOfHoles( block );
    else
      return containedMapper().numberOfHoles( block % containedMapper().numBlocks() );
  }


  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy > :: numBlocks () const
  {
    const int numContainedBlocks = containedMapper().numBlocks();
    if( policy == PointBased )
      return numContainedBlocks;
    else
      return numComponents * numContainedBlocks;
  }



  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: oldOffSet ( const int block ) const 
  {
    if( policy == PointBased )
      return numComponents * containedMapper().oldOffSet( block );
    else
    {
      const int numContainedBlocks = containedMapper().numBlocks();
      const int containedBlock = block % numContainedBlocks;
      const int component = block / numContainedBlocks;
      
      const int containedOffset = containedMapper().oldOffSet( containedBlock );
      return containedOffset + component * containedMapper().size();
    }
  }



  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: offSet ( const int block ) const
  {
    if( policy == PointBased )
      return numComponents * containedMapper().offSet( block );
    else
    {
      const int numContainedBlocks = containedMapper().numBlocks();
      const int containedBlock = block % numContainedBlocks;
      const int component = block / numContainedBlocks;
      
      const int containedOffset = containedMapper().offSet( containedBlock );
      return containedOffset + component * containedMapper().size();
    }
  }



  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline bool CombinedMapper< ContainedMapper, N, policy >
    :: consecutive () const
  {
    return containedMapper().consecutive();
  }


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline typename CombinedMapper< ContainedMapper, N, policy > :: ContainedMapperType &
  CombinedMapper< ContainedMapper, N, policy > :: containedMapper () const
  {
    return *containedMapper_;
  }

  
#if 0
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: chooseSize ( int pointBased,
                    int variableBased,
                    integral_constant<int, PointBased > )
  {
    return pointBased;
  }


  
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: chooseSize ( int pointBased,
                    int variableBased,
                    integral_constant<int, VariableBased > )
  {
    return variableBased;
  }
#endif

} // end namespace Dune

#endif //#ifndef DUNE_FEM_COMBINEDSPACE_MAPPER_HH

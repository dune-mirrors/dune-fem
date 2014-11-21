#ifndef DUNE_FEM_COMBINEDSPACE_MAPPER_HH
#define DUNE_FEM_COMBINEDSPACE_MAPPER_HH

//- Dune includes
#include <dune/common/math.hh>

//- Local includes
#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/space/combinedspace/combineddofstorage.hh>

namespace Dune
{

  namespace Fem
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

      typedef typename ContainedMapperType :: SizeType SizeType;
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
    : public AdaptiveDofMapper< CombinedMapperTraits< ContainedMapper, N, policy > >
    {
      typedef CombinedMapper< ContainedMapper, N, policy > ThisType;
      typedef AdaptiveDofMapper< CombinedMapperTraits< ContainedMapper, N, policy > > BaseType;

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

      // local dof ordering is always point based !!!
      typedef PointBasedDofConversionUtility< numComponents > LocalDofConversionUtilityType;

    protected:
      //- Data members
      ContainedMapperType *containedMapper_;

      const LocalDofConversionUtilityType utilLocal_;
      GlobalDofConversionUtilityType utilGlobal_;

      int containedSize_ ;
      int containedOldSize_ ;

    public:
      //! Constructor
      explicit CombinedMapper ( ContainedMapperType &mapper );

    private:
      // prohibit copying
      CombinedMapper ( const ThisType & );

    public:
      //! Total number of degrees of freedom
      int size () const;

      /** \copydoc Dune::DofMapper::contains( const int codim ) const */
      bool contains (const int codim ) const
      {
        return containedMapper().contains( codim );
      }

      /** \copydoc Dune::DofMapper::fixedDataSize( const int codim ) const */
      bool fixedDataSize( const int codim ) const
      {
        return containedMapper().fixedDataSize( codim );
      }

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
      template< class Entity, class Functor >
      void mapEachEntityDof ( const Entity &entity, Functor f ) const;

      /** \copydoc Dune::DofMapper::maxNumDofs() const */
      int maxNumDofs () const;

      /** \copydoc Dune::DofMapper::numDofs(const ElementType &element) const */
      int numDofs ( const ElementType &element ) const;

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

      //! return true if compress will affect data (always true for VariableBased)
      bool consecutive () const;

      template< class Entity >
      void insertEntity ( const Entity &entity ) { update(); }

      template< class Entity >
      void removeEntity ( const Entity &entity ) { }

      void resize () { update(); }

      bool compress () { update(); return true; }

      template <class StreamTraits>
      void write( OutStreamInterface< StreamTraits >& out ) const {}

      template <class StreamTraits>
      void read( InStreamInterface< StreamTraits >& in )
      {
        update();
      }

      void backup () const {}
      void restore () {}

    protected:
      void update ()
      {
        // store sizes for old-new index mapping
        containedOldSize_ = containedSize_ ;
        containedSize_    = containedMapper().size();
      }

      ContainedMapperType &containedMapper () const;

    }; // end class CombinedMapper

    /** @} **/



    // CombinedMapper::FunctorWrapper
    // ------------------------------

    template< class ContainedMapper, int N, DofStoragePolicy policy >
    template< class Functor >
    struct CombinedMapper< ContainedMapper, N, policy >::FunctorWrapper
    {
      FunctorWrapper ( const GlobalDofConversionUtilityType &dofUtil,
                       Functor functor )
      : dofUtil_( dofUtil ),
        functor_( functor )
      {}

      template< class GlobalKey >
      void operator() ( int localBlock, const GlobalKey &globalKey )
      {
        int localDof = localBlock*numComponents;
        for( int component = 0; component < numComponents; ++component, ++localDof )
          functor_( localDof, dofUtil_.combinedDof( globalKey, component ) );
      }

      template< class GlobalKey >
      void operator() ( const GlobalKey &globalKey )
      {
        for( int component = 0; component < numComponents; ++component )
          functor_( dofUtil_.combinedDof( globalKey, component ) );
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
      utilLocal_( numComponents ),
      utilGlobal_( containedMapper(), numComponents ),
      containedSize_( containedMapper().size() ),
      containedOldSize_( containedSize_ )
    {
    }


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
      containedMapper().mapEach( element,
          FunctorWrapper< Functor >( utilGlobal_, f ) );
    }


    template< class ContainedMapper, int N, DofStoragePolicy policy >
    inline int CombinedMapper< ContainedMapper, N, policy >
      :: mapToGlobal ( const ElementType &entity, const int localDof ) const
    {
      const int component      = utilLocal_.component( localDof );
      const int containedLocal = utilLocal_.containedDof( localDof );

      const int containedGlobal
        = containedMapper().mapToGlobal( entity, containedLocal );

      return utilGlobal_.combinedDof( containedGlobal, component );
    }


    template< class ContainedMapper, int N, DofStoragePolicy policy >
    template< class Entity, class Functor >
    inline void CombinedMapper< ContainedMapper, N, policy >
      ::mapEachEntityDof ( const Entity &entity, Functor f ) const
    {
      containedMapper().mapEachEntityDof( entity,
            FunctorWrapper< Functor >( utilGlobal_, f ) );
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
      :: newIndex ( const int hole, const int block ) const
    {
      if( policy == PointBased )
      {
        const int component      = utilGlobal_.component( hole );
        const int containedHole  = utilGlobal_.containedDof( hole );
        const int containedIndex = containedMapper().newIndex( containedHole, block );
        return utilGlobal_.combinedDof( containedIndex, component );
      }
      else
      {
        const int numContainedBlocks = containedMapper().numBlocks();
        const int containedBlock     = block % numContainedBlocks;
        const int component          = block / numContainedBlocks;

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
        const int component      = utilGlobal_.component( hole );
        const int containedHole  = utilGlobal_.containedDof( hole );
        const int containedIndex = containedMapper().oldIndex( containedHole, block );
        return utilGlobal_.combinedDof( containedIndex, component );
      }
      else
      {
        const int numContainedBlocks = containedMapper().numBlocks();
        const int containedBlock     = block % numContainedBlocks;
        const int component          = block / numContainedBlocks;

        const int containedOffIndex = containedMapper().oldIndex( hole, containedBlock );
        return containedOffIndex + component * containedOldSize_;
      }
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
        const int containedBlock     = block % numContainedBlocks;
        const int component          = block / numContainedBlocks;

        const int containedOffset = containedMapper().oldOffSet( containedBlock );
        return containedOffset + component * containedOldSize_;
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
      if( policy == PointBased )
        return containedMapper().consecutive();
      else
        // since we have different blocks it's true no matter what
        return true;
    }


    template< class ContainedMapper, int N, DofStoragePolicy policy >
    inline typename CombinedMapper< ContainedMapper, N, policy > :: ContainedMapperType &
    CombinedMapper< ContainedMapper, N, policy > :: containedMapper () const
    {
      assert( containedMapper_ );
      return *containedMapper_;
    }

  } // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_COMBINEDSPACE_MAPPER_HH

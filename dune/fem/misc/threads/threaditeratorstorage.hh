#ifndef DUNE_FEM_THREADITERATORSTORAGE_HH
#define DUNE_FEM_THREADITERATORSTORAGE_HH

#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>

#ifdef USE_SMP_PARALLEL
#include <dune/fem/misc/threads/threadpartitioner.hh>
#endif

namespace Dune {

  namespace Fem {

    /** \brief Storage of thread iterators using domain decomposition */
    template < class ThreadIterator >
    class ThreadIteratorStorageBase
    {
    public:
      typedef ThreadIterator ThreadIteratorType ;
      typedef typename ThreadIterator :: GridPartType GridPartType;
      typedef typename GridPartType :: IndexSetType   IndexSetType;

      typedef typename ThreadIteratorType :: FilterType    FilterType ;
      typedef typename ThreadIteratorType :: IteratorType  IteratorType;

      typedef typename IteratorType :: Entity EntityType ;

      static const PartitionIteratorType pitype = ThreadIteratorType :: pitype ;

    private:
      struct IteratorFactory
      {
        struct Key
        {
          const GridPartType& gridPart_;
          const IndexSetType& indexSet_;
          static const PartitionIteratorType ptype = pitype ;
          Key(const GridPartType& gridPart)
           : gridPart_( gridPart ),
             indexSet_( gridPart_.indexSet() )
          {}

          bool operator ==( const Key& other ) const
          {
            // compare grid pointers
            return (&indexSet_) == (& other.indexSet_ ) && ( ptype == other.ptype );
          }
          const GridPartType& gridPart() const { return gridPart_; }
        };

        typedef ThreadIteratorType ObjectType;
        typedef Key KeyType;

        inline static ObjectType *createObject ( const KeyType &key )
        {
          return new ObjectType( key.gridPart() );
        }

        inline static void deleteObject ( ObjectType *object )
        {
          delete object;
        }
      };


     typedef typename IteratorFactory :: KeyType KeyType;
     typedef SingletonList< KeyType,
             ThreadIteratorType, IteratorFactory > IteratorProviderType;

    protected:
      std::unique_ptr< ThreadIteratorType, typename IteratorProviderType::Deleter> iterators_;

    public:
      //! contructor creating thread iterators
      explicit ThreadIteratorStorageBase( const GridPartType& gridPart )
        : iterators_( &IteratorProviderType::getObject( KeyType( gridPart ) ) )
      {
        update();
      }

      ThreadIteratorType& iterators () const { assert( iterators_ ); return *iterators_; }

      //! return filter for given thread
      const FilterType& filter( const int thread ) const
      {
        return iterators().filter( thread );
      }

      //! update internal list of iterators
      void update()
      {
        iterators().update();
      }

      //! set ratio between master thread and other threads in comp time
      void setMasterRatio( const double ratio )
      {
        iterators().setMasterRatio( ratio );
      }

      //! return begin iterator for current thread
      IteratorType begin() const
      {
        return iterators().begin();
      }

      //! return end iterator for current thread
      IteratorType end() const
      {
        return iterators().end();
      }

      //! return thread number this entity belongs to
      int index(const EntityType& entity ) const
      {
        return iterators().index( entity );
      }

      //! return thread number this entity belongs to
      int thread(const EntityType& entity ) const
      {
        return iterators().thread( entity );
      }
    };
  } // end namespace Fem
} // end namespace Dune

#endif // #ifndef DUNE_FEM_DG_DOMAINTHREADITERATOR_HH

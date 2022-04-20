#ifndef DUNE_FEM_COMMINDEXMAP_HH
#define DUNE_FEM_COMMINDEXMAP_HH

#include <set>
#include <vector>

#include <dune/fem/storage/dynamicarray.hh>

namespace Dune
{

  namespace Fem
  {

    class CommunicationIndexMap
    {
    public:
      // this typedef will also change the type of index
      // for AuxiliaryDofs and AdaptiveLeafIndexSet
      // TODO: extract to a more common place
      typedef int IndexType ;
    private:
      DynamicArray< IndexType > indices_;

    public:
      //! constructor creating empty map
      CommunicationIndexMap()
      : indices_( 0 )
      {
        indices_.setMemoryFactor( 1.1 );
      }

      CommunicationIndexMap( const CommunicationIndexMap& ) = delete;

      //! return index map for entry i
      const IndexType& operator [] ( const size_t i ) const
      {
        assert( i < size() );
        return indices_[ i ];
      }

      //! clear index map
      void clear()
      {
        resize( 0 );
      }

      //! append index vector with idx
      //! result is unsorted
      template <class GlobalKey>
      void insert( const std :: vector< GlobalKey > &idx )
      {
        const size_t size = idx.size();
        size_t  count = indices_.size();

        // reserve memory
        resize( count + size );
        assert( indices_.size() == (count + size) );

        // copy indices to index vector
        for( size_t i = 0; i < size; ++i, ++count )
        {
          assert( idx[ i ] >= 0 );
          indices_[ count ] = idx[ i ];
        }
      }

      //! insert sorted set of indices
      template <class GlobalKey>
      void set( const std :: set< GlobalKey > &idxSet )
      {
        // resize to given new size
        resize( idxSet.size() );

        // copy all elements from set to array
        size_t count = 0;
        typedef typename std :: set< GlobalKey > :: const_iterator iterator;
        const iterator end = idxSet.end();
        for(iterator it = idxSet.begin(); it != end; ++it, ++count)
        {
          indices_[count] = *it;
        }
      }

      //! return size of map
      size_t size () const
      {
        return indices_.size();
      }

      //! print  map for debugging only
      void print( std :: ostream &s, int rank ) const
      {
        const size_t size = this->size();
        s << "Start print: size = " << size << std :: endl;
        for( size_t i = 0; i < size; ++i )
          s << rank << " idx[ " << i << " ] = " << indices_[ i ] << std :: endl;
        s << "End of Array" << std :: endl;
      }

      //! write all indices to buffer
      template <class CommBuffer>
      void writeToBuffer(CommBuffer& buffer) const
      {
        const size_t idxSize = indices_.size();
        buffer.write( idxSize );
        //std::cout << "P[" << MPIManager ::rank() << " write Buffer size " << idxSize << std::endl;
        for(size_t i=0; i<idxSize; ++i)
        {
          //std::cout << "P[" << MPIManager ::rank() << " write idx " << indices_[i] << std::endl;
          buffer.write( indices_[i] );
        }
      }

      //! read all indices from buffer
      template <class CommBuffer>
      void readFromBuffer(CommBuffer& buffer)
      {
        size_t idxSize;
        buffer.read( idxSize );
        //std::cout << "P[" << MPIManager ::rank() << " read Buffer size " << idxSize << std::endl;
        indices_.resize( idxSize );
        for(size_t i=0; i<idxSize; ++i)
        {
          buffer.read( indices_[i] );
          //std::cout << "P[" << MPIManager ::rank() << " read idx " << indices_[i] << std::endl;
        }
      }

    protected:
      //! resize map with size size
      inline void resize ( size_t size )
      {
        indices_.resize( size );
      }

      inline void reserve ( size_t size )
      {
        indices_.reserve( size );
      }

    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMINDEXMAP_HH

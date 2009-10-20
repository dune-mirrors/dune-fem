#ifndef DUNE_FEM_COMMINDEXMAP_HH
#define DUNE_FEM_COMMINDEXMAP_HH

//- system includes
#include <set>
#include <vector>

//- Dune includes
#include <dune/fem/space/common/arrays.hh>

namespace Dune
{

  class CommunicationIndexMap
  {
  private:
    MutableArray< int > indices_;

  public:
    //! constructor creating empty map
    CommunicationIndexMap()
    : indices_( 0 )
    {
      indices_.setMemoryFactor( 1.1 );
    }

  private:
    // prohibit copying
    CommunicationIndexMap( const CommunicationIndexMap &i );

  public:
    //! return index map for entry i
    const int operator [] ( const int i ) const
    {
      assert( (i >= 0) && (i < size()) );
      return indices_[ i ];
    }

    //! clear index map 
    void clear() 
    {
      resize( 0 );
    }

    //! append index vector with idx 
    //! result is unsorted 
    void insert( const std :: vector< int > &idx )
    {
      const int size = idx.size();
      int count = indices_.size();

      // reserve memory 
      resize( count + size );
      assert( indices_.size() == (count + size) );

      // copy indices to index vector 
      for( int i = 0; i < size; ++i, ++count )
      {
        assert( idx[ i ] >= 0 );
        indices_[ count ] = idx[ i ];
      }
    }

    //! insert sorted set of indices  
    void set( const std :: set<int> &idxSet )
    {
      // resize to given new size 
      resize( idxSet.size() );
      
      // copy all elements from set to array 
      int count = 0;
      typedef std :: set<int> :: const_iterator iterator; 
      const iterator end = idxSet.end();
      for(iterator it = idxSet.begin(); it != end; ++it, ++count) 
      {
        indices_[count] = *it;
      }
    }

    //! return size of map
    int size () const
    {
      return indices_.size();
    }

    //! print  map for debugging only 
    void print( std :: ostream &s, int rank ) const
    {
      const int size = this->size();
      s << "Start print: size = " << size << std :: endl;
      for( int i = 0; i < size; ++i )
        s << rank << " idx[ " << i << " ] = " << indices_[ i ] << std :: endl;
      s << "End of Array" << std :: endl;
    }

    //! write all indices to buffer
    template <class CommBuffer> 
    void writeToBuffer(CommBuffer& buffer) const 
    {
      const int idxSize = indices_.size(); 
      buffer.write( indices_.size() );
      for(int i=0; i<idxSize; ++i) 
      {
        buffer.write( indices_[i] );
      }
    }

    //! read all indices from buffer
    template <class CommBuffer> 
    void readFromBuffer(CommBuffer& buffer) 
    {
      int idxSize; 
      buffer.read( idxSize );
      indices_.resize( idxSize );
      for(int i=0; i<idxSize; ++i) 
      {
        buffer.read( indices_[i] );
      }
    }

  protected:  
    //! resize map with size size  
    inline void resize ( int size ) 
    {
      indices_.resize( size );
    }

    inline void reserve ( int size )
    {
      indices_.reserve( size );
    }
    
  };

}

#endif

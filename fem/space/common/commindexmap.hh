#ifndef DUNE_FEM_COMMINDEXMAP_HH
#define DUNE_FEM_COMMINDEXMAP_HH

//- system includes
#include <vector>

//- Dune includes
#include <dune/fem/space/common/arrays.hh>

namespace Dune
{

  class CommunicationIndexMap
  {
  protected:
    MutableArray< int > index_;

  public:
    //! constructor creating empty map
    CommunicationIndexMap()
    : index_( 0 )
    {
      index_.setMemoryFactor( 1.1 );
    }

  private:
    // prohibit copying
    CommunicationIndexMap( const CommunicationIndexMap &i );

  public:
    //! reserve memory 
    void reserve( int size ) 
    {
      // resize array, memory factor will be used 
      index_.resize( size );
    }

    //! clear index map 
    void clear() 
    {
      // resize 0 will free memory 
      index_.resize( 0 );
    }

    //! append index vector with idx 
    void insert( const std :: vector< int > &idx )
    {
      const int size = idx.size();
      int count = index_.size();
      
      // reserve memory 
      reserve( count + size );
      assert( index_.size() == (count + size) );

      // copy indices to index vector 
      for( int i = 0; i < size; ++i, ++count )
      { 
        assert( idx[ i ] >= 0 );
        index_[ count ] = idx[ i ];
      }
    }

    //! return index map for entry i
    const int operator [] ( const int i ) const
    {
      assert( (i >= 0) && (i < size()) );
      return index_[ i ];
    }

    //! return size of map
    int size () const
    {
      return index_.size();
    }

    //! print  map for debugging only 
    void print( std :: ostream &s, int rank ) const
    {
      const int size = index_.size();
      s << "Start print: size = " << size << std :: endl;
      for( int i = 0; i < size; ++i )
        s << rank << " idx[ " << i << " ] = " << index_[ i ] << std :: endl;
      s << "End of Array" << std :: endl;
    }

    void sort() 
    {
      std :: sort( index_.begin(), index_.end() );
    }
  };

}

#endif

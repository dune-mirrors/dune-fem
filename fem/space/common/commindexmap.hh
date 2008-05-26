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

    void insert ( const int &idx )
    {
      const int size = this->size();
      int pos = 0;
      for( ; pos < size; ++pos )
      {
        const int current = indices_[ pos ];
        if( current >= idx )
        {
          if( current == idx )
            return;
          else
            break;
        }
      }

      resize( size + 1 );
      for( int i = size; i > pos; --i )
        indices_[ i ] = indices_[ i-1 ];
      indices_[ pos ] = idx;
    }

    //! append index vector with idx 
    void insert( const std :: vector< int > &idx )
    {
      const int count = idx.size();
      reserve( count + size() );
      for( int i = 0; i < count; ++i )
        insert( idx[ i ] );
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

    void sort () DUNE_DEPRECATED
    {}

  private:
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

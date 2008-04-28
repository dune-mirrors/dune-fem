#ifndef DUNE_FEM_ADAPTIVEINDEXMAPPER_INLINE_HH
#define DUNE_FEM_ADAPTIVEINDEXMAPPER_INLINE_HH

#include "adaptiveindexmapper.hh"

namespace Dune
{

  bool AdaptiveIndexMapper :: compress ( HoleArrayType &holes )
  {
    const unsigned int domainSize = this->domainSize();
    std :: cout << "compressing..." << std :: endl;

    holes.resize( domainSize );
    unsigned int numHoles = 0;
    
    unsigned int newSize = 0;
    for( unsigned int i = 0; i < domainSize; ++i )
    {
      char &state = state_[ i ];

      std :: cout << "state: " << (int)state << ", index: " << index_[ i ]
                  << std :: endl;
      
      if( state == Deleted )
      {
        holes[ numHoles++ ].setNewIndex( index_[ i ] );
        state = Unused;
      }
      else if( state > 0 )
      {
        ++newSize;
        state = Used;
      }
    }
    assert( newSize + numHoles == size_ );
    std :: cout << "newSize: " << newSize << ", numHoles = " << numHoles << std :: endl;

    // remove holes with index greater newSize
    for( unsigned int i = 0; i < numHoles; )
    {
      if( holes[ i ].newIndex() >= newSize )
      {
        std :: cout << "removing holes[ " << i << " ]: " << holes[ i ].newIndex()
                    << std :: endl;
        holes[ i ] = holes[ --numHoles ];
      }
      else
        ++i;
    }
    std :: cout << "numHoles = " << numHoles << std :: endl;

    // fill the holes
    if( numHoles > 0 )
    {
      unsigned int hole = 0;
      for( unsigned int i = 0; i < domainSize; ++i )
      {
        const char state = state_[ i ];
        if( state == Unused )
          continue;
        
        unsigned int &index = index_[ i ];
        std :: cout << "current mapping for " << i << ": " << index << std :: endl;
        if( index < newSize )
          continue;

        assert( hole < numHoles );
        holes[ hole ].setOldIndex( index );
        index = holes[ hole ].newIndex();
        std :: cout << "mapping index: " << holes[ hole ].oldIndex() << "->"
                    << holes[ hole ].newIndex() << std :: endl;
        ++hole;
      }
    }
    size_ = newSize;

    for( unsigned int i = 0; i < domainSize; ++i )
      std :: cout << "state: " << (int)(state_[ i ]) << ", index: " << index_[ i ]
                  << std :: endl;

    // resize the holes array
    holes.resize( numHoles );
    return (numHoles != 0);
  }
  
}

#endif

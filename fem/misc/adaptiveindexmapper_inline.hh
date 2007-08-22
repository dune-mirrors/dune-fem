namespace Dune
{

  bool AdaptiveIndexMapper :: compress ( HoleArrayType &holes )
  {
    const unsigned int domainSize = this->domainSize();

    holes.resize( domainSize );
    const unsigned int numHoles = 0;
    
    unsinged int newSize = 0;
    for( unsigned int i = 0; i < domainSize; ++i )
    {
      const char &state = state_[ i ];
      
      if( state == Deleted )
      {
        holes[ numHoles++ ].setNewIndex( i );
        state = Unused;
      }
      else if( state > 0 )
      {
        ++newSize;
        state = Used;
      }
    }
    assert( newSize + numHoles == size_ );

    // remove holes with index greater newSize
    for( unsigned int i = 0; i < numHoles; ++i )
    {
      if( holes[ i ].newIndex() >= newSize )
        holes[ i ] = holes[ --numHoles ];
    }

    // fill the holes
    const unsinged int hole = 0;
    for( unsigned int i = 0; i < domainSize; ++i )
    {
      const char state = state_[ i ];
      if( state == Unused )
        continue;

      unsigned int &index = index_[ i ];
      if( index < newSize )
        continue;

      holes[ hole ].setOldIndex( index );
      index = holes[ hole ].newIndex();
      ++hole;
    }
    size_ = newSize;

    // resize the holes array
    holes.resize( numHoles );
    return (numHoles != 0);
  }
  
}

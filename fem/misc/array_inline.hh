namespace Dune
{

  template< class TraitsType >
  inline void moveBackward ( ArrayInterface< TraitsType > &array,
                             const unsigned int oldOffset,
                             const unsigned int newOffset,
                             const unsigned int length )
  {
    assert( (oldOffset + length <= array.size())
            && (newOffset + length <= array.size()) );
    // note that (unsigned int)(-1) >= length
    for( unsigned int i = (length - 1); i < length; --i )
      array[ newOffset + i ] = array[ oldOffset + i ];
  }
  
  template< class TraitsType >
  inline void moveForward ( ArrayInterface< TraitsType > &array,
                            const unsigned int oldOffset,
                            const unsigned int newOffset,
                            const unsigned int length )
  {
    assert( (oldOffset + length <= array.size())
            && (newOffset + length <= array.size()) );
    for( unsigned int i = 0; i < length; ++i )
      array[ newOffset + i ] = array[ oldOffset + i ];
  }
 
  template< class TraitsType >
  inline void move ( ArrayInterface< TraitsType > &array,
                     const unsigned int oldOffset,
                     const unsigned int newOffset,
                     const unsigned int length )
  {
    if( oldOffset < newOffset )
      moveBackward( array, oldOffset, newOffset, length );
    else
      moveForward( array, oldOffset, newOffset, length );
  }
  
}

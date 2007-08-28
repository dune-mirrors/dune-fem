namespace Dune
{

  //! Scalar product of two vectors
  template< class Traits1, class Traits2 >
  inline typename ExtractCommonFieldType< Traits1, Traits2 > :: FieldType
  operator* ( const VectorInterface< Traits1 > &u,
              const VectorInterface< Traits2 > &v )
  {
    typedef typename ExtractCommonFieldType< Traits1, Traits2 > :: FieldType FieldType;

    const unsigned int size = u.size();
    assert( size == v.size() );

    FieldType result = 0;
    for( unsigned int i = 0; i < size; ++i )
      result += u[ i ] * v[ i ];
    return result;
  }



  //! Print any Vector into a stream
  template< class Traits >
  inline std :: ostream &operator<< ( std :: ostream &out,
                                      const VectorInterface< Traits > &v )
  {
    const unsigned int size = v.size();
      
    if( size > 0 ) {
      out << "[ ";
      out << v[ 0 ];
      for( unsigned int i = 1; i < size; ++i )
        out << ", " << v[ i ];
      out << " ]";
    } else
      out << "[]";
    return out;
  }



  void match ( std :: istream &in, char excepted )
  {
    char c = 0;
    in >> c;
    if( c != excepted )
      in.clear( std :: ios_base :: badbit );
  }



  //! read any vector from a stream
  template< class Traits >
  inline std :: istream &operator>> ( std :: istream &in,
                                      VectorInterface< Traits > &v )
  {
    const unsigned int size = v.size();
    
    DynamicVector< typename Traits :: FieldType > w( size );

    char expected = '[';
    for( unsigned int i = 0; i < size; ++i )
    {
      match( in, expected );
      in >> w[ i ];
      expected = ',';
    }
    match( in, ']' );

    if( in )
      v.assign( w );
    return in;
  }

}

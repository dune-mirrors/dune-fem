namespace Dune
{

  //! Scalar product of two vectors
  template< class TraitsType, class TraitsOther >
  inline typename TraitsType :: FieldType
  operator* ( const VectorInterface< TraitsType > &u,
              const VectorInterface< TraitsOther > &v )
  {
    typedef typename TraitsType :: FieldType FieldType;
    typedef CompileTimeChecker<
      Conversion< FieldType, typename TraitsOther :: FieldType > :: sameType >
      __FieldType_Of_Both_Vectors_Must_Be_Identical__;

    const unsigned int size = u.size();
    assert( size == v.size() );

    FieldType result = 0;
    for( unsigned int i = 0; i < size; ++i )
      result += u[ i ] * v[ i ];
    return result;
  }



  //! Print any Vector into a stream
  template< class TraitsType >
  inline std :: ostream &operator<< ( std :: ostream &out,
                                      const VectorInterface< TraitsType > &v )
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
  template< class TraitsType >
  inline std :: istream &operator>> ( std :: istream &in,
                                      VectorInterface< TraitsType > &v )
  {
    const unsigned int size = v.size();
    
    DynamicVector< typename TraitsType :: FieldType > w( size );

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

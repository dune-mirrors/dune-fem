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


  
  // Stream Operators
  // ----------------

  /** \brief write a vector into an output stream
   *  \relates VectorInterface
   *  \relatesalso OutStreamInterface
   *
   *  \param[in]  out  stream to write to
   *  \param[in]  v    vector to write
   *
   *  \returns the ouput stream (for concatenation)
   */
  template< class StreamTraits, class VectorTraits >
  inline OutStreamInterface< StreamTraits > &
    operator<< ( OutStreamInterface< StreamTraits > &out,
                 const VectorInterface< VectorTraits > &v )
  {
    typedef typename VectorInterface< VectorTraits > :: ConstIteratorType
      ConstIteratorType;

    const ConstIteratorType end = v.end();
    for( ConstIteratorType it = v.begin(); it != end; ++it )
      out << *it;
    return out;
  }


  
  /** \brief read a vector from an input stream
   *  \relates VectorInterface
   *  \relatesalso InStreamInterface
   *
   *  \param[in]   in  stream to read from
   *  \param[out]  v   vector to read
   *
   *  \returns the input stream (for concatenation)
   */
  template< class StreamTraits, class VectorTraits >
  inline InStreamInterface< StreamTraits > &
    operator>> ( InStreamInterface< StreamTraits > &in,
                 VectorInterface< VectorTraits > &v )
  {
    typedef typename VectorInterface< VectorTraits > :: IteratorType
      IteratorType;

    const IteratorType end = v.end();
    for( IteratorType it = v.begin(); it != end; ++it )
      in >> *it;
    return in;
  }


  /** \brief write a vector into an STL stream
   *  \relates VectorInterface
   *
   *  \param[in]  out  STL stream to write to
   *  \param[in]  v    vector to write
   *
   *  \returns the STL stream (for concatenation)
   */
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



  /** \brief read a vector from an STL stream
   *  \relates VectorInterface
   *
   *  \param[in]   in  STL stream to read from
   *  \param[out]  v   vector to read
   *
   *  \note This method is STL compilant, i.e. the vector is only assigned
   *        after the read operation is successful. This means that a temporary
   *        vector is created.
   *
   *  \returns the STL stream (for concatenation)
   */
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

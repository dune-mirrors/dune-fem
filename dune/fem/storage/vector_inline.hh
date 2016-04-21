#ifndef DUNE_FEM_VECTOR_INLINE_HH
#define DUNE_FEM_VECTOR_INLINE_HH

#include <iostream>
#include <type_traits>

#include "vector.hh"

namespace Dune
{

  namespace Fem
  {
    //! Scalar product of two vectors
    template< class Traits1, class Traits2 >
    inline typename Traits1 :: FieldType
    operator* ( const VectorInterface< Traits1 > &u, const VectorInterface< Traits2 > &v )
    {
      static_assert( (std::is_same< typename Traits1::FieldType, typename Traits2::FieldType >::value),
                     "FieldType must be identical." );

      const auto size = u.size();
      assert( size == v.size() );

      typename Traits1::FieldType result = 0;
      for( unsigned int i = 0; i < size; ++i )
        result += u[ i ] * v[ i ];
      return result;
    }



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
      operator<< ( OutStreamInterface< StreamTraits > &out, const VectorInterface< VectorTraits > &v )
    {
      for( const auto& entry : v )
        out << entry;
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
      operator>> ( InStreamInterface< StreamTraits > &in, VectorInterface< VectorTraits > &v )
    {
      for( auto& entry : v )
        in >> v;
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
    inline std :: ostream &operator<< ( std :: ostream &out, const VectorInterface< Traits > &v )
    {
      out << "[ ";
      auto it = v.begin();
      const auto itEnd = v.end();
      if( it != itEnd )
      {
        out << *it;
        for( ; it != itEnd ; ++it )
          out << ", " << *it;
      }
      out << " ]";
      return out;
    }



    inline void match ( std :: istream &in, char excepted )
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
    inline std :: istream &operator>> ( std :: istream &in, VectorInterface< Traits > &v )
    {
      DynamicVector< typename Traits :: FieldType > w( v.size() );

      char expected = '[';
      for( auto& entry : w )
      {
        match( in, expected );
        in >> entry;
        expected = ',';
      }
      match( in, ']' );

      if( in )
        v.assign( w );
      return in;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_VECTOR_INLINE_HH

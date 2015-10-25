#ifndef EIGENVECTOR_HH
#define EIGENVECTOR_HH

#if HAVE_EIGEN

#include <dune/fem/storage/vector.hh>
#include <Eigen/Dense>

namespace Dune
{
  namespace Fem
  {
    template< class Field >
    class EigenVector
    : public VectorDefault< Field, EigenVector< Field > >
    {
      typedef EigenVector< Field > ThisType;
      typedef VectorDefault< Field, ThisType > BaseType;
    public:
      //! field type of the vector
      typedef Field FieldType;

      using BaseType :: assign;

    protected:
      typedef Eigen::Matrix<Field, Eigen::Dynamic, 1> DofStorageType;
      DofStorageType fields_;

    public:
      //! Constructor setting up a vector of a specified size
      inline explicit EigenVector ( unsigned int size = 0 )
      : fields_( size )
      {}

      //! Constructor setting up a vector iniitialized with a constant value
      inline EigenVector ( unsigned int size,
                               const FieldType s )
      : fields_( size )
      {
        assign( s );
      }

      //! Copy constructor setting up a vector with the data of another one
      template< class T >
      inline EigenVector ( const VectorInterface< T > &v )
      : fields_()
      {
        assign( v );
      }

      //! Copy constructor setting up a vector with the data of another one (of the same type)
      inline EigenVector ( const ThisType &v )
      : fields_()
      {
        assign( v );
      }

      //! Assign another vector to this one
      template< class T >
      inline ThisType &operator= ( const VectorInterface< T > &v )
      {
        assign( v );
        return *this;
      }

      //! Assign another vector (of the same type) to this one
      inline ThisType &operator= ( const ThisType &v )
      {
        assign( v );
        return *this;
      }

      //! Initialize all fields of this vector with a scalar
      inline ThisType &operator= ( const FieldType s )
      {
        assign( s );
        return *this;
      }

      inline const FieldType &operator[] ( unsigned int index ) const
      {
        return fields_( index );
      }

      inline FieldType &operator[] ( unsigned int index )
      {
        return fields_( index );
      }

      /*
      template< class T >
      inline void assign ( const VectorInterface< T > &v )
      {
        fields_.assign( v );
      }
      */

      inline const DofStorageType &coefficients () const
      {
        return fields_;
      }
      inline DofStorageType &coefficients ()
      {
        return fields_;
      }

      inline const FieldType *leakPointer () const
      {
        return fields_.memptr();
      }

      inline FieldType *leakPointer ()
      {
        return fields_.memptr();
      }

      inline void reserve ( unsigned int newSize )
      {
        fields_.resize( newSize );
      }

      inline void resize ( unsigned int newSize )
      {
        fields_.resize( newSize );
      }

      inline void resize ( unsigned int newSize,
                           const FieldType defaultValue )
      {
        fields_.resize( newSize ); assign( defaultValue );
      }

      inline unsigned int size () const
      {
        return fields_.size();
      }
    };
  }
}
#endif

#endif // EIGENVECTOR_HH


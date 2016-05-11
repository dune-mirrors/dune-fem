#ifndef EIGENVECTOR_HH
#define EIGENVECTOR_HH

#if HAVE_EIGEN

#include <dune/fem/misc/metaprogramming.hh>
#include <dune/fem/storage/vector.hh>

#include <Eigen/Dense>


/*! @addtogroup VectorClasses
    @{
*/

namespace Dune
{
  namespace Fem
  {
    /** \class EigenVector
     *  \ingroup Vector
     *  \brief An implementation of VectorInterface which uses Eigen::Matrix<Field, Eigen::Dynamic, 1> to provide the fields.
     */
    template< class Field >
    class EigenVector
    : public VectorDefault< Field, EigenVector< Field > >
    {
      typedef EigenVector< Field > ThisType;
      typedef VectorDefault< Field, ThisType > BaseType;

    public:
      //! Field type of the vector
      typedef Field FieldType;

      //! DOFs storage type
      typedef Eigen::Matrix<Field, Eigen::Dynamic, 1> DofStorageType;

      using BaseType :: assign;

      //! Constructor setting up a vector of a specified size
      explicit EigenVector ( unsigned int size = 0 )
      : fields_( size )
      {}

      //! Constructor setting up a vector iniitialized with a constant value
      EigenVector ( unsigned int size, const FieldType &s )
      : fields_( size )
      {
        assign( s );
      }

      //! Copy constructor setting up a vector with the data of another one
      template< class T >
      EigenVector ( const VectorInterface< T > &v )
      : fields_()
      {
        assign( v );
      }

      //! Copy constructor setting up a vector with the data of another one (of the same type)
      EigenVector ( const ThisType &v )
      : fields_()
      {
        assign( v );
      }

      /** \copydoc Dune::Fem::VectorInterface::operator=(const ThisType &v) */
      ThisType &operator= ( const ThisType &v )
      {
        assign( v );
        return *this;
      }

      /** \copydoc Dune::Fem::VectorInterface::operator[](unsigned int index) */
      const FieldType &operator[] ( unsigned int index ) const
      {
        return fields_( index );
      }

      /** \copydoc Dune::Fem::VectorInterface::operator[](unsigned int index) */
      FieldType &operator[] ( unsigned int index )
      {
        return fields_( index );
      }

      //! Obtain coefficients
      const DofStorageType &coefficients () const
      {
        return fields_;
      }

      //! Obtain coefficients
      DofStorageType &coefficients ()
      {
        return fields_;
      }

      //! Obtain pointer to data
      const FieldType *leakPointer () const
      {
        return fields_.memptr();
      }

      //! Obtain pointer to data
      FieldType *leakPointer ()
      {
        return fields_.memptr();
      }

      //! Allocate memory
      void reserve ( unsigned int newSize )
      {
        fields_.resize( newSize );
      }

      //! Resize vector
      void resize ( unsigned int newSize )
      {
        fields_.resize( newSize );
      }

      //! Resize vector and fill with the defaultValue
      void resize ( unsigned int newSize, const FieldType &defaultValue )
      {
        fields_.resize( newSize );
        assign( defaultValue );
      }

      /** \copydoc Dune::Fem::VectorInterface::size() */
      unsigned int size () const
      {
        return fields_.size();
      }

    private:
      DofStorageType fields_;
    };

  }
}
#endif

//! @}

#endif // EIGENVECTOR_HH


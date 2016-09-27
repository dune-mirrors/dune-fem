#ifndef EIGENVECTOR_HH
#define EIGENVECTOR_HH

#ifdef HAVE_EIGEN

#include <algorithm>
#include <iostream>

#include <dune/common/densevector.hh>
#include <dune/common/ftraits.hh>

#include <Eigen/Dense>

namespace Dune
{

  namespace Fem
  {
    // forward declaration
    template< class K > class EigenVector;
  }

  // specialization of DenseMatVecTraits for EigenVector
  template< class K >
  struct DenseMatVecTraits< Fem::EigenVector< K > >
  {
    typedef Fem::EigenVector< K > derived_type;
    typedef Eigen::Matrix< K, Eigen::Dynamic, 1 > container_type;
    typedef K value_type;
    typedef unsigned int size_type;
  };

  template< class K >
  struct FieldTraits< Fem::EigenVector< K > >
  {
    typedef typename FieldTraits< K >::field_type field_type;
    typedef typename FieldTraits< K >::real_type real_type;
  };


  namespace Fem
  {

    /** \brief An implementation of DenseVector which uses Eigen::Matrix<K, Eigen::Dynamic, 1> to provide the fields
     *
     * \tparam K is the field type (use float, double, complex, etc)
     */
    template< class K >
    class EigenVector : public DenseVector< EigenVector< K > >
    {
      typedef EigenVector< K > ThisType;
      typedef DenseVector< ThisType > BaseType;

    public:
      typedef typename BaseType::size_type size_type;
      typedef size_type SizeType;
      typedef typename BaseType::value_type value_type;
      typedef value_type FieldType;
      typedef typename DenseMatVecTraits< ThisType >::container_type container_type;
      typedef container_type DofStorageType;

      //! Constructor setting up a vector of a specified size
      explicit EigenVector( size_type size = 0 )
      : data_( size )
      {}

      //! Constructor setting up a vector initialized with a constant value
      EigenVector( size_type size, const value_type& s )
      : data_( size )
      {
        std::fill( data_.begin(), data_.end(), s );
      }

      //! Copy constructor setting up a vector with the data of another one
      template< class T >
      EigenVector( const DenseVector< T >& v )
      : data_( v.size() )
      {
        std::copy( v.begin(), v.end(), data_.begin() );
      }

      //! Copy assignment operator
      EigenVector &operator=( const EigenVector& other )
      {
        data_.resize( other.size() );
        std::copy( other.begin(), other.end(), data_.begin() );
        return *this;
      }

      const value_type& operator[]( size_type index ) const
      {
        return data_( index );
      }

      value_type& operator[]( size_type index )
      {
        return data_( index );
      }

      //! Obtain coefficients
      const DofStorageType& coefficients() const
      {
        return data_;
      }

      //! Obtain coefficients
      DofStorageType& coefficients()
      {
        return data_;
      }

      //! Obtain pointer to data
      const value_type* data() const
      {
        return data_.data();
      }

      //! Obtain pointer to data
      value_type* data()
      {
        return data_.data();
      }

      //! Allocate memory
      void reserve( size_type newSize )
      {
        data_.resize( newSize );
      }

      //! Resize vector
      void resize( size_type newSize )
      {
        data_.resize( newSize );
      }

      //! Resize vector and initialize
      void resize( size_type newSize, const value_type& s )
      {
        data_.resize( newSize );
        std::fill( data_.begin(), data_.end(), s );
      }

      size_type size() const
      {
        return data_.size();
      }

    private:
      DofStorageType data_;
    };

    /** \brief Read a EigenVector from an input stream
     *  \relates EigenVector
     *
     *  \note This operator is STL compilant, i.e., the content of v is only
     *        changed if the read operation is successful.
     *
     *  \param[in]  in  std::istream to read from
     *  \param[out] v   EigenVector to be read
     *
     *  \returns the input stream (in)
     */
    template< class K >
    inline std::istream& operator>>( std::istream& in, EigenVector< K >& v )
    {
      EigenVector< K > w(v);
      for( typename EigenVector< K >::size_type i = 0; i < w.size(); ++i )
        in >> w[ i ];
      if(in)
        v = w;
      return in;
    }

  } // namespace Fem

} // namespace Dune

#endif

#endif // #ifndef EIGENVECTOR_HH

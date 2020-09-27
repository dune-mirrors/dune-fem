#ifndef DUNE_FEM_REFERENCEVECTOR_HH
#define DUNE_FEM_REFERENCEVECTOR_HH

#include <algorithm>
#include <cassert>
#include <memory>
#include <utility>
#include <vector>

#include <dune/common/densevector.hh>
#include <dune/common/ftraits.hh>

#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune
{

  namespace Fem
  {
    // forward declaration
    template< class K, class A = std::allocator< K* > >
    class DynamicReferenceVector;
  }

  // specialization of DenseMatVecTraits for DynamicReferenceVector
  template< class K, class A >
  struct DenseMatVecTraits< Fem::DynamicReferenceVector< K, A > >
  {
    typedef Fem::DynamicReferenceVector< K, A > derived_type;
    typedef std::vector< K*, typename A::template rebind< K* >::other > container_type;

    typedef K value_type;
    typedef typename container_type::size_type size_type;
  };

  template< class K, class A >
  struct FieldTraits< Fem::DynamicReferenceVector< K, A > >
  {
    typedef typename FieldTraits< K >::field_type field_type;
    typedef typename FieldTraits< K >::real_type real_type;
  };

  namespace Fem
  {

    /** \brief An implementation of DenseVector which uses a std::vector of references as storage
      *
      * \tparam T is the field type (use float, double, complex, etc)
      * \tparam A allocator type
      */
    template< class K, class A >
    class DynamicReferenceVector : public DenseVector< DynamicReferenceVector< K, A > >
    {
      typedef DynamicReferenceVector< K, A > ThisType;
      typedef DenseVector< ThisType > BaseType;

    public:
      typedef typename BaseType::size_type size_type;
      typedef typename BaseType::value_type value_type;
      typedef value_type FieldType;
      typedef typename DenseMatVecTraits< ThisType >::container_type DofStorageType;

      //! Constructor with uninitialized vector
      explicit DynamicReferenceVector ( const A &a = A() )
      : data_( a )
      {}

      //! Constructor with uninitialized vector of size n
      explicit DynamicReferenceVector ( size_type n, const A &a = A() )
      : data_( n, nullptr, a )
      {}

      //! Copy constructor
      DynamicReferenceVector ( const ThisType &other )
      : BaseType(), // DenseVector stores no data anyway
        data_( other.data_ )
      {}

      //! Move constructor
      DynamicReferenceVector ( ThisType &&other )
      : data_( std::move( other.data_ ) )
      {}

      // issue with gcc 5.5 and ambiguity of
      //      operator= ( const DenseVector< V > &other )
      // using BaseType::operator=;
      // to avoid this adding to missing operator= overload:
      ThisType& operator= (const value_type& k)
      {
        for (size_type i=0; i<size(); i++)
          data_[i] = k;
        return *this;
      }

      template< class V >
      ThisType & operator= ( const DenseVector< V > &other )
      {
        assert( data_.size() == other.size() );
        std::copy( other.begin(), other.end(), BaseType::begin() );
        return *this;
      }

      ThisType & operator= ( const ThisType &other )
      {
        assert( data_.size() == other.size() );
        std::copy( other.begin(), other.end(), BaseType::begin() );
        return *this;
      }

      ThisType & operator= ( ThisType &&other )
      {
        data_ = std::move( other.data_ );
        return *this;
      }

      size_type capacity () const
      {
        return data_.capacity();
      }

      void resize ( size_type n )
      {
        data_.resize( n, nullptr );
      }

      void reserve ( size_type n )
      {
        data_.reserve( n );
      }

      //! Bind i-th entry to a reference
      void bind ( size_type i, K& u )
      {
        assert( i < data_.size() );
        data_[ i ] = &u;
      }

       //! Unbind i-th entry
      void unbind ( size_type i )
      {
        asssert( i < data_.size() );
        data_[ i ] = nullptr;
      }

      size_type size () const
      {
        return data_.size();
      }

      value_type &operator[] ( size_type i )
      {
        return *data_[ i ];
      }

      const value_type &operator[] ( size_type i ) const
      {
        return *data_[ i ];
      }

      void clear()
      {
        for (std::size_t i=0;i<size();++i)
          (*this)[i] = 0;
      }

    private:
      DofStorageType data_;
    };

  } // namespace Fem


} // namespace Dune

#endif //#ifndef DUNE_FEM_REFERENCEVECTOR_HH

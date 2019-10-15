#ifndef DUNE_FEM_SUBVECTOR_HH
#define DUNE_FEM_SUBVECTOR_HH

#include <algorithm>

#include <dune/common/densevector.hh>
#include <dune/common/ftraits.hh>

#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune
{

  namespace Fem
  {
    // forward declaration
    template< class K, class M > class SubVector;
  }

  // specialization of DenseMatVecTraits for SubVector
  template< class K, class M >
  struct DenseMatVecTraits< Fem::SubVector< K, M > >
  {
    typedef Fem::SubVector< K, M > derived_type;
    typedef K container_type;

    typedef std::decay_t< decltype( std::declval< K >().size() ) > size_type;
    typedef std::decay_t< decltype( std::declval< K >()[ std::declval< size_type >() ] ) > value_type;
  };

  template< class K, class M >
  struct FieldTraits< Fem::SubVector< K, M > >
  {
    typedef typename FieldTraits< typename DenseMatVecTraits< Fem::SubVector< K, M > >::value_type >::field_type field_type;
    typedef typename FieldTraits< typename DenseMatVecTraits< Fem::SubVector< K, M > >::value_type >::real_type real_type;
  };

  namespace Fem
  {

    //! Abstract index mapper interface
    template< class IM >
    class IndexMapperInterface
    : public BartonNackmanInterface< IndexMapperInterface< IM >, IM >
    {
      typedef IndexMapperInterface< IM > ThisType;
      typedef BartonNackmanInterface< ThisType, IM > BaseType;

    public:
      //! Type of the implementation (Barton-Nackman)
      typedef IM IndexMapperType;

      //! Type of the interface
      typedef ThisType IndexMapperInterfaceType;

      //! Maps an index onto another one
      unsigned int operator[] ( unsigned int index ) const
      {
        return asImp().operator[]( index );
      }

      //! Returns the map's range
      unsigned int range () const
      {
        return asImp().range();
      }

      //! Returns the map's size
      unsigned int size () const
      {
        return asImp().size();
      }

    protected:
      using BaseType::asImp;
    };



    //! Index mapper which simply adds an offset to the index
    class OffsetSubMapper
    : public IndexMapperInterface< OffsetSubMapper >
    {
      typedef OffsetSubMapper ThisType;
      typedef IndexMapperInterface< OffsetSubMapper > BaseType;

    public:
      OffsetSubMapper( unsigned int size, unsigned int offset )
      : size_( size ), offset_( offset )
      {}

      OffsetSubMapper( const ThisType& ) = default;
      OffsetSubMapper( ThisType&& ) = default;

      unsigned int size() const
      {
        return size_;
      }

      unsigned int range() const
      {
        return size_;
      }

      unsigned int operator[]( unsigned int i) const
      {
        return i+offset_;
      }

    private:
      const unsigned int size_;
      const unsigned int offset_;
    };



    //! Index mapper with static size which simply adds an offset to the index
    template<unsigned int dim>
    class StaticOffsetSubMapper
    : public IndexMapperInterface< StaticOffsetSubMapper< dim > >
    {
      typedef StaticOffsetSubMapper< dim > ThisType;
      typedef IndexMapperInterface< StaticOffsetSubMapper< dim > > BaseType;

    public:
      StaticOffsetSubMapper( unsigned int offset )
      : offset_( offset )
      {}

      StaticOffsetSubMapper( const ThisType& ) = default;
      StaticOffsetSubMapper( ThisType&& ) = default;

      static constexpr unsigned int size()
      {
        return dim;
      }

      static constexpr unsigned int range()
      {
        return dim;
      }

      unsigned int operator[]( unsigned int i) const
      {
        return i+offset_;
      }

    private:
      const unsigned int offset_;
    };



    /** \brief An implementation of DenseVector to extract a portion, not necessarly contiguos, of a vector
     *
     * \tparam BaseVectorImp The base vector
     * \tparam IndexMapperImp The index mapper
     */
    template< class BaseVectorImp, class IndexMapperImp >
    class SubVector : public DenseVector< SubVector< BaseVectorImp, IndexMapperImp > >
    {
      typedef SubVector< BaseVectorImp, IndexMapperImp > ThisType;
      typedef DenseVector< ThisType > BaseType;

    public:
      typedef typename BaseType::size_type size_type;
      typedef typename BaseType::value_type value_type;

      using BaseType::operator=;

      //! Type of the base vector
      typedef BaseVectorImp BaseVectorType;

      //! Type of the index mapper
      typedef IndexMapperImp IndexMapperType;

      //! Type of vector elements
      typedef value_type FieldType;

      //! Constructor
      explicit SubVector( BaseVectorType& baseVector, IndexMapperType&& indexMapper )
      : baseVector_( baseVector ), indexMapper_( indexMapper )
      {}

      SubVector( const ThisType & other )
        : baseVector_( other.baseVector_ ), indexMapper_( other.indexMapper_ )
      {}

      //! Copy entries
      ThisType& operator=( const ThisType & other)
      {
        std::copy( other.begin(), other.end(), this->begin() );
        return *this;
      }

      void clear()
      {
        std::fill( this->begin(), this->end(), FieldType(0) );
      }

      void resize( size_type )
      {}

      const value_type& operator[]( size_type i ) const
      {
        return baseVector_[ indexMapper_[ i ] ];
      }

      value_type& operator[]( size_type i )
      {
        return baseVector_[ indexMapper_[ i ] ];
      }

      size_type size() const
      {
        return indexMapper_.size();
      }

    private:
      BaseVectorType& baseVector_;
      IndexMapperType indexMapper_;
    };


  } // namespace Fem

} // namespace Dune

#endif

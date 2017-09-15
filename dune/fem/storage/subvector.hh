#ifndef DUNE_FEM_SUBVECTOR_HH
#define DUNE_FEM_SUBVECTOR_HH

#include <cassert>

#include <algorithm>

#include <dune/common/densevector.hh>
#include <dune/common/ftraits.hh>

#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class K, class M >
    class SubVector;

    template< class BaseMatrix, class IndexMapper >
    class SubRowMatrix;

  } // namespace Fem



  // DenseMatVecTraits for SubVector
  // -------------------------------

  template< class K, class M >
  struct DenseMatVecTraits< Fem::SubVector< K, M > >
  {
    typedef Fem::SubVector< K, M > derived_type;
    typedef K container_type;

    typedef std::decay_t< decltype( std::declval< K >().size() ) > size_type;
    typedef std::decay_t< decltype( std::declval< K >()[ std::declval< size_type >() ] ) > value_type;
  };



  // DenseMatVecTraits for SubRowMatrix
  // ----------------------------------

  template< class BaseMatrix, class IndexMapper >
  struct DenseMatVecTraits< Fem::SubRowMatrix< BaseMatrix, IndexMapper > >
  {
    typedef Fem::SubRowMatrix< BaseMatrix, IndexMapper > derived_type;
    typedef BaseMatrix container_type;

    typedef typename DenseMatVecTraits< BaseMatrix >::row_type row_type;

    typedef typename DenseMatVecTraits< BaseMatrix >::row_reference row_reference;
    typedef typename DenseMatVecTraits< BaseMatrix >::const_row_reference const_row_reference;

    typedef typename DenseMatVecTraits< BaseMatrix >::size_type size_type;
    typedef typename DenseMatVecTraits< BaseMatrix >::value_type value_type;
  };

  template< class BaseMatrix, class IndexMapper >
  struct DenseMatVecTraits< Fem::SubRowMatrix< const BaseMatrix, IndexMapper > >
  {
    typedef Fem::SubRowMatrix< BaseMatrix, IndexMapper > derived_type;
    typedef BaseMatrix container_type;

    typedef typename DenseMatVecTraits< BaseMatrix >::row_type row_type;

    typedef typename DenseMatVecTraits< BaseMatrix >::const_row_reference row_reference;
    typedef typename DenseMatVecTraits< BaseMatrix >::const_row_reference const_row_reference;

    typedef typename DenseMatVecTraits< BaseMatrix >::size_type size_type;
    typedef typename DenseMatVecTraits< BaseMatrix >::value_type value_type;
  };



  // FieldTraits for SubVector
  // -------------------------

  template< class K, class M >
  struct FieldTraits< Fem::SubVector< K, M > >
  {
    typedef typename FieldTraits< typename DenseMatVecTraits< Fem::SubVector< K, M > >::value_type >::field_type field_type;
    typedef typename FieldTraits< typename DenseMatVecTraits< Fem::SubVector< K, M > >::value_type >::real_type real_type;
  };



  // FieldTraits for SubRowMatrix
  // ----------------------------

  template< class BaseMatrix, class IndexMapper >
  struct FieldTraits< Fem::SubRowMatrix< BaseMatrix, IndexMapper > >
  {
    typedef typename FieldTraits< typename DenseMatVecTraits< Fem::SubRowMatrix< BaseMatrix, IndexMapper > >::value_type >::field_type field_type;
    typedef typename FieldTraits< typename DenseMatVecTraits< Fem::SubRowMatrix< BaseMatrix, IndexMapper > >::value_type >::real_type real_type;
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
      ThisType& operator=( const ThisType& ) = default;
      ThisType& operator=( ThisType&& ) = default;

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
      ThisType& operator=( const ThisType& ) = default;
      ThisType& operator=( ThisType&& ) = default;

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

      using BaseType::operator=;

      //! Copy entries
      ThisType& operator=( const ThisType & other)
      {
        assert( size() == other.size() );
        std::copy( other.begin(), other.end(), this->begin() );
        return *this;
      }

      template< class T >
      ThisType &operator= ( const DenseVector< T > &other )
      {
        assert( size() == other.size() );
        std::copy( other.begin(), other.end(), this->begin() );
        return *this;
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



    // SubRowMatrix
    // ------------

    template< class BaseMatrix, class IndexMapper >
    class SubRowMatrix
      : public DenseMatrix< SubRowMatrix< BaseMatrix, IndexMapper > >
    {
      typedef SubVector< BaseMatrix, IndexMapper > ThisType;
      typedef DenseMatrix< SubRowMatrix< BaseMatrix, IndexMapper > > BaseType;

    public:
      typedef typename BaseType::size_type size_type;
      typedef typename BaseType::value_type value_type;

      typedef typename BaseType::const_row_reference const_row_reference;
      typedef typename BaseType::row_reference row_reference;

      typedef BaseMatrix BaseMatrixType;
      typedef IndexMapper IndexMapperType;

      explicit SubRowMatrix ( BaseMatrixType &baseMatrix, IndexMapperType &&indexMapper )
        : baseMatrix_( baseMatrix ), indexMapper_( std::move( indexMapper ) )
      {}

      using BaseType::operator=;

      ThisType &operator= ( const ThisType &other )
      {
        assert( mat_rows() == other.N() );
        std::copy( other.begin(), other.end(), this->begin() );
        return *this;
      }

      template< class T >
      ThisType &operator= ( const DenseMatrix< T > &other )
      {
        assert( mat_rows() == other.N() );
        std::copy( other.begin(), other.end(), this->begin() );
        return *this;
      }

      const_row_reference mat_access ( size_type i ) const { return baseMatrix_[ indexMapper_[ i ] ]; }
      row_reference mat_access ( size_type i ) { return baseMatrix_[ indexMapper_[ i ] ]; }

      size_type mat_rows () const { return indexMapper_.size(); }
      size_type mat_cols () const { return baseMatrix_.N(); }

    private:
      BaseMatrixType &baseMatrix_;
      IndexMapperType indexMapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif

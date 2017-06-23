#ifndef DUNE_FEM_STORAGE_ROWREFERENCEVECTOR_HH
#define DUNE_FEM_STORAGE_ROWREFERENCEVECTOR_HH

#include <dune/common/densevector.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class K >
    class RowReferenceVector;

  } // namespace Fem



  // DenseMatVecTraits for RowReferenceVector
  // ----------------------------------------

  template< class K >
  struct DenseMatVecTraits< Fem::RowReferenceVector< K > >
  {
    typedef Fem::RowReferenceVector< K > derived_type;
    typedef K value_type;
    typedef std::size_t size_type;
  };



  // FieldTraits for RowReferenceVector
  // ----------------------------------

  template< class K >
  struct FieldTraits< Fem::RowReferenceVector< K > >
  {
    typedef typename FieldTraits< K >::field_type field_type;
    typedef typename FieldTraits< K >::real_type real_type;
  };



  namespace Fem
  {

    // RowReferenceVector
    // ------------------

    template< class K >
    class RowReferenceVector
      : public Dune::DenseVector< RowReferenceVector< K > >
    {
      typedef Dune::DenseVector< RowReferenceVector< K > > Base;

    public:
      typedef typename Base::size_type size_type;
      typedef typename Base::value_type value_type;

      RowReferenceVector ( K *data, size_type size )
        : data_( data ), size_( size )
      {}

      RowReferenceVector ( const RowReferenceVector &other )
        : data_( other.data_ ), size_( other.size_ )
      {}

      using Base::operator=;

      const K &operator[] ( size_type i ) const { return data_[ i ]; }
      K &operator[] ( size_type i ) { return data_[ i ]; }

      size_type size () const { return size_; }

      const K *data () const { return data_; }
      K *data () { return data_; }

    private:
      K *data_;
      size_type size_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_STORAGE_ROWREFERENCEVECTOR_HH

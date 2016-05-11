#ifndef DUNE_FEM_COMMON_SUBVECTOR_HH
#define DUNE_FEM_COMMON_SUBVECTOR_HH

#include <dune/common/densevector.hh>
#include <dune/common/ftraits.hh>

namespace Dune
{

  namespace Fem
  {
    template< class V >
    class DenseSubVector;
  }


  template< class V >
  struct DenseMatVecTraits< Fem::DenseSubVector< V > >
  {
    typedef Fem::DenseSubVector< V > derived_type;
    typedef V container_type;

    typedef typename V::value_type value_type;
    typedef typename V::size_type size_type;
  };

  template< class V >
  struct FieldTraits< Fem::DenseSubVector< V > >
  {
    typedef typename FieldTraits< typename V::value_type >::field_type field_type;
    typedef typename FieldTraits< typename V::value_type >::real_type real_type;
  };

  namespace Fem
  {

    template< class V >
    class DenseSubVector : public DenseVector< DenseSubVector< V > >
    {
      typedef DenseSubVector< V > This;
      typedef DenseVector< DenseSubVector < V > > Base;

    public:
      typedef typename Base::size_type size_type;
      typedef typename Base::value_type value_type;

      //! Constructor making uninitialized vector
      explicit DenseSubVector ( V& v, size_type size, size_type offset )
      : v_( v ),
        size_( size ),
        offset_( offset )
      {}

      DenseSubVector ( const This &other )
      : v_( other.v_ ),
        size_( other.size_ ),
        offset_( other.offset_ )
      {}

      using Base::operator=;

      void resize ( size_type ) {}

      size_type size () const { return size_; }

      value_type &operator[] ( size_type i ) { return v_[ i + offset_ ]; }
      const value_type &operator[] ( size_type i ) const { return v_[ i + offset_ ]; }

    private:
      typename DenseMatVecTraits< This > :: container_type &v_;
      const size_type size_, offset_;
    };

  } // namespace Fem
} // namespace Dune

#endif //#ifndef DUNE_FEM_COMMON_SUBVECTOR_HH

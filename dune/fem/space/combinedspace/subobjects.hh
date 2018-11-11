#ifndef DUNE_FEM_SUBOBJECTS_HH
#define DUNE_FEM_SUBOBJECTS_HH

#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/fem/common/explicitfieldvector.hh>


namespace Dune
{

  namespace Fem
  {

    template< class T >
    struct RowType;

    template< class T >
    struct RowType< const T>
    {
      typedef const typename RowType<T> :: Type Type;
      static const int size = RowType<T> :: size;
    };

    template< class K, int SIZE >
    struct RowType< FieldVector< K, SIZE > >
    {
      typedef K Type;
      static const int size = SIZE;
    };

    template< class K, int SIZE >
    struct RowType< ExplicitFieldVector< K, SIZE > >
    {
      typedef K Type;
      static const int size = SIZE;
    };

    template< class K, int ROWS, int COLS >
    struct RowType< FieldMatrix< K, ROWS, COLS > >
    {
      typedef FieldVector<K, COLS> Type;
      static const int size = ROWS;
    };



    template <class DomainObject, class RangeObject, int offset >
    class SubObject
    {
      typedef DomainObject DomainObjectType;
      typedef RangeObject RangeObjectType;

      typedef typename RowType< RangeObject > :: Type RowType;

    public:
      SubObject( DomainObjectType &host )
      : host_( host )
      {}

      const RowType &operator[] ( const int i ) const
      {
        assert( (i >=0 ) && (i < size()) );
        return host_[ i + offset ];
      }

      RowType& operator[] ( const int i )
      {
        assert( (i >=0 ) && (i < size()) );
        return host_[ i + offset ];
      }

      int size () const
      {
        return Dune::Fem::RowType< RangeObject > :: size;
      }

      operator typename std::remove_const< RangeObjectType >::type () const
      {
        typename std::remove_const< RangeObjectType >::type y;
        for( int i = 0; i < size(); ++i )
          y[ i ] = (*this)[ i ];
        return y;
      }

    private:
      DomainObjectType &host_;
    };

  } // namespace Fem


  // cast into fieldMatrix
  template< class DomianObj, class RangeObj, int offset >
  struct DenseMatrixAssigner< typename std::remove_const< RangeObj >::type, Fem::SubObject< DomianObj, RangeObj, offset > >
  {
    static void apply ( typename std::remove_const< RangeObj >::type &fm, const Fem::SubObject< DomianObj, RangeObj, offset > &s )
    {
      for( int i = 0; i < s.size(); ++i )
        fm[ i ] = s[ i ];
    }
  };

} //  namespace Dune

#endif // #ifndef DUNE_FEM_SUBOBJECTS_HH

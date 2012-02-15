#ifndef DUNE_FEM_SUBOBJECTS_HH
#define DUNE_FEM_SUBOBJECTS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>


/**************************************/

namespace Dune
{

  template< class DofVector, class Dof > 
  class SubDofVector 
  {
    typedef DofVector DofVectorType;
    typedef Dof DofType;

    public:
      SubDofVector( DofVectorType &dofs, int size, int offset ) :
        dofs_( dofs ),
        offset_ ( offset ),
        size_( size )
      {}

      const DofType  &operator[] ( const int i ) const
      {
        assert( (i < size_ )&& (i >= 0 ) );
        return dofs_[ i + offset_ ];
      }

      DofType &operator[] ( const int i ) 
      {
        assert( (i < size_ )&& (i >= 0 ) );
        return dofs_[ i + offset_ ];
      }

      int size() const
      {
        return size_;
      }


    private:
      DofVectorType &dofs_;
      const int offset_;
      const int size_;
  };


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

  template< class K, int ROWS, int COLS >
  struct RowType< FieldMatrix< K, ROWS, COLS > >
  {
    typedef FieldVector<K, COLS> Type;
    static const int size = ROWS;
  };



  template <class DomainObject, class RangeObject, int offset >
  class SubObject 
  {
    typedef DomainObject DomainObjectType
    typedef RangeObject RangeObjectType;

    typedef typename Dune::RowType< RangeObject > :: Type RowType;

  public:
    SubObject( DomainObjectType &host ):
      host_( host )
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

    int size() const
    {
      return Dune::RowType< RangeObject > :: size;
    }

    operator typename remove_const< RangeObjectType >::type () const
    {
      typename remove_const< RangeObjectType >::type y;
      for( int i = 0; i < size(); ++i )
        y[ i ] = (*this)[ i ];
      return y;
    }

  private:
    DomainObjectType &host_;
  };

} //  end namespace Dune
#endif//#ifndef DUNE_FEM_SUBOBJECTS_HH

#ifndef DUNE_FEM_RECTARRAY_HH
#define DUNE_FEM_RECTARRAY_HH

#include <dune/fem/storage/array.hh>

namespace Dune
{

  template< class RowImp, class RectArrayImp >
  class RectArrayInterface
  {
  public:
    typedef RowImp RowType;
    typedef typename RowType :: ElementType ElementType;

  private:
    typedef RectArrayInterface< RowType, RectArrayImp > ThisType;

  public:
    inline RectArrayImp& operator= ( const ElementType &element )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().operator=( element ) );
      return asImp();
    }
    
    inline const RowType& operator[] ( unsigned int row ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp()[ row ] );
      return asImp()[ row ];
    }

    inline RowType& operator[] ( unsigned int row )
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp()[ row ] );
      return asImp()[ row ];
    }

    inline unsigned int columns () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().columns() );
      return asImp().columns();
    }

    inline unsigned int rows () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().rows() );
      return asImp().rows();
    }

  protected:
    inline const RectArrayImp& asImp () const
    {
      return static_cast< const RectArrayImp& >( *this );
    }

    inline RectArrayImp& asImp ()
    {
      return static_cast< RectArrayImp& >( *this );
    }
  };


  
  template< class RowImp, class RectArrayImp >
  class RectArrayDefault
  : public RectArrayInterface< RowImp, RectArrayImp >
  {
  public:
    typedef RowImp RowType;
    typedef typename RowType :: ElementType ElementType;

  private:
    typedef RectArrayDefault< RowType, RectArrayImp > ThisType;
    typedef RectArrayInterface< RowType, RectArrayImp > BaseType;

  public:
    inline RectArrayImp& operator= ( const ElementType &element )
    {
      RectArrayImp &imp = this->asImp();
      const unsigned int rows = imp.rows();
      for( unsigned int i = 0; i < rows; ++i )
        imp[ i ] = element;
      return imp;
    }
  };



  template< class ElementImp >
  class DynamicRectArray
  : public RectArrayDefault< DynamicArray< ElementImp >,
                             DynamicRectArray< ElementImp > >
  {
  public:
    typedef ElementImp ElementType;
    typedef DynamicArray< ElementType > RowType;

  private:
    typedef DynamicRectArray< ElementType > ThisType;
    typedef RectArrayDefault< RowType, ThisType > BaseType;

  private:
    const unsigned int rows_, columns_;
    RowType *rowArray_;

  public:
    inline DynamicRectArray ( unsigned int rows, unsigned int columns )
    : rows_( rows ),
      columns_( columns )
    {
      unsigned int defaultSize = RowType :: defaultSize( columns_ );
      rowArray_ = new RowType[ rows_ ];
      RowType :: defaultSize( defaultSize );
      assert( rowArray_ != NULL );
    }

    inline ~DynamicRectArray ()
    {
      delete[] rowArray_;
    }
    
    inline ThisType& operator= ( const ElementType &element )
    {
      return BaseType :: operator=( element );
    }
    
    inline const RowType& operator[] ( unsigned int row ) const
    {
      assert( row < rows_ );
      return rowArray_[ row ];
    }

    inline RowType& operator[] ( unsigned int row )
    {
      assert( row < rows_ );
      return rowArray_[ row ];
    }

    inline unsigned int columns ()
    {
      return columns_;
    }

    inline unsigned int rows ()
    {
      return rows_;
    }
  };

}

#endif

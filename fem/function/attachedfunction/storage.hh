#ifndef DUNE_FEM_ATTACHEDFUNCTION_STORAGE_HH
#define DUNE_FEM_ATTACHEDFUNCTION_STORAGE_HH

#include <dune/fem/storage/array.hh>
#include <dune/fem/io/parameter.hh>

namespace Dune
{

  template< class Dof >
  class AttachedDiscreteFunctionStorage
  {
    typedef AttachedDiscreteFunctionStorage< Dof > ThisType;

  public:
    typedef Dof DofType;

  private:
    class DofArray;

    typedef DynamicArray< DofArray > DofStorageType;

  public:
    typedef DofArray value_type;

    typedef typename DofStorageType :: IteratorType IteratorType;
    typedef typename DofStorageType :: ConstIteratorType ConstIteratorType;

  private:
    const unsigned int granularity_;
    unsigned int dofsUsed_;
    DynamicArray< DofArray > dofArrays_;
    DynamicArray< unsigned int > freeDofs_;

  public:
    inline explicit AttachedDiscreteFunctionStorage ( unsigned int size )
    : granularity_( Parameter :: getValidValue
        ( "fem.attacheddiscretefunction.dofgranularity", (unsigned int)4,
          ValidateNotLess< unsigned int >( 1 ) ) ),
      dofsUsed_( 0 ),
      dofArrays_( size ),
      freeDofs_( granularity_ )
    {
      resizeDofArrays( 0, granularity_ );

      for( unsigned int i = 0; i < granularity_; ++i )
        freeDofs_[ i ] = i;
    }

    inline ~AttachedDiscreteFunctionStorage ()
    {
      assert( dofsUsed_ == 0 );
    }

  private:
    // prohibit copying
    AttachedDiscreteFunctionStorage ( const ThisType & );
    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  public:
    inline const DofArray &operator[] ( unsigned int index ) const
    {
      return dofArrays_[ index ];
    }

    inline DofArray &operator[] ( unsigned int index )
    {
      return dofArrays_[ index ];
    }

    inline ConstIteratorType begin () const
    {
      return dofArrays_.begin();
    }

    inline IteratorType begin ()
    {
      return dofArrays_.begin();
    }

    inline ConstIteratorType end () const
    {
      return dofArrays_.end();
    }

    inline IteratorType end ()
    {
      return dofArrays_.end();
    }

    inline unsigned int allocDof ()
    {
      reserveDofs( 1 );
      return freeDofs_[ dofsUsed_++ ];
    }

    inline void freeDof ( unsigned int dof )
    {
      assert( dofsUsed_ > 0 );
      freeDofs_[ --dofsUsed_ ] = dof;
    }

    inline void reserve ( unsigned int newSize )
    {
      dofArrays_.reserve( newSize );
    }

    inline void resize ( unsigned int newSize )
    {
      const unsigned int capacity = freeDofs_.size();
      const unsigned int oldSize = dofArrays_.size();

      dofArrays_.resize( newSize );
      for( unsigned int i = oldSize; i < newSize; ++i )
        dofArrays_[ i ].resize( 0, capacity );
    }

    inline void reserveDofs ( unsigned int additionalSize )
    {
      const unsigned int capacity = freeDofs_.size();
      additionalSize -= (capacity - dofsUsed_);
      if( additionalSize > 0 )
      {
        const unsigned int newCapacity = roundUp( capacity + additionalSize );
        resizeDofArrays( capacity, newCapacity );

        freeDofs_.resize( newCapacity );
        for( unsigned int i = capacity; i < newCapacity; ++i )
          freeDofs_[ i ] = i;
      }
    }

    inline unsigned int size () const
    {
      return dofArrays_.size();
    }

  private:
    inline void resizeDofArrays ( unsigned int oldSize,
                                  unsigned int newSize )
    {
      const unsigned int numDofArrays = dofArrays_.size();
      for( unsigned int i = 0; i < numDofArrays; ++i )
        dofArrays_[ i ].resize( oldSize, newSize );
    }

    inline unsigned int roundUp ( unsigned int size )
    {
      size += granularity_ - 1;
      size -= size % granularity_;
      return size;
    }
  };



  template< class Dof >
  class AttachedDiscreteFunctionStorage< Dof > :: DofArray
  {
    typedef DofArray ThisType;

  public:
    typedef DofType ElementType;

  protected:
    ElementType *dofs_;

  public:
    inline DofArray ()
    : dofs_( 0 )
    {}

    inline explicit DofArray ( unsigned int size )
    : dofs_( new ElementType[ size ] )
    {
      assert( dofs_ != 0 );
    }

  private:
    DofArray ( const ThisType & );

  public:
    inline ~DofArray ()
    {
      if( dofs_ != 0 )
        delete dofs_;
    }

  private:
    ThisType &operator= ( const ThisType & );

  public:
    inline ThisType &operator= ( ThisType &other )
    {
      ElementType *dofs = dofs_;
      dofs_ = other.dofs_;
      other.dofs_ = dofs;
      return *this;
    }

    inline const ElementType &operator[] ( unsigned int index ) const
    {
      return dofs_[ index ];
    }

    inline ElementType &operator[] ( unsigned int index )
    {
      return dofs_[ index ];
    }

    inline void resize ( unsigned int oldSize, unsigned int newSize )
    {
      ElementType *newdofs = new ElementType[ newSize ];
      if( dofs_ != 0 )
      {
        const unsigned int copySize = std :: min( oldSize, newSize );
        for( unsigned int i = 0; i < copySize; ++i )
          newdofs[ i ] = dofs_[ i ];
        delete dofs_;
      }
      dofs_ = newdofs;
    }
  };

}

#endif

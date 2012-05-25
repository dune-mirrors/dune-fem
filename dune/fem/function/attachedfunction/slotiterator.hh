#ifndef DUNE_FEM_SLOTITERATOR_HH
#define DUNE_FEM_SLOTITERATOR_HH

namespace Dune
{

  namespace
  {

    template< class Container >
    struct ExtractElement
    {
      typedef typename Container :: ElementType Type;
    };

    template< class Container >
    struct ExtractElement< const Container >
    {
      typedef const typename Container :: ElementType Type;
    };

  }



  template< class StorageIterator >
  class SlotIterator
  {
    typedef SlotIterator< StorageIterator > ThisType;

  public:
    typedef StorageIterator StorageIteratorType;
    typedef typename StorageIteratorType :: ElementType ElementArrayType;
    typedef typename ExtractElement< ElementArrayType > :: Type ElementType;

  protected:
    StorageIteratorType storageIterator_;
    unsigned int slot_;

  public:
    inline SlotIterator ( const StorageIteratorType &storageIterator,
                          unsigned int slot )
    : storageIterator_( storageIterator ),
      slot_( slot )
    {}

    inline SlotIterator ( const ThisType &other )
    : storageIterator_( other.storageIterator_ ),
      slot_( other.slot_ )
    {}

    inline ThisType &operator= ( const ThisType &other )
    {
      storageIterator_ = other.storageIterator_;
      slot_ = other.slot_;
    }

    inline ThisType &operator++ ()
    {
      ++storageIterator_;
      return *this;
    }

    inline ElementType &operator* ()
    {
      return (*storageIterator_)[ slot_ ];
    }

    inline ElementType *operator-> ()
    {
      return &(operator*());
    }

    inline bool operator== ( const ThisType &other ) const
    {
      return
        (storageIterator_ == other.storageIterator_) && (slot_ == other.slot_);
    }

    inline bool operator!= ( const ThisType &other ) const
    {
      return !(*this == other);
    }
  };

}

#endif

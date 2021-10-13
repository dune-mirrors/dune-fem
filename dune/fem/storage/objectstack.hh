#ifndef DUNE_FEM_OBJECTSTACK_HH
#define DUNE_FEM_OBJECTSTACK_HH

#include <dune/fem/storage/referencecounter.hh>
#include <dune/fem/misc/threads/threadsafevalue.hh>

namespace Dune
{

  namespace Fem
  {

    template< class ObjectFactoryImp >
    class ObjectStack;

    template< class ObjectFactoryImp >
    class ObjectStackEntry;


    template< class ObjectFactoryImp >
    struct ObjectStackEntryTraits
    {
      typedef ObjectStackEntry< ObjectFactoryImp > ReferenceCounterType;

      typedef typename ObjectFactoryImp :: ObjectType ObjectType;
    };


    template< class ObjectFactoryImp >
    class ObjectStackEntry
    : public ReferenceCounterDefault< ObjectStackEntryTraits< ObjectFactoryImp > >
    {
    public:
      typedef ObjectFactoryImp ObjectFactoryType;

      typedef ObjectStackEntryTraits< ObjectFactoryType > Traits;

    private:
      typedef ObjectStackEntry< ObjectFactoryType > ThisType;
      typedef ReferenceCounterDefault< Traits > BaseType;

      friend class ObjectStack< ObjectFactoryType >;

    protected:
      typedef ObjectStack< ObjectFactoryType > ObjectStackType;

    public:
      typedef typename ObjectFactoryType :: ObjectType ObjectType;

    protected:
      // reference to the stack
      ObjectStackType &stack_;

      // pointer to the actual object
      ObjectType *const object_;

      // next object on the stack
      ThisType *next_;

      inline explicit ObjectStackEntry ( ObjectStackType &stack )
      : BaseType( 0 ),
        stack_( stack ),
        object_( stack_.factory().newObject() )
      {
      }

    public:
      ObjectStackEntry ( const ThisType& ) = delete;

      inline ~ObjectStackEntry ()
      {
        delete object_;
      }

      ThisType& operator= ( const ThisType& ) = delete;

      inline operator const ObjectType& () const
      {
        return *object_;
      }

      inline operator ObjectType& ()
      {
        return *object_;
      }

      inline void deleteObject ()
      {
        stack_.push( this );
      }

      inline const ObjectType &getObject () const
      {
        return *object_;
      }

      inline ObjectType &getObject ()
      {
        return *object_;
      }
    };



    //! \ingroup HelperClasses
    //! Stores pointers to a given class in a stack
    //! used for local functions and for basefunctionsets
    template< class ObjectFactoryImp >
    class ObjectStack
    {
      typedef ObjectFactoryImp ObjectFactoryType;

    private:
      typedef ObjectStack< ObjectFactoryType > ThisType;

      friend class ObjectStackEntry< ObjectFactoryType >;

    public:
      //! type of the stored objects
      typedef typename ObjectFactoryType :: ObjectType ObjectType;

      //! type of the storage objects
      typedef ObjectStackEntry< ObjectFactoryType > StackEntryType;

      //! type of object pointers
      typedef ObjectPointer< StackEntryType > PointerType;

    protected:
      const ObjectFactoryType &factory_;

      // store stack such that entries are thread safe
      typedef StackEntryType*  StackEntryPtrType;
      typedef ThreadSafeValue< StackEntryPtrType > ThreadSafeValuesType;

      // thread safe stack entries (in multi thread mode a vector)
      ThreadSafeValuesType stackEntries_;
    public:
      //! constructor
      ObjectStack ( const ObjectFactoryType &factory )
      : factory_( factory ),
        stackEntries_( StackEntryPtrType(0) )
      {
      }

      ObjectStack ( const ThisType& ) = delete;

      //! delete all objects on stack
      ~ObjectStack ()
      {
        // make sure this is only called in single thread mode
        // because the master thread is taking care of all object delete
        assert( MPIManager::singleThreadMode() );
        const size_t threadSize = stackEntries_.size();
        for( size_t i=0; i<threadSize; ++i )
        {
          StackEntryPtrType& stackEntry = stackEntries_[ i ];
          while ( stackEntry != 0 )
          {
            StackEntryType *obj = stackEntry;
            stackEntry = obj->next_;
            delete obj;
          }
        }
      }

      ThisType& operator= ( const ThisType& ) = delete;

      //! get an object pointer to a storage object
      inline PointerType getObject ()
      {
        return PointerType( pop() );
      }

    protected:
      inline const ObjectFactoryType &factory() const
      {
        return factory_;
      }

      // push storage object to the stack
      inline void push ( StackEntryType *obj )
      {
        // get thread private value
        push( obj, *stackEntries_ );
      }

      // pop a storage object from the stack
      inline StackEntryType *pop ()
      {
        // get thread private value
        return pop( *stackEntries_ );
      }

    private:
      // push storage object to the stack
      inline void push ( StackEntryType *obj, StackEntryPtrType& stackEntry )
      {
        obj->next_ = stackEntry;
        stackEntry = obj;
      }

      // pop a storage object from the stack
      inline StackEntryType *pop ( StackEntryPtrType& stackEntry )
      {
        StackEntryType *ptr = stackEntry;
        if( ptr != 0 )
          stackEntry = stackEntry->next_;
        else
          ptr = new StackEntryType( *this );
        return ptr;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OBJECTSTACK_HH

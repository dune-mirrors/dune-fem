#ifndef DUNE_FEM_OBJECTSTACK_HH
#define DUNE_FEM_OBJECTSTACK_HH

#include <dune/fem/misc/debug.hh>
#include <dune/fem/storage/referencecounter.hh>
#include <dune/fem/misc/threads/threadsafevalues.hh>

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

      template< class, class >
      friend class Conversion;

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

    protected:
      inline explicit ObjectStackEntry ( ObjectStackType &stack )
      : BaseType( 0 ),
        stack_( stack ),
        object_( stack_.factory().newObject() )
      {
      }

    private:
      // prohibit copying
      ObjectStackEntry ( const ThisType & );

    public:
      inline ~ObjectStackEntry ()
      {
        delete object_;
      }
      
    private:
      // prohivit assignment
      ThisType &operator= ( const ThisType & );

    public:
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

      typedef std::pair< StackEntryType*, DebugCounter<> > StackEntryPairType ;
      typedef ThreadSafeValue< StackEntryPairType > ThreadSafeValuesType;

      // thread safe stack entries (in multi thread mode a vector)
      ThreadSafeValuesType stackEntries_;
    public:
      //! constructor 
      ObjectStack ( const ObjectFactoryType &factory )
      : factory_( factory ),
        stackEntries_( StackEntryPairType( (StackEntryType* )0, DebugCounter<> () ) )
      {
      }

    private:
      // Disallow copying
      ObjectStack ( const ThisType & );

    public:
      //! delete all objects on stack 
      ~ObjectStack ()
      {
        // make sure this is only called in single thread mode 
        // because the master thread is taking care of all object delete
        assert( ThreadManager::singleThreadMode() );
        const size_t threadSize = stackEntries_.size();
        for( size_t i=0; i<threadSize; ++i )
        {
          StackEntryPairType& stackEntry = stackEntries_[ i ];
#ifndef BASEFUNCTIONSET_CODEGEN_GENERATE
          // this assertion will fail during code generation 
          assert( stackEntry.second == 0 );
#endif
          while ( stackEntry.first != 0 )
          {
            StackEntryType *obj = stackEntry.first;
            stackEntry.first = obj->next_;
            delete obj;
          }
        }
      }
      
    private:
      // Disallow copying
      ThisType &operator= ( const ThisType & );

    public:
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
        push( obj, *stackEntries_ );
      }

      // pop a storage object from the stack
      inline StackEntryType *pop ()
      {
        return pop( *stackEntries_ );
      }

    private:  
      // push storage object to the stack
      inline void push ( StackEntryType *obj, StackEntryPairType& stackEntry )
      {
        --stackEntry.second;
        obj->next_ = stackEntry.first;
        stackEntry.first = obj;
      }

      // pop a storage object from the stack
      inline StackEntryType *pop ( StackEntryPairType& stackEntry )
      {
        ++stackEntry.second;

        StackEntryType *ptr = stackEntry.first;
        if( ptr != 0 )
          stackEntry.first = stackEntry.first->next_;
        else
          ptr = new StackEntryType( *this );
        return ptr;
      }
    };

  } // namespace Fem
   
} // namespace Dune

#endif // #ifndef DUNE_FEM_OBJECTSTACK_HH

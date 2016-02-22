#ifndef DUNE_FEM_REFERENCECOUNTER_HH
#define DUNE_FEM_REFERENCECOUNTER_HH

#include <cassert>
#include <type_traits>

#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class ReferenceCounterInterface
     *  \brief interface for objects capable of reference counting
     *
     *  In many situations it has to be decided, whether an object is still
     *  required or must be freed. A simple approach to this decision is to
     *  count all references to an object. When the last reference is released,
     *  the object is freed.
     *
     *  The ReferenceCounterInterface provides methods, an object needs to
     *  support reference counting. Classes implementing this interface can be
     *  used with ObjectReference and ObjectPointer. Both will implicitly call
     *  addReference when they start pointing to an object and removeReference
     *  when they stop pointing to it.
     *
     *  \note Reference counting is not a good mehtod for garbage collection,
     *        since it does not detect cycles. So be careful when using reference
     *        counters!
     */
    template< class RCT >
    class ReferenceCounterInterface
    : public BartonNackmanInterface< ReferenceCounterInterface< RCT >, typename RCT::ReferenceCounterType >
    {
      typedef ReferenceCounterInterface< RCT > ThisType;
      typedef BartonNackmanInterface< ThisType, typename RCT::ReferenceCounterType > BaseType;

    public:
      //! type of the traits
      typedef RCT Traits;

      //! type of the implementation (Barton-Nackman)
      typedef typename Traits::ReferenceCounterType ReferenceCounterType;

      //! type of the reference counter interface
      typedef ThisType ReferenceCounterInterfaceType;

      //! type of the object, this is a reference counter for
      typedef typename Traits::ObjectType ObjectType;

    protected:
      using BaseType::asImp;

    public:
      /** \brief add a reference to this object
       *
       *  This method should be called whenever a permanent reference to this
       *  object is established (a pointer counts as a reference, too)
       *
       *  \note This method is declared const, since we want to be able to add
       *        references to const objects, too. Hence, the reference coutner
       *        must be declared mutable.
       */
      void addReference () const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().addReference() );
      }

      /** \brief delete to object
       *
       *  This method is used by the default implementation to free the object
       *  when the reference counter becomes zero.
       *
       *  \note For the implementor's convenience, this method is not declared
       *        const. The caller (usually the default implementation's
       *        removeReference method) is responsible for performing a necessary
       *        const_cast.
       */
      void deleteObject ()
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().deleteObject() );
      }

      /** \brief access the real object (const version)
       *
       *  Sometimes the reference counter just wraps the object for which it
       *  counts references. In such cases, this method returns the wrapped
       *  object. Otherwise, the object itself may be returned.
       *
       *  \returns a constant reference to the real object
       */
      const ObjectType &getObject () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().getObject() );
        return asImp().getObject();
      }

      /** \brief access the real object (non-const version)
       *
       *  Sometimes the reference counter just wraps the object for which it
       *  counts references. In such cases, this method returns the wrapped
       *  object. Otherwise, the object itself may be returned.
       *
       *  \returns a reference to the real object
       */
      ObjectType &getObject ()
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().getObject() );
        return asImp().getObject();
      }

      /** \brief remove a reference to this object
       *
       *  This method should be called whenever a previously added reference
       *  to this object is released.
       *
       *  \note This method is declared const, since we want to be able to remove
       *        references to const objects, too. Hence, the reference coutner
       *        must be declared mutable.
       */
      void removeReference () const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().removeReference() );
      }
    };



    template< class ReferenceCounter >
    struct SupportsReferenceCounterInterface
    {
      typedef ReferenceCounterInterface< typename ReferenceCounter::Traits > ReferenceCounterInterfaceType;
      static const bool v = std::is_convertible< ReferenceCounter, ReferenceCounterInterfaceType >::value;
    };



    /** \class ReferenceCounterDefault
     *  \brief default implementation of ReferenceCounterInterface
     *
     *  Reference counting usually uses a class variable to hold the number of
     *  references to the instance. To make reference counting simple to use,
     *  this mechanism is default implemented.
     *
     *  \note The deleteObject method is used to delete the object when the
     *        reference counter reaches zero. To perform any action other than
     *        deleting the object, simply overwrite this method.
     */
    template< class RCT >
    class ReferenceCounterDefault
    : public ReferenceCounterInterface< RCT >
    {
      typedef ReferenceCounterDefault< RCT > ThisType;
      typedef ReferenceCounterInterface< RCT > BaseType;

    public:
      //! type of the implementation (Barton-Nackman)
      typedef typename BaseType::ReferenceCounterType ReferenceCounterType;

    protected:
      using BaseType::asImp;

    public:
      /** \brief constructor initializing the reference counter
       *
       *  \note Since we want to be able to count referenced to stack objects,
       *        we initialize the reference counter to 1 by default. This means
       *        that after calling new we already have one reference to the
       *        object. In such cases, just call removeReference directly when
       *        the reference is no longer used (instead of simply deleting it).
       *
       *  \param[in]  refcount  start value for the reference counter; this
       *              value defalts to 1
       */
      explicit ReferenceCounterDefault ( unsigned int refcount = 1 )
      : refcount_( refcount )
      {}

      ReferenceCounterDefault ( const ThisType& ) = delete;
      ThisType& operator= ( const ThisType& ) = delete;

      /** \copydoc Dune :: ReferenceCounterInterface :: addReference */
      void addReference () const
      {
        ++refcount_;
      }

      /** \copydoc Dune :: ReferenceCounterInterface :: deleteObject */
      void deleteObject ()
      {
        delete this;
      }

      /** \copydoc Dune :: ReferenceCounterInterface :: removeReference */
      void removeReference () const
      {
        assert( refcount_ > 0 );
        --refcount_;
        if( refcount_ == 0 )
          const_cast< ReferenceCounterType& >( asImp() ).deleteObject();
      }

      /** \brief return current reference count */
      unsigned int referenceCounter () const { return refcount_; }

    protected:
      mutable unsigned int refcount_;
    };



    /** \class ObjectPointer
     *  \brief models a pointer to a reference countable object
     *
     *  ObjectPointer tries to behave like a normal pointer to an object
     *  implementing the ReferenceCounterInterface. Internally, however, it
     *  calls the object's addReference and removeReference methods whenever the
     *  pointer is created, assign or deleted.
     */
    template< class ReferenceCounter >
    class ObjectPointer
    {
      typedef ObjectPointer< ReferenceCounter > ThisType;

      static_assert( SupportsReferenceCounterInterface< ReferenceCounter >::v, "ObjectPointer can only point to reference counting types." );

    public:
      //! type of the object, this pointer points to
      typedef ReferenceCounter ReferenceCounterType;

      typedef typename ReferenceCounterType::ObjectType ObjectType;

      /** \brief initialize a pointer (with a standard C++ pointer)
       *
       *  \param[in]  object  C++ pointer to initialize this pointer with; the
       *                      default value is 0
       */
      explicit ObjectPointer ( ReferenceCounterType *const object = 0 )
      : object_( object )
      {
        if( object_ != 0 )
          object_->addReference();
      }

      /** \brief copy constructor
       *
       *  Copying an ObjectPointer will also increase the reference counter of
       *  the object pointed to.
       *
       *  \param[in]  other  pointer to assign to this one
       */
      ObjectPointer ( const ThisType &other )
      : object_( other.object_ )
      {
        if( object_ != 0 )
          object_->addReference();
      }

      /** \brief destructor
       *
       *  When the pointer is deleted, the reference counter of the object pointed
       *  to is automatically decreased.
       */
      ~ObjectPointer ()
      {
        if( object_ != 0 )
          object_->removeReference();
      }

      /** \brief assign another pointer to this one.
       */
      ThisType &operator= ( const ThisType &other )
      {
        // Note that it is safe to remove the reference first. If other holds
        // a reference to the same object, the reference counter cannot reach
        // zero.
        if( object_ != 0 )
          object_->removeReference();
        object_ = other.object_;
        if( object_ != 0 )
          object_->addReference();

        return *this;
      }

      /** \brief dereference the ObjectPointer
       *
       *  \note This method asserts that the pointer is not 0.
       *
       *  \returns a reference to the object pointed to.
       */
      ObjectType &operator* () const
      {
        assert( object_ != 0 );
        return object_->getObject();
      }

      /** \brief return current reference count */
      unsigned int referenceCounter () const { return object_->referenceCounter(); }

    protected:
      ReferenceCounterType *object_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_REFERENCECOUNTER_HH

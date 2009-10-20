#ifndef DUNE_FEM_DEBUG_HH
#define DUNE_FEM_DEBUG_HH

#include <cassert>

namespace Dune
{

  /** \class DebugCounter
   *  \brief A counter only present if NDEBUG is not defined
   *
   *  There are several cases, where we need a counter for debugging
   *  purposes that should only be present, if NDEBUG is not defined. 
   *
   *  In debug mode, this counter wraps a standard integer type, 
   *  otherwise its size is zero.
   *
   *  \note The comparison operators always return true, 
   *  if NDEBUG is defined!
   */
  template< class CounterImp = unsigned int >
  class DebugCounter
  {
  public:
    //! integral type for the actual counting
    typedef CounterImp CounterType;
    
  private:
    typedef DebugCounter< CounterType > ThisType;

  protected:
#ifndef NDEBUG
    CounterType count_;
#endif
    
  public:
    /** \brief constructor
     *
     *  \note This constructor implicitly defines a conversion from CounterType
     *        to DebugCounter. This is very useful in comparison statements.
     *
     *  \param[in]  count  value to initialize the counter with (defaults to 0)
     */
    inline DebugCounter ( const CounterType count = 0 )
#ifndef NDEBUG
    : count_( count )
#endif
    {
    }

    /** \brief copy constructor
     */
    inline DebugCounter ( const ThisType &other )
#ifndef NDEBUG
    : count_( other.count_ )
#endif
    {
    }

    /** \brief increment operator
     *
     *  If NDEBUG is not defined, the counter is incremented by 1. Otherwise
     *  nothing happens (and the entire call will be removed during
     *  oprimization).
     */
    inline ThisType &operator++ ()
    {
#ifndef NDEBUG
      ++count_;
#endif
      return *this;
    }

    /** \brief decrement operator
     *
     *  If NDEBUG is not defined, the counter is decremented by 1. Otherwise
     *  nothing happens (and the entire call will be removed during
     *  oprimization).
     */
    inline ThisType &operator-- ()
    {
#ifndef NDEBUG
      --count_;
#endif
      return *this;
    }

    /** \brief comparison for equality
     *
     *  Compares to DebugCounters for equality. If NDEBUG is defined, the
     *  result will be true.
     *
     *  \note Due to the implicit conversion, the second argument may also be
     *        of CounterType.
     *
     *  \param[in]  other  DebugCounter to compare this one to
     *
     *  \returns true, if the counters equal or NDEBUG is defined
     */
    inline bool operator== ( const ThisType &other )
    {
#ifndef NDEBUG
      return count_ == other.count_;
#else
      return true;
#endif
    }

    /** \brief comparison for inequality
     *
     *  Compares to DebugCounters for inequality. If NDEBUG is defined, the
     *  result will be true.
     *
     *  \note Due to the implicit conversion, the second argument may also be
     *        of CounterType.
     *
     *  \param[in]  other  DebugCounter to compare this one to
     *
     *  \returns true, if the counters differ or NDEBUG is defined
     */
    inline bool operator!= ( const ThisType &other )
    {
#ifndef NDEBUG
      return count_ != other.count_;
#else
      return true;
#endif
    }
  };



  class DebugLock
  {
  private:
    typedef DebugLock ThisType;
    
  protected:
#ifndef NDEBUG
    bool lock_;
#endif

  public:
    inline DebugLock ()
#ifndef NDEBUG
    : lock_( false )
#endif
    {
    }
    
  private:
    // prohibit copying
    DebugLock ( const ThisType & );

    // prohibit copying
    ThisType &operator= ( const ThisType & );

  public:
    inline bool operator ! () const
    {
#ifndef NDEBUG
      return !lock_;
#else
      return true;
#endif
    }
    
    inline void lock ()
    {
#ifndef NDEBUG
      assert( !lock_ );
      lock_ = true;
#endif
    }

    inline void unlock ()
    {
#ifndef NDEBUG
      assert( lock_ );
      lock_ = false;
#endif
    }
  };

};

#endif

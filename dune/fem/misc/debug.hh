#ifndef DUNE_FEM_DEBUG_HH
#define DUNE_FEM_DEBUG_HH

#include <cassert>

namespace Dune
{

  namespace Fem
  {

#if not defined NDEBUG
#define USE_DEBUG_CNT
#endif

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
#ifdef USE_DEBUG_CNT
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
#ifdef USE_DEBUG_CNT
      : count_( count )
#endif
      {
      }

      /** \brief copy constructor
       */
      inline DebugCounter ( const ThisType &other )
#ifdef USE_DEBUG_CNT
      : count_( other.count_ )
#endif
      {
      }

      /** \brief increment operator
       *
       *  If USE_DEBUG_CNT is not defined, the counter is incremented by 1. Otherwise
       *  nothing happens (and the entire call will be removed during
       *  oprimization).
       */
      inline ThisType &operator++ ()
      {
#ifdef USE_DEBUG_CNT
        ++count_;
#endif
        return *this;
      }

      /** \brief decrement operator
       *
       *  If USE_DEBUG_CNT is not defined, the counter is decremented by 1. Otherwise
       *  nothing happens (and the entire call will be removed during
       *  oprimization).
       */
      inline ThisType &operator-- ()
      {
#ifdef USE_DEBUG_CNT
        --count_;
#endif
        return *this;
      }

      /** \brief comparison for equality
       *
       *  Compares to DebugCounters for equality. If USE_DEBUG_CNT is defined, the
       *  result will be true.
       *
       *  \note Due to the implicit conversion, the second argument may also be
       *        of CounterType.
       *
       *  \param[in]  other  DebugCounter to compare this one to
       *
       *  \returns true, if the counters equal or USE_DEBUG_CNT is defined
       */
      inline bool operator== ( const ThisType &other )
      {
#ifdef USE_DEBUG_CNT
        return count_ == other.count_;
#else
        return true;
#endif
      }

      /** \brief comparison for inequality
       *
       *  Compares to DebugCounters for inequality. If USE_DEBUG_CNT is defined, the
       *  result will be true.
       *
       *  \note Due to the implicit conversion, the second argument may also be
       *        of CounterType.
       *
       *  \param[in]  other  DebugCounter to compare this one to
       *
       *  \returns true, if the counters differ or USE_DEBUG_CNT is defined
       */
      inline bool operator!= ( const ThisType &other )
      {
#ifdef USE_DEBUG_CNT
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
#ifdef USE_DEBUG_CNT
      bool lock_;
#endif

    public:
      inline DebugLock ()
#ifdef USE_DEBUG_CNT
      : lock_( false )
#endif
      {
      }

      DebugLock ( const ThisType& ) = delete;
      ThisType& operator= ( const ThisType& ) = delete;

      inline bool operator ! () const
      {
#ifdef USE_DEBUG_CNT
        return !lock_;
#else
        return true;
#endif
      }

      inline void lock ()
      {
#ifdef USE_DEBUG_CNT
        assert( !lock_ );
        lock_ = true;
#endif
      }

      inline void unlock ()
      {
#ifdef USE_DEBUG_CNT
        assert( lock_ );
        lock_ = false;
#endif
      }
    };

  } // namespace Fem

} // namespace Dune

#endif

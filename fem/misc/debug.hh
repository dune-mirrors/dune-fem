#ifndef DUNE_FEM_DEBUG_HH
#define DUNE_FEM_DEBUG_HH

namespace Dune
{

  template< class CounterImp = unsigned int >
  class DebugCounter
  {
  public:
    typedef CounterImp CounterType;
    
  private:
    typedef DebugCounter< CounterType > ThisType;

  protected:
#ifndef NDEBUG
    CounterType count_;
#endif
    
  public:
    inline DebugCounter ( const CounterType count = 0 )
#ifndef NDEBUG
    : count_( count )
#endif
    {
    }

    inline DebugCounter ( const ThisType &other )
#ifndef NDEBUG
    : count_( other.count_ )
#endif
    {
    }

    inline ThisType &operator++ ()
    {
#ifndef NDEBUG
      ++count_;
#endif
      return *this;
    }

    inline ThisType &operator-- ()
    {
#ifndef NDEBUG
      --count_;
#endif
      return *this;
    }

    inline bool operator== ( const ThisType &other )
    {
#ifndef NDEBUG
      return count_ == other.count_;
#else
      return true;
#endif
    }

    inline bool operator!= ( const ThisType &other )
    {
#ifndef NDEBUG
      return count_ != other.count_;
#else
      return true;
#endif
    }
  };

  
};

#endif

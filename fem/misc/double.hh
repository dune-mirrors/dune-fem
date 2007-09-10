#ifndef DUNE_FEM_DOUBLE_HH
#define DUNE_FEM_DOUBLE_HH

#include <iostream>

#include <dune/fem/io/file/xdrio.hh>

namespace Dune
{

#ifdef COUNT_FLOPS

  template< class FloatImp >
  class FlOpCounter;



  template<>
  class FlOpCounter< void >
  {
  private:
    typedef FlOpCounter< void > ThisType;

  protected:
    unsigned long count_;

  protected:
    inline FlOpCounter ()
    : count_( 0 )
    {
    }

  public:
    inline ~FlOpCounter ()
    {
      std :: cout << "Total number of floating point operations: "
                  << count_ << std :: endl;
    }


    inline ThisType &operator++ ()
    {
      ++count_;
      return *this;
    }

    inline static ThisType &instance ()
    {
      static ThisType instance;
      return instance;
    }
  };



  template< class FloatImp >
  class FlOpCounter
  {
  public:
    typedef FloatImp FloatType;

  private:
    typedef FlOpCounter< FloatType > ThisType;

  protected:
    unsigned long count_;

  protected:
    inline FlOpCounter ()
    : count_( 0 )
    {
    }

  public:
    inline ~FlOpCounter ()
    {
      std :: cout << "Number of floating point operations for "
                  << FloatType :: typeName() << ": "
                  << count_ << std :: endl;
    }

    inline ThisType &operator++ ()
    {
      ++count_;
      ++(FlOpCounter< void > :: instance());
      return *this;
    }

    inline static ThisType &instance ()
    {
      static ThisType instance;
      return instance;
    }
  };
  
#else

  template< class FloatImp >
  class FlOpCounter
  {
  public:
    typedef FloatImp FloatType;

  private:
    typedef FlOpCounter< FloatType > ThisType;

  protected:
    inline FlOpCounter ()
    {
    }

  public:
    inline ThisType &operator++ ()
    {
      return *this;
    }

    inline static ThisType &instance ()
    {
      static ThisType instance;
      return instance;
    }
  };

#endif

  //- forward declaration 
  class Double; 
  // wrap of std power 
  static double pow (const Double& v, const double p);
  // wrap of std log 
  static double log (const Double& v);
  // wrap of std sqrt 
  static double sqrt(const Double& v);

  class Double
  {
  private:
    typedef Double ThisType;

    friend Double operator+ ( const Double&, const Double& );
    friend Double operator+ ( const Double&, const double );
    friend Double operator+ ( const double, const Double& );
    friend Double operator+ ( const Double&, const int );
    friend Double operator+ ( const int, const Double& );
    friend Double operator+ ( const Double&, const unsigned int );
    friend Double operator+ ( const unsigned int, const Double& );

    friend Double operator- ( const Double&, const Double& );
    friend Double operator- ( const Double&, const double );
    friend Double operator- ( const double, const Double& );
    friend Double operator- ( const Double&, const int );
    friend Double operator- ( const int, const Double& );
    friend Double operator- ( const Double&, const unsigned int );
    friend Double operator- ( const unsigned int, const Double& );
    
    friend Double operator* ( const Double&, const Double& );
    friend Double operator* ( const Double&, const double );
    friend Double operator* ( const double, const Double& );
    friend Double operator* ( const Double&, const int );
    friend Double operator* ( const int, const Double& );
    friend Double operator* ( const Double&, const unsigned int );
    friend Double operator* ( const unsigned int, const Double& );

    friend Double operator/ ( const Double&, const Double& );
    friend Double operator/ ( const double, const Double& );
    friend Double operator/ ( const Double&, const double );
    friend Double operator/ ( const int, const Double& );
    friend Double operator/ ( const Double&, const int );
    friend Double operator/ ( const unsigned int, const Double& );
    friend Double operator/ ( const Double&, const unsigned int );

    friend bool operator== ( const Double&, const Double& );
    friend bool operator== ( const double, const Double& );
    friend bool operator== ( const Double&, const double );
    friend bool operator== ( const int, const Double& );
    friend bool operator== ( const Double&, const int );
    friend bool operator== ( const unsigned int, const Double& );
    friend bool operator== ( const Double&, const unsigned int );
    
    friend bool operator!= ( const Double&, const Double& );
    friend bool operator!= ( const double, const Double& );
    friend bool operator!= ( const Double&, const double );
    friend bool operator!= ( const int, const Double& );
    friend bool operator!= ( const Double&, const int );
    friend bool operator!= ( const unsigned int, const Double& );
    friend bool operator!= ( const Double&, const unsigned int );
    
    friend bool operator< ( const Double&, const Double& );
    friend bool operator< ( const double, const Double& );
    friend bool operator< ( const Double&, const double );
    friend bool operator< ( const int, const Double& );
    friend bool operator< ( const Double&, const int );
    friend bool operator< ( const unsigned int, const Double& );
    friend bool operator< ( const Double&, const unsigned int );

    friend bool operator<= ( const Double&, const Double& );
    friend bool operator<= ( const double, const Double& );
    friend bool operator<= ( const Double&, const double );
    friend bool operator<= ( const int, const Double& );
    friend bool operator<= ( const Double&, const int );
    friend bool operator<= ( const unsigned int, const Double& );
    friend bool operator<= ( const Double&, const unsigned int );
    
    friend bool operator> ( const Double&, const Double& );
    friend bool operator> ( const double, const Double& );
    friend bool operator> ( const Double&, const double );
    friend bool operator> ( const int, const Double& );
    friend bool operator> ( const Double&, const int );
    friend bool operator> ( const unsigned int, const Double& );
    friend bool operator> ( const Double&, const unsigned int );

    friend bool operator>= ( const Double&, const Double& );
    friend bool operator>= ( const double, const Double& );
    friend bool operator>= ( const Double&, const double );
    friend bool operator>= ( const int, const Double& );
    friend bool operator>= ( const Double&, const int );
    friend bool operator>= ( const unsigned int, const Double& );
    friend bool operator>= ( const Double&, const unsigned int );

    friend std :: ostream &operator<< ( std :: ostream&, const Double& );
    friend std :: istream &operator>> ( std :: istream&, Double& );

    friend double pow (const Double& v, const double p);
    friend double log (const Double& v);
    friend double sqrt(const Double& v);


    friend struct XdrIO< Double >;

  protected:
    typedef FlOpCounter< ThisType > FlOpCounterType;

  protected:
    double value_;

  public:
    inline Double ()
    {
    }

    inline Double ( const double value )
    : value_( value )
    {
    }

    inline Double ( const ThisType &other )
    : value_( other.value_ )
    {
    }

    inline operator int () const
    {
      return (int)value_;
    }

    inline ThisType &operator= ( const ThisType other )
    {
      value_ = other.value_;
      return *this;
    }
    
    inline ThisType &operator+= ( const ThisType other )
    {
      flOp();
      value_ += other.value_;
      return *this;
    }
    
    inline ThisType &operator-= ( const ThisType other )
    {
      flOp();
      value_ -= other.value_;
      return *this;
    }
    
    inline ThisType &operator*= ( const ThisType other )
    {
      flOp();
      value_ *= other.value_;
      return *this;
    }
    
    inline ThisType &operator/= ( const ThisType other )
    {
      flOp();
      value_ /= other.value_;
      return *this;
    }

    static std :: string typeName ()
    {
      return "Double";
    }

  protected:
    static inline void flOp ()
    {
      ++(FlOpCounterType :: instance());
    }
  };

  static inline double pow (const Double& v, const double p)
  {
    return std::pow(v.value_,p);
  }

  static inline double log (const Double& v)
  {
    return std::log(v.value_);
  }

  static inline double sqrt(const Double& v)
  {
    return std::sqrt(v.value_);
  }

  template<>
  struct XdrIO< Double >
  {
    static inline int io( XDR *xdrs, Double &a )
    {
      return XdrIO< double > :: io( xdrs, a.value_ );
    }
  };

}

#include "double_inline.hh"

#endif

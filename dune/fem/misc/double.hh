#ifndef DUNE_FEM_DOUBLE_HH
#define DUNE_FEM_DOUBLE_HH

//- system includes
#include <iostream>
#include <cmath>
#include <limits>
#include <string>

//- Dune includes
#include <dune/common/version.hh>
#include <dune/common/typetraits.hh>

#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/threads/threadsafevalue.hh>
#include <dune/fem/storage/singleton.hh>

namespace Dune
{

  namespace Fem
  {

#ifdef COUNT_FLOPS

    template< class FloatImp >
    class FlOpCounter;


    template <class FloatImp>
    class FlOpSummary
    {
      typedef FlOpSummary< FloatImp >   ThisType;

      std::string name_;
      ThreadSafeValue< unsigned long > count_;

      FlOpSummary(const std::string& name) : name_(name), count_( 0 ) {}
    public:
      ~FlOpSummary ()
      {
        unsigned long totalCount = 0;
        for( size_t i=0; i<count_.size(); ++i )
        {
          std::cout << name_ << " thread[ " << i << " ]: " << count_[ i ] << std::endl;
          totalCount += count_[ i ];
        }

        std :: cout << "Total number of floating point operations (" << name_ << "): "
                    << totalCount << std :: endl;
      }

      inline void add( const unsigned long count, const int thread )
      {
        count_[ thread ] += count ;
      }

      friend class Dune::Fem::Singleton< ThisType >;

      static ThisType& instance( const std::string& name )
      {
        return Singleton< ThisType >::instance( name );
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
      const int thread_;

    protected:
      inline FlOpCounter ()
      : count_( 0 ), thread_( MPIManager::thread() )
      {
        FlOpSummary< FloatImp > :: instance( FloatImp::typeName() ) ;
      }

    public:
      inline ~FlOpCounter ()
      {
        FlOpSummary< FloatImp > :: instance( FloatImp::typeName() ).add( count_, thread_ );
      }

      inline ThisType &operator++ ()
      {
        ++(count_);
        return *this;
      }

      inline ThisType &operator += ( const unsigned long i )
      {
        count_ += i;
        return *this;
      }

      friend class Dune::Fem::Singleton< ThisType >;
      static ThisType &instance ()
      {
#ifdef HAVE_PTHREAD
        static thread_local ThisType instance;
        return instance;
#else
        return Singleton< ThisType >::instance();
#endif
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

      inline ThisType &operator += ( const unsigned long i )
      {
        return *this;
      }

      DUNE_EXPORT static ThisType &instance ()
      {
#ifdef HAVE_PTHREAD
        static thread_local ThisType instance;
#else
        static ThisType instance;
#endif
        return instance;
      }
    };

#endif

    //- forward declaration
    class Double;
    // wrap of std log
    static double log (const Double& v);
    // wrap of std sqrt
    static double sqrt(const Double& v);
    // wrap of std sin
    static double cos (const Double& v);
    // wrap of std cos
    static double sin(const Double& v);

    // wrap of std min
    static inline double min (const Double& v, const double p);
    // wrap of std min
    static inline double min (const double v, const Double& p);
    // wrap of std max
    static inline double max (const Double& v, const double p);
    // wrap of std max
    static inline double max (const double v, const Double& p);


    // numeric limits
    // --------------

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

      template< class Traits >
      friend OutStreamInterface< Traits > &
        operator<< ( OutStreamInterface< Traits > &, const Double );
      template< class Traits >
      friend InStreamInterface< Traits > &
        operator>> ( InStreamInterface< Traits > &, Double & );

      friend double pow (const Double& v, const double p);
      friend double pow (const Double& v, const Double& p);
      friend double log (const Double& v);
      friend double sqrt(const Double& v);
      friend double sin(const Double& v);
      friend double cos(const Double& v);

      friend Double abs ( const Double & );
      friend double min(const Double&, const double);
      friend double min(const double,  const Double&);
      friend double max(const Double&, const double);
      friend double max(const double,  const Double&);

      friend double real( const std::complex<Double>& );
      friend double real( const Double& );
      friend double imag( const std::complex<Double>& );
      friend double imag( const Double& );

#if DUNE_FEM_COMPATIBILITY
      friend struct XdrIO< Double >;
#endif

      friend void field_cast ( const Double &, double & );

    protected:
      typedef FlOpCounter< ThisType > FlOpCounterType;

    protected:
      double value_;

    public:
      operator double () const
      {
        return value_;
      }

      inline Double ()
        // : value_( 0 )
//#ifndef NDEBUG
//        : value_( std::numeric_limits< double >::signaling_NaN() )
//#endif
      {
      }

      inline Double ( const double value )
      : value_( value )
      {}

      inline Double ( const ThisType &other )
      : value_( other.value_ )
      {}

      inline ThisType &operator= ( const ThisType other )
      {
        // flOp();
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

      Double operator- () const
      {
        flOp();
        return Double( -value_ );
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

    // min/max
    // ---------

    // wrap of std min
    static inline double min (const Double& v, const double p)
    {
      return (v.value_ > p) ? p : v.value_;
    }

    // wrap of std min
    static inline double min (const double v, const Double& p)
    {
      return (v > p.value_) ? p.value_ : v;
    }

    // wrap of std max
    static inline double max (const Double& v, const double p)
    {
      return (v.value_ < p) ? p : v.value_;
    }

    // wrap of std max
    static inline double max (const double v, const Double& p)
    {
      return (v < p.value_) ? p.value_ : v;
    }

    // operator+
    // ---------

    inline Double operator+ ( const Double &a, const Double &b )
    {
      Double :: flOp();
      return Double( a.value_ + b.value_ );
    }

    inline Double operator+ ( const double a, const Double &b )
    {
      Double :: flOp();
      return Double( a + b.value_ );
    }

    inline Double operator+ ( const Double &a, const double b )
    {
      Double :: flOp();
      return Double( a.value_ + b );
    }

    inline Double operator+ ( const int a, const Double &b )
    {
      Double :: flOp();
      return Double( a + b.value_ );
    }

    inline Double operator+ ( const Double &a, const int b )
    {
      Double :: flOp();
      return Double( a.value_ + b );
    }

    inline Double operator+ ( const unsigned int a, const Double &b )
    {
      Double :: flOp();
      return Double( a + b.value_ );
    }

    inline Double operator+ ( const Double &a, const unsigned int b )
    {
      Double :: flOp();
      return Double( a.value_ + b );
    }



    // operator-
    // ---------

    inline Double operator- ( const Double &a, const Double &b )
    {
      Double :: flOp();
      return Double( a.value_ - b.value_ );
    }

    inline Double operator- ( const double a, const Double &b )
    {
      Double :: flOp();
      return Double( a - b.value_ );
    }

    inline Double operator- ( const Double &a, const double b )
    {
      Double :: flOp();
      return Double( a.value_ - b );
    }

    inline Double operator- ( const int a, const Double &b )
    {
      Double :: flOp();
      return Double( a - b.value_ );
    }

    inline Double operator- ( const Double &a, const int b )
    {
      Double :: flOp();
      return Double( a.value_ - b );
    }

    inline Double operator- ( const unsigned int a, const Double &b )
    {
      Double :: flOp();
      return Double( a - b.value_ );
    }

    inline Double operator- ( const Double &a, const unsigned int b )
    {
      Double :: flOp();
      return Double( a.value_ - b );
    }



    // operator*
    // ---------

    inline Double operator* ( const Double &a, const Double &b )
    {
      Double :: flOp();
      return Double( a.value_ * b.value_ );
    }

    inline Double operator* ( const double a, const Double &b )
    {
      Double :: flOp();
      return Double( a * b.value_ );
    }

    inline Double operator* ( const Double &a, const double b )
    {
      Double :: flOp();
      return Double( a.value_ * b );
    }

    inline Double operator* ( const int a, const Double &b )
    {
      Double :: flOp();
      return Double( a * b.value_ );
    }

    inline Double operator* ( const Double &a, const int b )
    {
      Double :: flOp();
      return Double( a.value_ * b );
    }

    inline Double operator* ( const unsigned int a, const Double &b )
    {
      Double :: flOp();
      return Double( a * b.value_ );
    }

    inline Double operator* ( const Double &a, const unsigned int b )
    {
      Double :: flOp();
      return Double( a.value_ * b );
    }



    // operator/
    // ---------

    inline Double operator/ ( const Double &a, const Double &b )
    {
      Double :: flOp();
      return Double( a.value_ / b.value_ );
    }

    inline Double operator/ ( const double a, const Double &b )
    {
      Double :: flOp();
      return Double( a / b.value_ );
    }

    inline Double operator/ ( const Double &a, const double b )
    {
      Double :: flOp();
      return Double( a.value_ / b );
    }

    inline Double operator/ ( const int a, const Double &b )
    {
      Double :: flOp();
      return Double( a / b.value_ );
    }

    inline Double operator/ ( const Double &a, const int b )
    {
      Double :: flOp();
      return Double( a.value_ / b );
    }

    inline Double operator/ ( const unsigned int a, const Double &b )
    {
      Double :: flOp();
      return Double( a / b.value_ );
    }

    inline Double operator/ ( const Double &a, const unsigned int b )
    {
      Double :: flOp();
      return Double( a.value_ / b );
    }



    // operator==
    // ----------

    inline bool operator== ( const Double &a, const Double &b )
    {
      return (a.value_ == b.value_);
    }

    inline bool operator== ( const double a, const Double &b )
    {
      return (a == b.value_);
    }

    inline bool operator== ( const Double &a, const double b )
    {
      return (a.value_ == b);
    }

    inline bool operator== ( const int a, const Double &b )
    {
      return (a == b.value_);
    }

    inline bool operator== ( const Double &a, const int b )
    {
      return (a.value_ == b);
    }

    inline bool operator== ( const unsigned int a, const Double &b )
    {
      return (a == b.value_);
    }

    inline bool operator== ( const Double &a, const unsigned int b )
    {
      return (a.value_ == b);
    }



    // operator!=
    // ----------

    inline bool operator!= ( const Double &a, const Double &b )
    {
      return (a.value_ != b.value_);
    }

    inline bool operator!= ( const double a, const Double &b )
    {
      return (a != b.value_);
    }

    inline bool operator!= ( const Double &a, const double b )
    {
      return (a.value_ != b);
    }

    inline bool operator!= ( const int a, const Double &b )
    {
      return (a != b.value_);
    }

    inline bool operator!= ( const Double &a, const int b )
    {
      return (a.value_ != b);
    }

    inline bool operator!= ( const unsigned int a, const Double &b )
    {
      return (a != b.value_);
    }

    inline bool operator!= ( const Double &a, const unsigned int b )
    {
      return (a.value_ != b);
    }



    // operator<
    // ---------

    inline bool operator< ( const Double &a, const Double &b )
    {
      return (a.value_ < b.value_);
    }

    inline bool operator< ( const double a, const Double &b )
    {
      return (a < b.value_);
    }

    inline bool operator< ( const Double &a, const double b )
    {
      return (a.value_ < b);
    }

    inline bool operator< ( const int a, const Double &b )
    {
      return (a < b.value_);
    }

    inline bool operator< ( const Double &a, const int b )
    {
      return (a.value_ < b);
    }

    inline bool operator< ( const unsigned int a, const Double &b )
    {
      return (a < b.value_);
    }

    inline bool operator< ( const Double &a, const unsigned int b )
    {
      return (a.value_ < b);
    }



    // operator<=
    // ----------

    inline bool operator<= ( const Double &a, const Double &b )
    {
      return (a.value_ <= b.value_);
    }

    inline bool operator<= ( const double a, const Double &b )
    {
      return (a <= b.value_);
    }

    inline bool operator<= ( const Double &a, const double b )
    {
      return (a.value_ <= b);
    }

    inline bool operator<= ( const int a, const Double &b )
    {
      return (a <= b.value_);
    }

    inline bool operator<= ( const Double &a, const int b )
    {
      return (a.value_ <= b);
    }

    inline bool operator<= ( const unsigned int a, const Double &b )
    {
      return (a <= b.value_);
    }

    inline bool operator<= ( const Double &a, const unsigned int b )
    {
      return (a.value_ <= b);
    }



    // operator>
    // ---------

    inline bool operator> ( const Double &a, const Double &b )
    {
      return (a.value_ > b.value_);
    }

    inline bool operator> ( const double a, const Double &b )
    {
      return (a > b.value_);
    }

    inline bool operator> ( const Double &a, const double b )
    {
      return (a.value_ > b);
    }

    inline bool operator> ( const int a, const Double &b )
    {
      return (a > b.value_);
    }

    inline bool operator> ( const Double &a, const int b )
    {
      return (a.value_ > b);
    }

    inline bool operator> ( const unsigned int a, const Double &b )
    {
      return (a > b.value_);
    }

    inline bool operator> ( const Double &a, const unsigned int b )
    {
      return (a.value_ > b);
    }



    // operator>=
    // ----------

    inline bool operator>= ( const Double &a, const Double &b )
    {
      return (a.value_ >= b.value_);
    }

    inline bool operator>= ( const double a, const Double &b )
    {
      return (a >= b.value_);
    }

    inline bool operator>= ( const Double &a, const double b )
    {
      return (a.value_ >= b);
    }

    inline bool operator>= ( const int a, const Double &b )
    {
      return (a >= b.value_);
    }

    inline bool operator>= ( const Double &a, const int b )
    {
      return (a.value_ >= b);
    }

    inline bool operator>= ( const unsigned int a, const Double &b )
    {
      return (a >= b.value_);
    }

    inline bool operator>= ( const Double &a, const unsigned int b )
    {
      return (a.value_ >= b);
    }



    // stream operators
    // ----------------

    inline std :: ostream &operator<< ( std :: ostream &out, const Double &a )
    {
      return out << a.value_;
    }

    inline std :: istream &operator>> ( std :: istream &in, Double &a )
    {
      return in >> a.value_;
    }

    template< class Traits >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const Double a )
    {
      return out << a.value_;
    }

    template< class Traits >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   Double &a )
    {
      return in >> a.value_;
    }



    // standard functions
    // ------------------

    inline Double abs ( const Double &a )
    {
      return Double( std::abs( a.value_ ) );
    }

    static inline double log (const Double& v)
    {
      return std::log(v.value_);
    }

    inline double pow (const Double& a, const double b )
    {
      return std::pow(a.value_, b);
    }

    static inline double sqrt(const Double& v)
    {
      return std::sqrt(v.value_);
    }

    static inline double sin (const Double& v)
    {
      return std::sin(v.value_);
    }

    static inline double cos(const Double& v)
    {
      return std::cos(v.value_);
    }

    inline void field_cast ( const Double &f1, double &f2 )
    {
      f2 = f1.value_;
    }

    inline double real (const std::complex<Double>& x)
    {
      return x.real().value_;
    }

    inline double real (const Double& x)
    {
      return x.value_;
    }

    inline double imag (const std::complex<Double>& x)
    {
      return x.imag().value_;
    }

    inline double imag (const Double& x)
    {
      return x.value_;
    }

  } // namespace Fem

  using Fem :: Double ;

  template <>
  struct IsNumber< Double > : public IsNumber< double > {};
#if DUNE_VERSION_GT(DUNE_COMMON, 2, 6)
  template <>
  struct HasNaN < Double > : public HasNaN < double > {};
#endif

} // namespace Dune


namespace std
{
  inline Dune::Fem::Double abs ( const Dune::Fem::Double &a )
  {
    return Dune::Fem::abs( a );
  }

  inline Dune::Fem::Double pow( const Dune::Fem::Double &a, const Dune::Fem::Double &b )
  {
    return Dune::Fem::pow( a, b );
  }

  // wrap of std min
  inline double min (const Dune::Fem::Double& v, const double p)
  {
    return Dune::Fem::min(v,p);
  }

  // wrap of std min
  inline double min (const double v, const Dune::Fem::Double& p)
  {
    return Dune::Fem::min(v,p);
  }

  // wrap of std max
  inline double max (const Dune::Fem::Double& v, const double p)
  {
    return Dune::Fem::max(v,p);
  }

  // wrap of std max
  inline double max (const double v, const Dune::Fem::Double& p)
  {
    return Dune::Fem::max(v,p);
  }

  // wrap of std sqrt
  inline double sqrt( const Dune::Fem::Double& v )
  {
    return Dune::Fem::sqrt( v );
  }

  // wrap of std real
  inline double real (const complex<Dune::Fem::Double>& x)
  {
    return Dune::Fem::real( x );
  }

  // wrap of std real
  inline double real (const Dune::Fem::Double& x)
  {
    return Dune::Fem::real( x );
  }

  // wrap of std imag
  inline double imag (const complex<Dune::Fem::Double>& x)
  {
    return Dune::Fem::imag( x );
  }

  // wrap of std imag
  inline double imag (const Dune::Fem::Double& x)
  {
    return Dune::Fem::imag( x );
  }



  // numeric limits
  // --------------

  template<>
  struct numeric_limits< Dune::Fem::Double >
  {
    static const bool is_specialized = true;

    static const int radix = numeric_limits< double > :: radix;
    static const int digits = numeric_limits< double > :: digits;
    static const int digits10 = numeric_limits< double > :: digits10;

    static const bool is_signed = numeric_limits< double > :: is_signed;
    static const bool is_integer = numeric_limits< double > :: is_integer;
    static const bool is_exact = numeric_limits< double > :: is_exact;

    inline static Dune::Fem::Double min () throw ()
    {
      return Dune::Fem::Double( numeric_limits< double > :: min() );
    }

    inline static Dune::Fem::Double max () throw ()
    {
      return Dune::Fem::Double( numeric_limits< double > :: max() );
    }

    inline static Dune::Fem::Double epsilon () throw ()
    {
      return Dune::Fem::Double( numeric_limits< double > :: epsilon() );
    }

    inline static Dune::Fem::Double round_error () throw ()
    {
      return Dune::Fem::Double( numeric_limits< double > :: round_error() );
    }

    inline static Dune::Fem::Double infinity () throw ()
    {
      return Dune::Fem::Double( numeric_limits< double > :: infinity() );
    }

    inline static Dune::Fem::Double quiet_NaN () throw ()
    {
      return Dune::Fem::Double( numeric_limits< double > :: quiet_NaN() );
    }

    inline static Dune::Fem::Double signaling_NaN () throw ()
    {
      return Dune::Fem::Double( numeric_limits< double > :: signaling_NaN() );
    }

    inline static Dune::Fem::Double denorm_min () throw ()
    {
      return Dune::Fem::Double( numeric_limits< double > :: denorm_min() );
    }

    static const int min_exponent = numeric_limits< double > :: min_exponent;
    static const int max_exponent = numeric_limits< double > :: max_exponent;
    static const int min_exponent10 = numeric_limits< double > :: min_exponent10;
    static const int max_exponent10 = numeric_limits< double > :: max_exponent10;

    static const bool has_infinity = numeric_limits< double > :: has_infinity;
    static const bool has_quiet_NaN = numeric_limits< double > :: has_quiet_NaN;
    static const bool has_signaling_NaN = numeric_limits< double > :: has_signaling_NaN;
    static const float_denorm_style has_denorm = numeric_limits< double > :: has_denorm;
    static const bool has_denorm_loss = numeric_limits< double > :: has_denorm_loss;

    static const bool is_iec559 = numeric_limits< double > :: is_iec559;
    static const bool is_bounded = numeric_limits< double > :: is_bounded;
    static const bool is_modulo = numeric_limits< double > :: is_modulo;

    static const bool traps = numeric_limits< double > :: traps;
    static const bool tinyness_before = numeric_limits< double > :: tinyness_before;
    static const float_round_style round_style
      = numeric_limits< double > :: round_style;

  };

  template <>
  struct is_floating_point< Dune::Fem::Double > : public is_floating_point< double > {};

} // namespace std

#endif // #ifndef DUNE_FEM_DOUBLE_HH

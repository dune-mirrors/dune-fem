#ifndef DUNE_FEM_VALIDATOR_HH
#define DUNE_FEM_VALIDATOR_HH

namespace Dune
{
  
  template< class T, class Impl >
  class ValidatorDefault;



  template< class T, class Impl >
  class ValidatorInterface
  {
    typedef ValidatorInterface< T, Impl > ThisType;

    friend class ValidatorDefault< T, Impl >;

  private:
    inline ValidatorInterface () {}

    ValidatorInterface ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    inline bool operator () ( const T &value ) const
    {
      return asImp()( value );
    }
    
  protected:
    inline const Impl &asImp () const
    {
      return static_cast< const Impl & >( *this );
    }

    inline Impl &asImp ()
    {
      return static_cast< Impl & >( *this );
    }
  };



  template< class T, class Impl >
  class ValidatorDefault
  {
    typedef ValidatorDefault< T, Impl > ThisType;
    typedef ValidatorInterface< T, Impl > BaseType;

  protected:
    inline ValidatorDefault () {}

  private:
    ValidatorDefault ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  private:
    bool operator () ( const T &value ) const;
  };


  
  template< class T >
  class ValidateGreater
  : public ValidatorDefault< T, ValidateGreater< T > >
  {
    typedef ValidateGreater< T > ThisType;
    typedef ValidatorDefault< T, ThisType > BaseType;

  protected:
    const T threshold_;

  public:
    inline ValidateGreater ( const T &threshold )
    : threshold_( threshold )
    {}

    inline ValidateGreater ( const ThisType &other )
    : threshold_( other.threshold_ )
    {}

  public:
    inline bool operator() ( const T &value ) const
    {
      return value > threshold_;
    }
  };



  template< class T >
  class ValidateLess
  : public ValidatorDefault< T, ValidateLess< T > >
  {
    typedef ValidateLess< T > ThisType;
    typedef ValidatorDefault< T, ThisType > BaseType;

  protected:
    const T threshold_;

  public:
    inline ValidateLess ( const T &threshold )
    : threshold_( threshold )
    {}

    inline ValidateLess ( const ThisType &other )
    : threshold_( other.threshold_ )
    {}

    inline bool operator() ( const T &value ) const
    {
      return value < threshold_;
    }
  };



  template< class T >
  class ValidateNotGreater
  : public ValidatorDefault< T, ValidateNotGreater< T > >
  {
    typedef ValidateNotGreater< T > ThisType;
    typedef ValidatorDefault< T, ThisType > BaseType;

  protected:
    const T threshold_;

  public:
    inline ValidateNotGreater ( const T &threshold )
    : threshold_( threshold )
    {}

    inline ValidateNotGreater ( const ThisType &other )
    : threshold_( other.threshold_ )
    {}

    inline bool operator() ( const T &value ) const
    {
      return value <= threshold_;
    }
  };



  template< class T >
  class ValidateNotLess
  : public ValidatorDefault< T, ValidateNotLess< T > >
  {
    typedef ValidateNotLess< T > ThisType;
    typedef ValidatorDefault< T, ThisType > BaseType;

  protected:
    const T threshold_;

  public:
    inline ValidateNotLess ( const T &threshold )
    : threshold_( threshold )
    {}

    inline ValidateNotLess ( const ThisType &other )
    : threshold_( other.threshold_ )
    {}

    inline bool operator() ( const T &value ) const
    {
      return value >= threshold_;
    }
  };



  template< class T, bool leftClosed, bool rightClosed >
  class ValidateInterval
  : public ValidatorDefault< T, ValidateInterval< T, leftClosed, rightClosed > >
  {
    typedef ValidateInterval< T, leftClosed, rightClosed > ThisType;
    typedef ValidatorDefault< T, ThisType > BaseType;

  protected:
    const T lThreshold_,rThreshold_;

  public:
    inline ValidateInterval ( const T &lThreshold, const T &rThreshold )
    : lThreshold_( lThreshold ),
      rThreshold_( rThreshold )
    {}

    inline ValidateInterval ( const ThisType &other )
    : lThreshold_( other.lThreshold_ ),
      rThreshold_( other.rThreshold_ )
    {}

    inline bool operator() ( const T &value ) const
    {
      bool ret = true;
      ret &= (leftClosed  ? value >= lThreshold_ : value > lThreshold_);
      ret &= (rightClosed ? value <= rThreshold_ : value < rThreshold_);
      return ret;
    }
  };

}

#endif

#ifndef DUNE_FEM_VALIDATOR_HH
#define DUNE_FEM_VALIDATOR_HH

#include <iostream>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class T, class Impl >
    class ValidatorDefault;



    // ValidatorInterface
    // ------------------

    template< class T, class Impl >
    class ValidatorInterface
    {
      typedef ValidatorInterface< T, Impl > ThisType;

      friend class ValidatorDefault< T, Impl >;

    private:
      ValidatorInterface () {}

      ValidatorInterface ( const ThisType & );
      ThisType &operator= ( const ThisType & );

    public:
      bool operator () ( const T &value ) const
      {
        return asImp()( value );
      }

      void print(std::ostream& s) const
      {
        asImp().print( s );
      }
      
    protected:
      const Impl &asImp () const
      {
        return static_cast< const Impl & >( *this );
      }

      Impl &asImp ()
      {
        return static_cast< Impl & >( *this );
      }
    };



    // ValidatorDefault
    // ----------------

    template< class T, class Impl >
    class ValidatorDefault
    {
      typedef ValidatorDefault< T, Impl > ThisType;
      typedef ValidatorInterface< T, Impl > BaseType;

    protected:
      ValidatorDefault () {}

    private:
      ValidatorDefault ( const ThisType & );
      ThisType &operator= ( const ThisType & );

    private:
      bool operator () ( const T &value ) const;
      void print ( std::ostream &s ) const;
    };


    
    template< class T >
    class ValidateGreater
    : public ValidatorDefault< T, ValidateGreater< T > >
    {
      typedef ValidateGreater< T > ThisType;
      typedef ValidatorDefault< T, ThisType > BaseType;

    public:
      ValidateGreater ( const T &threshold )
      : threshold_( threshold )
      {}

      ValidateGreater ( const ThisType &other )
      : threshold_( other.threshold_ )
      {}

    public:
      bool operator() ( const T &value ) const
      {
        return value > threshold_;
      }

      void print(std::ostream& s) const
      {
        s << "ValidateLess: valid values are: > " << threshold_ << std::endl << std::endl;
      }

    protected:
      const T threshold_;
    };



    template< class T >
    class ValidateLess
    : public ValidatorDefault< T, ValidateLess< T > >
    {
      typedef ValidateLess< T > ThisType;
      typedef ValidatorDefault< T, ThisType > BaseType;

    public:
      ValidateLess ( const T &threshold )
      : threshold_( threshold )
      {}

      ValidateLess ( const ThisType &other )
      : threshold_( other.threshold_ )
      {}

      bool operator() ( const T &value ) const
      {
        return value < threshold_;
      }

      void print(std::ostream& s) const
      {
        s << "ValidateLess: valid values are: < " << threshold_ << std::endl << std::endl;
      }

    protected:
      const T threshold_;
    };



    template< class T >
    class ValidateNotGreater
    : public ValidatorDefault< T, ValidateNotGreater< T > >
    {
      typedef ValidateNotGreater< T > ThisType;
      typedef ValidatorDefault< T, ThisType > BaseType;

    public:
      ValidateNotGreater ( const T &threshold )
      : threshold_( threshold )
      {}

      ValidateNotGreater ( const ThisType &other )
      : threshold_( other.threshold_ )
      {}

      bool operator() ( const T &value ) const
      {
        return value <= threshold_;
      }

      void print(std::ostream& s) const
      {
        s << "ValidateNotGreater: valid values are: <= " << threshold_ << std::endl << std::endl;
      }

    protected:
      const T threshold_;
    };



    template< class T >
    class ValidateNotLess
    : public ValidatorDefault< T, ValidateNotLess< T > >
    {
      typedef ValidateNotLess< T > ThisType;
      typedef ValidatorDefault< T, ThisType > BaseType;

    public:
      ValidateNotLess ( const T &threshold )
      : threshold_( threshold )
      {}

      ValidateNotLess ( const ThisType &other )
      : threshold_( other.threshold_ )
      {}

      bool operator() ( const T &value ) const
      {
        return value >= threshold_;
      }

      void print(std::ostream& s) const
      {
        s << "ValidateNotLess: valid values are: >= " << threshold_ << std::endl << std::endl;
      }

    protected:
      const T threshold_;
    };



    template< class T, bool leftClosed, bool rightClosed >
    class ValidateInterval
    : public ValidatorDefault< T, ValidateInterval< T, leftClosed, rightClosed > >
    {
      typedef ValidateInterval< T, leftClosed, rightClosed > ThisType;
      typedef ValidatorDefault< T, ThisType > BaseType;

    public:
      ValidateInterval ( const T &lThreshold, const T &rThreshold )
      : lThreshold_( lThreshold ),
        rThreshold_( rThreshold )
      {}

      ValidateInterval ( const ThisType &other )
      : lThreshold_( other.lThreshold_ ),
        rThreshold_( other.rThreshold_ )
      {}

      bool operator() ( const T &value ) const
      {
        bool ret = true;
        ret &= (leftClosed  ? value >= lThreshold_ : value > lThreshold_);
        ret &= (rightClosed ? value <= rThreshold_ : value < rThreshold_);
        return ret;
      }

      void print(std::ostream& s) const
      {
        const char* left  = (leftClosed)  ? "[" : "(";
        const char* right = (rightClosed) ? "]" : ")";
        s << "ValidateInterval: valid values are " << left << lThreshold_ << "," <<
                rThreshold_ << right << std::endl << std::endl;
      }

    protected:
      const T lThreshold_, rThreshold_;
    };



    // NoWhiteSpaceValidator
    // ---------------------

    class NoWhiteSpaceValidator
    : public ValidatorDefault< std::string, NoWhiteSpaceValidator >
    {
      typedef NoWhiteSpaceValidator ThisType;
      typedef ValidatorDefault< std::string, ThisType > BaseType;

    public:
      NoWhiteSpaceValidator ()
      {}

      NoWhiteSpaceValidator ( const ThisType &other )
      {}

      bool operator() ( const std::string &value ) const
      {
        return (value.find_first_of( " \t" ) == std::string::npos);
      }

      void print ( std::ostream &s ) const
      {
        s << "NoWhiteSpaceValidator" << std::endl;
      }
    };

  } // namespace Fem

} // namespace Dune 

#endif // #ifndef DUNE_FEM_VALIDATOR_HH

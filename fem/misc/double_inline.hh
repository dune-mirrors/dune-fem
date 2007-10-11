namespace Dune
{

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
 
}

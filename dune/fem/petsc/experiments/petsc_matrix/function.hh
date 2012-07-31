#ifndef FUNCTION_HH
#define FUNCTION_HH

#include <cmath>
#include <assert.h>

/*
 * Definition of special functions
 */

// Abstract base class for functions on [a,b]. S is the scalar type
template< typename S >
class Function {
public:
    typedef S ScalarType;
    
    // ctor/dtor
    Function ( const ScalarType &a, const ScalarType &b )
    : a_( a ),
      b_( b )
    {
        assert( a_ < b_ );
    }

    virtual ~Function () {}

    // eval the function at x
    virtual ScalarType evaluate ( const ScalarType &x ) const = 0;

    // get the lower bound a and the upper bound b of [a,b]
    ScalarType lowerBound () const { return a_; }
    ScalarType upperBound () const { return b_; }

protected:
    // data fields
    ScalarType a_;
    ScalarType b_;

private:
    // prohibited methods
    Function ();
};


// A sine which is zero on the boundary of [a, b]
template< typename S > 
class FunctionSine : public Function< S > {
    typedef Function< S > BaseType;

public:
    typedef typename BaseType::ScalarType ScalarType;

    FunctionSine ( const ScalarType &a, const ScalarType &b )
    : BaseType( a, b ),
      factor_( 2*M_PI/( b - a ) ),
      offset_( a*2*M_PI/( b - a ) )
    {}
    
    ScalarType evaluate ( const ScalarType &x ) const {
        assert( this->a_ < x + 1.e-12 );
        assert( this->b_ > x - 1.e-12 );

        return std::sin( x*factor_ - offset_ );
    }

    ScalarType getFactor () const { return factor_; }
    ScalarType getOffset () const { return offset_; }

private:
    // prohibit standard methods
    FunctionSine ();


    ScalarType factor_;
    ScalarType offset_;
};

template< typename S >
class SolutionSine : public Function< S > {
    typedef Function< S > BaseType;

public:
    typedef S ScalarType;


    SolutionSine ( const ScalarType &a, const ScalarType &b )
    : BaseType( a, b ),
      sine_( a, b ),
      scalingFactor_( 1./ (sine_.getFactor() * sine_.getFactor() ) )
    {
        assert( a < b );
    }

    ScalarType evaluate ( const ScalarType &x ) const { return scalingFactor_ * sine_.evaluate( x ); }

private:
    // prohibited methods
    SolutionSine ();


    FunctionSine< S > sine_;
    ScalarType scalingFactor_;
};


#endif // FUNCTION_HH

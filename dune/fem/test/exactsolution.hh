#ifndef DUNE_FEM_TEST_EXACTSOLUTION_HH
#define DUNE_FEM_TEST_EXACTSOLUTION_HH

#include <cmath>

#include <complex>

#include <dune/fem/function/common/function.hh>

namespace Dune
{

  namespace Fem
  {

    // ExactSolution
    // -------------

    template< class FunctionSpaceImp >
    class ExactSolution
      : public Fem::Function< FunctionSpaceImp, ExactSolution< FunctionSpaceImp > >
    {
      typedef ExactSolution< FunctionSpaceImp > ThisType;
      typedef Fem::Function< FunctionSpaceImp, ThisType > BaseType;

    public:
      typedef FunctionSpaceImp FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

    public:
      void evaluate( const DomainType &x, RangeType &phi ) const
      {
        using std::sin; using std::pow;
        phi = 1;
        for( int r = 0; r < RangeType::dimension; ++r )
        {
          for( int i = 0; i < DomainType::dimension; ++i )
          {
            const double s = sin( M_PI * x[ i ] );
            phi[ r ] += pow( s, double( r+1 ) );
          }
        }
        phi *= factor( RangeFieldType( 1 ) );
      }

      template< class T >
      static std::complex< T > factor ( std::complex< T > alpha ) noexcept
      {
        return alpha * std::complex< T >( 1, -2 );
      }

      template< class T >
      static T factor ( T alpha ) noexcept
      {
        return alpha;
      }

      void evaluate( const DomainType &x, RangeFieldType t, RangeType &phi ) const
      {
        evaluate( x, phi );
      }

      void jacobian( const DomainType &x, JacobianRangeType &Dphi ) const
      {
        using std::cos; using std::sin; using std::pow;
        Dphi = 0;
        for( int r = 0; r < RangeType::dimension; ++r )
        {
          for( int i = 0; i < DomainType::dimension; ++i )
          {
            const double s = sin( M_PI * x[ i ] );
            const double c = cos( M_PI * x[ i ] );
            Dphi[ r ][ i ] += M_PI * double( r+1 ) * pow( s, double( r ) ) * c;
          }
        }
        Dphi *= factor( RangeFieldType( 1 ) );
      }

      void jacobian( const DomainType &x, RangeFieldType t, JacobianRangeType &Dphi ) const
      {
        jacobian( x, Dphi );
      }

      void hessian( const DomainType &x, HessianRangeType &H ) const
      {
        using std::cos; using std::sin; using std::pow;

        for( int r = 0; r < RangeType::dimension; ++r )
        {
          H[ r ] = RangeFieldType( 0 );
          for( int i = 0; i < DomainType::dimension; ++i )
          {
            const double s = sin( M_PI * x[ i ] );
            const double c = cos( M_PI * x[ i ] );
            H[ r ][ i ][ i ] = M_PI * M_PI * double( r+1 ) * (r * std::pow( s, double( r-1 ) ) * c*c - std::pow( s, double( r+1 ) ));
          }
        }
        H *= factor( RangeFieldType( 1 ) );
      }

      void hessian( const DomainType &x, RangeFieldType t, HessianRangeType &H ) const
      {
        hessian( x, H );
      }
    };

    // ExactSolution
    // -------------

    template< class FunctionSpaceImp>
    class Polynomial
      : public Fem::Function< FunctionSpaceImp, ExactSolution< FunctionSpaceImp > >
    {
      typedef Polynomial< FunctionSpaceImp > ThisType;
      typedef Fem::Function< FunctionSpaceImp, ThisType > BaseType;

    public:
      typedef FunctionSpaceImp FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;
    protected:
      int order_;

    public:
      Polynomial( const int order ) : order_( order ) {}

      void evaluate( const DomainType &x, RangeType &phi ) const
      {
        phi = 1;
        for( int i=0; i<DomainType::dimension; ++i )
        {
          RangeType val( 1. );
          for( int p=0; p<order_; ++p )
            val *= x[ i ];

          phi += val;
        }
      }

      void evaluate( const DomainType &x, RangeFieldType t, RangeType &phi ) const
      {
        evaluate( x, phi );
      }

      void jacobian( const DomainType &x, JacobianRangeType &Dphi ) const
      {
        Dphi = 0;
        for( int i=0; i<DomainType::dimension; ++i )
        {
          JacobianRangeType val = order_;
          for( int p=0; p<order_-1; ++p )
            val *= x[ i ];

          Dphi += val;
        }
      }

      void jacobian( const DomainType &x, RangeFieldType t, JacobianRangeType &Dphi ) const
      {
        jacobian( x, Dphi );
      }

      void hessian( const DomainType &x, HessianRangeType &H ) const
      {
        H = 0 ;
      }

      void hessian( const DomainType &x, RangeFieldType t, HessianRangeType &H ) const
      {
        hessian( x, H );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_TEST_EXACTSOLUTION_HH

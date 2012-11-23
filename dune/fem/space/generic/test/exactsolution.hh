#ifndef DUNE_FEM_SPACE_GENERICDISCRETE_TEST_EXACTSOLUTION_HH
#define DUNE_FEM_SPACE_GENERICDISCRETE_TEST_EXACTSOLUTION_HH


namespace Dune
{

  namespace Fem
  {

    // ExactSolution
    // -------------

    template< class FunctionSpace >
    struct ExactSolution
    {

      typedef FunctionSpace FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

      static void evaluate ( const DomainType &x, RangeType &val )
      {
        val = RangeFieldType( 1 );
        for( int r = 0; r < RangeType::dimension; ++r )
          for( int i = 0; i < DomainType::dimension; ++i )
            val[ r ] += pow( cos( M_PI * x[ i ] ), double(r+1) );
      }

      static void evaluate ( const DomainType &x, double t, RangeType &val )
      {
        return evaluate( x, val );
      }

      static void jacobian ( const DomainType &x, JacobianRangeType &jac )
      {
        jac = 1;
        for (int r = 0; r < RangeType::dimension; ++r)
        {
          for( int i = 0; i < DomainType::dimension; ++i )
          {
            for( int j = 0; j < DomainType::dimension; ++j )
            {
              jac[ r ][ j ] -= double(r+1)*pow(cos( M_PI * x[ i ] ),double(r))*
                ((i != j) ? 0 : M_PI * sin( M_PI * x[ i ] ));
            }
          }
        }
      }

      void jacobian ( const DomainType &x, double t, JacobianRangeType &jac ) const
      {
        jacobian( x, jac );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_GENERICDISCRETE_TEST_EXACTSOLUTION_HH

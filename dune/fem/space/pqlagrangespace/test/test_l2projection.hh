#ifndef  DUNE_FEM_P12DSPACE_TEST_L2PROJECTION_HH__
#define  DUNE_FEM_P12DSPACE_TEST_L2PROJECTION_HH__


// Includes from DUNE-FEM
// ----------------------

// include support for versioning
#include <dune/fem/version.hh>
// include parameter handling
#include <dune/fem/io/parameter.hh>
// include basic grid parts
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
// include Lagrange discrete function space
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/pqlagrangespace/pqlagrangespace.hh>
#include <dune/fem/space/pqlagrangespace/genericlagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/misc/l2error.hh>
#include "massoperator.hh"

namespace Dune
{

template<class GridType, int polOrder, bool isCube=false>
class L2Projection_Test
  : public Test
{
public:
  static const unsigned int dimension = GridType :: dimension;
  // select grid part to use
  //typedef Dune::LeafGridPart< GridType > GridPartType;
  typedef Dune::AdaptiveLeafGridPart< GridType > GridPartType;

  typedef FunctionSpace< double, double, dimension, 1 >                FunctionSpaceType;

/*  typedef Dune :: P12DSpace< FunctionSpaceType, GridPartType >       DiscreteFunctionSpaceType;*/

#if 0 // USE_PQSPACE
  /*
  typedef typename SelectType< isCube,
                               QLagrangeSpace< FunctionSpaceType,
                                               GridPartType,
                                               polOrder >,
                               PLagrangeSpace< FunctionSpaceType,
                                               GridPartType,
                                               polOrder > > :: Type  DiscreteFunctionSpaceType; 
  */
  typedef class PQ22DLagrangeSpace< FunctionSpaceType, GridPartType > DiscreteFunctionSpaceType;
#else
  typedef LagrangeSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;
#endif
  typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >      DiscreteFunctionType;
private:
  template< class FunctionSpace >
  struct MyFunction
  : public Fem::Function< FunctionSpace, MyFunction< FunctionSpace > >
  {
    typedef typename FunctionSpace :: DomainType                     DomainType;
    typedef typename FunctionSpace :: RangeType                      RangeType;
  
    void evaluate ( const DomainType &x, RangeType &value ) const
    {
      value = 1.;
      for( int i = 0; i < FunctionSpace::dimDomain; ++i )
        value *= sin( M_PI * x[ i ] );
    }
  
    void evaluate ( const DomainType &x, const double &t, RangeType &value ) const
    {
      evaluate( x, value );
    }
  };

  template< class Function >
  double algorithm ( const Function &function, GridType &grid, int repeat )
  {
    std::cout << "starting algorithm run no. " << repeat << std::endl;
    typedef Dune :: MassOperator< DiscreteFunctionType >             MassOperatorType;
    typedef Dune :: CGInverseOp< DiscreteFunctionType,
                                 MassOperatorType >                  InverseOperator;

    GridPartType gridPart( grid );
    DiscreteFunctionSpaceType dfSpace( gridPart );

    MassOperatorType massOperator( dfSpace );
    DiscreteFunctionType functional( "functional", dfSpace );
    massOperator( function, functional );

    DiscreteFunctionType solution( "solution", dfSpace );
    solution.clear();
    InverseOperator inverseOperator( massOperator, 1e-10, 1e-10, 50000, 0 );
    inverseOperator( functional, solution );

    Dune::L2Error< DiscreteFunctionType > l2error;
    return l2error.norm( function, solution );
  }

public:
  L2Projection_Test(std::string gridFile)
    : gridFile_(gridFile), level_(0), repeats_(3)
  { }

  void run() 
  {
    try
    {
      GridPtr< GridType > gridPtr( gridFile_ );
      GridType& grid = *gridPtr;
      grid.globalRefine( level_ );

      std::cout << std::endl << std::endl;
      std::cout << "=============================================================" << std::endl;
      std::cout << "Starting L2 projection test with EOC measurements for" << std::endl;
      std::cout << "base functions of polynomial order " << polOrder << std::endl;
      std::cout << "on a grid of dimension             " << GridType::dimension << std::endl;
      std::cout << "=============================================================" << std::endl;

      MyFunction< FunctionSpaceType > function;
      GridPartType gridPart( grid );

      double error = algorithm( function, grid, 0 );
      std::cout << error << std::endl;
      double eoc = 0;
      for( int i = 0; i < repeats_; ++i )
      {
        const double prevError = error;
        Dune::GlobalRefine::apply( grid, Dune::DGFGridInfo< GridType >::refineStepsForHalf() );
        error = algorithm( function, grid, i+1 );
        eoc   = log( prevError / error ) / M_LN2;

        std::cout << "error: " << error << " EOC: " << eoc << std::endl;
      }
      _test(eoc > 0.5 + polOrder);
    }
    catch( const Dune::Exception &e )
    {
      std::cerr << e << std::endl;
    }
    catch(...)
    {
      std::cerr << "Generic exception!" << std::endl;
    }       
  }

private:
  std::string gridFile_;
  const int   level_;
  const int   repeats_;
};

};

#endif  /*__TEST_L2PROJECTION_HH__*/

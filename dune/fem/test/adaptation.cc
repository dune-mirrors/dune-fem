#include <config.h>

// iostream includes
#include <iostream>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/hierarchicgridpart.hh>

// include output
#include <dune/fem/io/file/dataoutput.hh>

// include discrete function space
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/space/fvspace.hh>

// adaptation ...
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/space/common/adaptmanager.hh>

// Scheme
// ------

template < class GridPart, class FunctionSpace > 
struct Scheme
{
  typedef GridPart GridPartType;
  typedef typename GridPartType::GridType GridType;
  
  typedef FunctionSpace FunctionSpaceType;
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 > DiscreteFunctionSpaceType;
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

  typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >  RestrictionProlongationType;

  typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType > AdaptationManagerType;

  typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
  typedef typename GridType :: template Codim< 0 > :: Entity ElementType;
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType; 

  Scheme( GridPartType &gridPart )
    : gridPart_( gridPart ),
      grid_( gridPart_.grid() ),
      discreteSpace_( gridPart_ ),
      solution_( "solution", discreteSpace_ ),
      restrictProlong_( solution_ ),
      adaptationManager_( gridPart_.grid(), restrictProlong_ )
  {
    solution_.clear();
  }

  const DiscreteFunctionType &solution() const
  {
    return solution_;
  }

  //! mark elements for adaptation 
  bool mark ( const double time ) const
  {
    int marked = 0;
    int count = 0;
    int total = 0;

	  // grid_.mark( 1, *(discreteSpace_.begin()) ); ++marked;
    // loop over all elements 
    const IteratorType end = discreteSpace_.end();
    for( IteratorType it = discreteSpace_.begin(); it != end; ++it )
    {
	    const ElementType &entity = *it;

      //    find center
      // DomainType center = entity.geometry().center();
      //    x = center - ( t, t, t )
      // DomainType x(time);
      // x -= center;

      DomainType center(0);
      for( int i = 0; i < entity.geometry().corners(); ++i )
      {
        center += entity.geometry().corner( i );
      }
      DomainType x;
      for( int i = 0; i < 3; ++i )
        x[ i ] = center[ i ] - time;

      // refine if 0.3 < |x| < 1.0, otherwise (possibly) coarsen
      if( x.two_norm() > 0.3 && x.two_norm() < 1.0 && entity.level() <= 9 )
      {
        grid_.mark( 1, entity );
        marked = 1;
        count ++;
      }
      else
      {
        // grid_.mark( -1, entity );
      }

      total++;
    }

    // get global max 
    marked = grid_.comm().max( marked );

    // print info
    if( bool( marked ) )
      std::cout << "P" << Dune::Fem::MPIManager::rank() << ": " 
                << "marked (" << count << " of " << total << ")" << std::endl;
    return (marked >= 0);
  }

  //! do the adaptation for a given marking 
  void adapt() 
  {
    // apply adaptation and load balancing 
    std::cout << "P" << Dune::Fem::MPIManager::rank() << ": " 
              << "Calling adapt:" << std::endl;
    adaptationManager_.adapt();
  }

protected:  
  GridPartType &gridPart_; // grid part(view), e.g. here the leaf grid the discrete space is build with 
  GridType &grid_;
  DiscreteFunctionSpaceType discreteSpace_; // discrete function space 
  DiscreteFunctionType solution_;   // the unknown 
  RestrictionProlongationType restrictProlong_ ; // local restriction/prolongation object
  AdaptationManagerType  adaptationManager_ ;    // adaptation manager handling adaptation
};

template< class FunctionSpace >
struct Function : Dune::Fem::Function< FunctionSpace, Function< FunctionSpace > >
{
  void evaluate( const typename FunctionSpace::DomainType &x,
		 typename FunctionSpace::RangeType &y ) const
  {
    y[ 0 ] = 0.0;
  }
};

// algorithm
// ---------

template <class HGridType>
double algorithm ( HGridType &grid )
{
  // we want to solve the problem on the leaf elements of the grid
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType, Dune::InteriorBorder_Partition > GridPartType;
  GridPartType gridPart(grid);

  // use a scalar function space
  typedef Dune::Fem::FunctionSpace< double, double, 
              HGridType :: dimensionworld, 1 > FunctionSpaceType;

  typedef Scheme< GridPartType, FunctionSpaceType > SchemeType;
  SchemeType scheme( gridPart );

  for( double time = 0; time <= 1.0; time += 0.05 )
  {
    if( Dune::Fem::MPIManager::rank() == 0 )
     std::cout << "time: " << time << std::endl;

    // mark element for adaptation 
    int max = 0;
    while( scheme.mark( time ) ) 
    {
      // adapt grid 
      scheme.adapt();
      max++;
      if( max > 10 )
      break;
    }
  }

  return 0.0;
}

// main
// ----

int main ( int argc, char **argv )
try
{
  // initialize MPI, if necessary
  Dune::Fem::MPIManager::initialize( argc, argv );

  // append overloaded parameters from the command line 
  Dune::Fem::Parameter::append( argc, argv );

  // append possible given parameter files 
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append( argv[ i ] );
  // make sure output parameters are added
  Dune::Fem::Parameter::append( "fem.prefix","output" );
  Dune::Fem::Parameter::append( "fem.io.savetime", "0.0" );
  Dune::Fem::Parameter::append( "fem.io.savecount", "1" );
  Dune::Fem::Parameter::append( "fem.io.outputformat", "vtk-cell" );
  Dune::Fem::Parameter::append( "fem.io.partitioning", "true" );

  // type of hierarchical grid 
  typedef Dune :: GridSelector :: GridType  HGridType ;

  // create grid from DGF file
  std::stringstream gridfilestr;
  if( HGridType :: dimension == 3 ) 
    gridfilestr << "unitcube-3d.dgf";
  else   
    gridfilestr << HGridType :: dimension << "dgrid.dgf";

  const std::string gridfile ( gridfilestr.str() );

  // the method rank and size from MPIManager are static 
  if( Dune::Fem::MPIManager::rank() == 0 )
    std::cout << "Loading macro grid: " << gridfile << std::endl;

  // construct macro using the DGF Parser 
  Dune::GridPtr< HGridType > gridPtr( gridfile );
  HGridType& grid = *gridPtr ;

  // initial load balance
  grid.loadBalance();

  // initial grid refinement
  const int level = Dune::Fem::Parameter::getValue< int >( "level", 0 );

  // number of global refinements to bisect grid width 
  const int refineStepsForHalf = Dune::DGFGridInfo< HGridType >::refineStepsForHalf();

  // refine grid 
  Dune :: Fem :: GlobalRefine::apply( grid, level * refineStepsForHalf );

  // calculate first step  
  algorithm( grid );

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}

#ifndef INTERPOLATIONSCHEME_HH
#define INTERPOLATIONSCHEME_HH

// iostream includes
#include <iostream>

// include discrete function space
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/schemes/solver.hh>

// interpolate function
#include <dune/fem/space/common/interpolate.hh>

// InterpolationScheme
//----------
template < class GridFunction, int polOrder, SolverType solver >
class InterpolationScheme
{
public:
  //! type of GridFunction
  typedef GridFunction GridFunctionType;
  typedef GridFunction ModelType;
  typedef GridFunction ExactSolutionType;

  //! grid view (e.g. leaf grid view) provided in the template argument list
  typedef typename GridFunctionType::GridPartType GridPartType;

  //! type of underyling hierarchical grid needed for data output
  typedef typename GridPartType::GridType GridType;

  //! type of function space (scalar functions, \f$ f: \Omega -> R) \f$
  typedef typename GridFunctionType :: FunctionSpaceType   FunctionSpaceType;

  //! choose type of discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;

  // choose type of discrete function, Matrix implementation and solver implementation
  typedef Solvers<DiscreteFunctionSpaceType,solver,false> UsedSolverType;
  static_assert( UsedSolverType::solverConfigured, "chosen solver is not configured" );

  typedef typename UsedSolverType::DiscreteFunctionType DiscreteFunctionType;

  static const int dimRange = FunctionSpaceType::dimRange;
  /*********************************************************/

  InterpolationScheme( GridPartType &gridPart,
             const GridFunctionType &gridFunc,
             const std::string &prefix);
  ~InterpolationScheme();
  DiscreteFunctionType &solution() { return solution_; }
  const DiscreteFunctionType &solution() const { return solution_; }
  const ExactSolutionType& exactSolution() { return gridFunc_; }

  void prepare();
  void solve( bool assemble = false );

protected:
  GridPartType  &gridPart_;         // grid part(view), e.g. here the leaf grid the discrete space is build with
  const GridFunctionType& gridFunc_;
  DiscreteFunctionSpaceType discreteSpace_; // discrete function space
  DiscreteFunctionType solution_;   // the unknown
};

// DataOutputParameters
// --------------------

template < class GridFunction, int polOrder, SolverType solver >
InterpolationScheme<GridFunction, polOrder, solver>::
InterpolationScheme( GridPartType &gridPart,
           const GridFunctionType &gridFunc,
           const std::string &prefix)
    : gridPart_( gridPart ),
      gridFunc_( gridFunc ),
      discreteSpace_( gridPart_ ),
      solution_( prefix.c_str(), discreteSpace_ )
{
  // set all DoF to zero
  solution_.clear();
}
template < class GridFunction, int polOrder, SolverType solver >
InterpolationScheme<GridFunction, polOrder, solver>::
~InterpolationScheme()
{
  std::cout << "in InterpolationScheme destructor" << std::endl;
}

  //! setup the right hand side
template < class GridFunction, int polOrder, SolverType solver >
void InterpolationScheme<GridFunction, polOrder, solver>::
prepare()
{
}

template < class GridFunction, int polOrder, SolverType solver >
void InterpolationScheme<GridFunction, polOrder, solver>::
solve ( bool assemble )
{
  interpolate( gridFunc_, solution_ );
}

#endif // end #if INTERPOLATIONSCHEME_HH

#ifndef  DUNE_FEM_P12DSPACE_TEST_POISSON_HH__
#define  DUNE_FEM_P12DSPACE_TEST_POISSON_HH__

#define VERBOSE false

// System Includes
// ---------------

#include <iostream>
#include <sstream>


// DUNE Core Includes
// ------------------
#include <dune/common/version.hh>
#include  <dune/grid/common/gridinfo.hh>

// if this value is not defined, then we have version 1.1.1
#ifndef DUNE_VERSION_HH
#define OLD_DUNE_GRID_VERSION
#endif

#include <dune/common/stdstreams.hh>
#include <dune/common/timer.hh>

#include <dgfgridtype.hh>


// DUNE-FEM includes
// -----------------

// grid parts 
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// adaptation classes 
#include <dune/fem/space/common/adaptmanager.hh>
// lagrange space 
#include <dune/fem/space/pqlagrangespace/pqlagrangespace.hh>
/*#include <dune/fem/space/lagrangespace.hh>*/

// discrete functions 
#include <dune/fem/function/adaptivefunction.hh>

#include <dune/fem/operator/discreteoperatorimp.hh>

// matrix implementations 
#include <dune/fem/operator/matrix/spmatrix.hh>

// linear solvers 
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>

// l2 norm and h1 norm 
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// Local Includes
// --------------

// the problem definition
#include "problem.hh"

// matrix assembler 
#include "laplace.hh"


namespace Dune
{

template<class GridType, unsigned int polOrder, bool isCube=false>
class Poisson_Test
  : public Test
{
public:
  /*  choose the grid partition (and hence the index set) to use
   *
   *  \note Not all index sets are continuous. The LeafIndexSet for AlbertaGrid,
   *        for example, is not. If you want to use OEM solvers, the index set
   *        must be continuous. In such a case use AdaptiveLeafGridPart.
   */
  //---- GridParts -----------------------------------------------------------
  //---- best choice ---------------------------------------------------------
  typedef LeafGridPart< GridType >                                   GridPartType;
  static const int dimension = GridType :: dimension;
  //---- other choices--------------------------------------------------------
  //typedef AdaptiveLeafGridPart< GridType, InteriorBorder_Partition > GridPartType;

  typedef FunctionSpace< double, double, dimension, 1 >              FunctionSpaceType;

  // The data functions (as defined in problem.hh)
  //---- Right Hand Side, Exact Solution, and Stiffness tensor ---------------
  typedef RHSFunction< FunctionSpaceType >                           RHSFunctionType;
  typedef ExactSolution< FunctionSpaceType >                         ExactSolutionType;
  typedef Tensor< FunctionSpaceType >                                TensorType;

  //---- Adapter for exact solution ------------------------------------------
  typedef DiscreteFunctionAdapter< ExactSolutionType, GridPartType > GridExactSolutionType;

  //---- DiscreteFunctionSpace -----------------------------------------------
  //! define the discrete function space our unkown belongs to
  typedef LagrangeSpace < FunctionSpaceType,
                          GridPartType,
                          polOrder > DiscreteSpaceType;
  /*
  typedef typename SelectType< isCube,
                               QLagrangeSpace< FunctionSpaceType,
                                               GridPartType,
                                               polOrder >,
                               PLagrangeSpace< FunctionSpaceType,
                                               GridPartType,
                                               polOrder > > :: Type  DiscreteSpaceType;
    */
/*  typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType,
 *                                         GridPartType, polOrder,
 *                                         CachingStorage >            DiscreteSpaceType;*/

  //---- DiscreteFunction ----------------------------------------------------
  //---- good choice for adaptive simulations using OEM solver ---------------
  //! define the type of discrete function we are using
  typedef AdaptiveDiscreteFunction< DiscreteSpaceType >              DiscreteFunctionType;
  //---- other possible choices, use BlockVectorDiscreteFunction for ISTL ----
  //typedef BlockVectorDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;
  //typedef ManagedDiscreteFunction< VectorDiscreteFunction< DiscreteSpaceType, DynamicVector< double > > > DiscreteFunctionType;

  //---- MatrixObjects -------------------------------------------------------
  //---- good choice for build in solvers ------------------------------------
  //! define the type of the system matrix object
  typedef SparseRowMatrixTraits< DiscreteSpaceType,
                                 DiscreteSpaceType >                 MatrixObjectTraits;

  //---- other choices, ISTLMatrixTraits for BCRSMatrix from DUNE-ISTL -------
  //typedef ISTLMatrixTraits < DiscreteSpaceType, DiscreteSpaceType > MatrixObjectTraits;
  //typedef BlockMatrixTraits < DiscreteSpaceType, DiscreteSpaceType > MatrixObjectTraits;
  //typedef OnTheFlyMatrixTraits < DiscreteSpaceType, DiscreteSpaceType > MatrixObjectTraits;

  //! define the discrete laplace operator, see ./fem.cc
  typedef LaplaceFEOp< DiscreteFunctionType, MatrixObjectTraits,
                       TensorType >                                  LaplaceOperatorType;

  //---- InverseOperator ----------------------------------------------------
  //---- good choice for build in CG solver ---------------------------------
  //! define the inverse operator we are using to solve the system 
  typedef CGInverseOp< DiscreteFunctionType, LaplaceOperatorType >   InverseOperatorType;
    //---- other choices ------------------------------------------------------ 
    //typedef ISTLBICGSTABOp< DiscreteFunctionType, LaplaceOperatorType >
    //typedef ISTLCGOp< DiscreteFunctionType, LaplaceOperatorType >
    //typedef OEMCGOp<DiscreteFunctionType,LaplaceOperatorType>
    //typedef OEMBICGSTABOp<DiscreteFunctionType,LaplaceOperatorType>
    //typedef OEMBICGSQOp<DiscreteFunctionType,LaplaceOperatorType>
    //typedef OEMGMRESOp<DiscreteFunctionType,LaplaceOperatorType>
    //    InverseOperatorType;

private:
  //! set the dirichlet points to exact values
  template< class EntityType, class GridFunctionType, class DiscreteFunctionType >
    void boundaryTreatment( const EntityType &entity,
                            const GridFunctionType &exactSolution,
                            DiscreteFunctionType &rhs,
                            DiscreteFunctionType &solution )
    {
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

      typedef typename DiscreteFunctionSpaceType::LocalFiniteElementType
        LocalFiniteElementType;
      typedef typename GridFunctionType::LocalFunctionType LocalExactSolutionType;

/*      typedef typename DiscreteFunctionSpaceType::LagrangePointSetType LagrangePointSetType;*/
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

/*      const int faceCodim = 1;*/
/*      typedef typename GridPartType :: IntersectionIteratorType      IntersectionIteratorType;*/
/*      typedef typename LagrangePointSetType
 *                :: template Codim< faceCodim >
 *                :: SubEntityIteratorType                             FaceDofIteratorType;*/

      if( ! entity.hasBoundaryIntersections() )
        return;

      const DiscreteFunctionSpaceType &dfSpace = rhs.space();
      const LocalFiniteElementType & fem = dfSpace.localFiniteElement();

      LocalExactSolutionType exactLocal = exactSolution.localFunction( entity );
      LocalFunctionType rhsLocal        = rhs.localFunction( entity );
      LocalFunctionType solutionLocal   = solution.localFunction( entity );

      std::vector< typename LocalExactSolutionType :: RangeType > out;
      fem.localInterpolation().interpolate(exactLocal, out);

      // This is a hack: it returns a vector of ones and zeros where a one
      // marks a boundary dof
      // Note: It is not possible to use the intersection's method boundaryID()
      // with this concept, is it?
      std::vector< typename LocalExactSolutionType :: RangeType > dofOut;
      BoundaryCheck<EntityType> bc(entity);
      fem.localInterpolation().interpolate(bc, dofOut);

      for( size_t i = 0 ; i < dofOut.size() ; ++i )
      {
        if( dofOut[i][0] == 1 )
        {
          rhsLocal[i]      = out[i][0];
          solutionLocal[i] = out[i][0];
        }
      }
/*      const LagrangePointSetType &lagrangePointSet = dfSpace.lagrangePointSet( entity );*/

/*      IntersectionIteratorType it          = gridPart.ibegin( entity );
 *      const IntersectionIteratorType endit = gridPart.iend( entity );
 *      for( ; it != endit; ++it )
 *      {
 *        if( !it->boundary() )
 *          continue;

 *        const int face = it->numberInSelf();
 *        FaceDofIteratorType faceIt
 *          = lagrangePointSet.template beginSubEntity< faceCodim >( face );
 *        const FaceDofIteratorType faceEndIt
 *          = lagrangePointSet.template endSubEntity< faceCodim >( face );
 *        for( ; faceIt != faceEndIt; ++faceIt )
 *        {
 *          const int localDof = *faceIt;

 *          typename LocalExactSolutionType :: RangeType phi;
 *          exactLocal.evaluate( lagrangePointSet.point( localDof ), phi );

 *          rhsLocal[ localDof ]      = phi[ 0 ];
 *          solutionLocal[ localDof ] = phi[ 0 ];
 *        }
 *      }*/
    }

private:
  //! solve the resulting linear system
  void solve ( LaplaceOperatorType &laplace,
               const DiscreteFunctionType &rhs,
               DiscreteFunctionType &solution )
  {
    // solve the linear system (with CG)
    double dummy = 12345.67890;
    double solverEps = 1e-8 ;

    // create inverse operator 
    InverseOperatorType cg( laplace, dummy, solverEps, 20000, VERBOSE );

    // solve the system 
    cg( rhs, solution );
  }

public:
  explicit Poisson_Test ( const std::string &gridFileName )
  : gridPtr_(gridFileName),
    grid_(*gridPtr_),
    prevErrors_(),
    eocState_( true )
  {}

  //! setup of the linear system, solution contains values from previous runs
  // (in case of adaptive scheme)
  void algorithm ( const unsigned int repeat )
  {
    std::cout << "starting algorithm run no. " << repeat << std::endl;

    GridPartType gridPart( grid_ );
    DiscreteSpaceType space( gridPart );
    RHSFunctionType exactRhs;
    DiscreteFunctionType u( "solution", space );
    u.clear();

    std::cout << "size of u: " << u.size() << std::endl;

    // type of grid iterator 
    typedef typename DiscreteSpaceType :: IteratorType               IteratorType;
    // type of entity 
    typedef typename IteratorType :: Entity                          EntityType;

    // create exact solution
    ExactSolutionType uexact; 

    // create adapter (for visualization with grape)
    GridExactSolutionType ugrid( "exact solution", uexact, gridPart,
                                 DiscreteSpaceType :: polynomialOrder + 1 );

    // create laplace assembler (is assembled on call of systemMatrix by solver)
    LaplaceOperatorType laplace(space);
    laplace.assemble();

    // create right hand side
    DiscreteFunctionType rhs( "rhs", space );

    // initialize as zero 
    rhs.clear();

    // setup right hand side 
    RightHandSideAssembler< DiscreteFunctionType >
      :: template assemble< 2 * DiscreteSpaceType :: polynomialOrder > ( exactRhs , rhs );

    // we're having some trouble with NaNs lately, so check the right hand side for NaNs
    if( !rhs.dofsValid() )
      std :: cout << "right hand side invalid before boundary treatment." << std :: endl;

    // set Dirichlet Boundary to exact solution
    IteratorType endit = space.end();
    for( IteratorType it = space.begin(); it != endit; ++it )
    {
      const EntityType &entity = *it;
      // in entity has intersections with the boundary adjust dirichlet nodes 
      if( entity.hasBoundaryIntersections() )
        boundaryTreatment( entity, ugrid, rhs, u );
/*        boundaryTreatment( entity, uexact, rhs, u );*/
    }

    // check the right hand side for NaNs again
    if( ! rhs.dofsValid() )
      std :: cout << "right hand side invalid after boundary treatment." << std :: endl;

    // solve the linear system
    solve( laplace, rhs, u );

    std :: vector<double> errors;
    // calculate L2 - Norm 
    L2Norm< GridPartType > l2norm( gridPart );
    errors.push_back( l2norm.distance( ugrid, u ) );
    // calculate H1 - Norm 
    H1Norm< GridPartType > h1norm( gridPart );
    errors.push_back( h1norm.distance( ugrid, u ) );

    // in verbose mode print errors  
    std :: cout << "L2-Error: " << errors[0]  << "\n\n";
    std :: cout << "H1-Error: " << errors[1]  << "\n" << std :: endl;

    if( prevErrors_.size() != 0 ) 
    {
      double l2eoc = log(prevErrors_[0]/errors[0]) / M_LN2;
      double h1eoc = log(prevErrors_[1]/errors[1]) / M_LN2;
      std :: cout << "L2 - EOC: " << l2eoc << "\n\n";
      std :: cout << "H1 - EOC: " << h1eoc << "\n\n";
      if(l2eoc < 0.8 + polOrder || h1eoc < -0.2 + polOrder)
        eocState_ = false;
    }
    prevErrors_ = errors;
  }

  bool eocState()
  {
    return eocState_;
  }

  void run()
  {
    std::cout << std::endl << std::endl;
    std::cout << "=============================================================" << std::endl;
    std::cout << "Starting Poisson test with EOC measurements for" << std::endl;
    std::cout << "base functions of polynomial order " << polOrder << std::endl;
    std::cout << "on a grid of dimension             " << GridType::dimension << std::endl;
    std::cout << "=============================================================" << std::endl;

    const int startLevel = 0;

    // refine the grid until the startLevel is reached
    grid_.globalRefine(startLevel);

    // set some parameters
    const int eocSteps   = 3;

    algorithm( 0 );

    for( int eocloop = 1; eocloop < eocSteps; ++eocloop )
    {
      // refine the grid
      grid_.globalRefine( DGFGridInfo<GridType>::refineStepsForHalf() );

      algorithm( eocloop );
    } /***** END of EOC Loop *****/
    _test(eocState_);
  }

private:
  GridPtr<GridType>     gridPtr_;
  GridType             &grid_;
  std::vector<double>   prevErrors_;
  bool                  eocState_;
};




/*// main programm, run algorithm twice to calc EOC
 *int main( int argc, char **argv )
 *{
 *  // initialize MPI 
 *  MPIManager :: initialize ( argc, argv );
 *  const int rank = MPIManager :: rank ();
 *  
 *  try
 *  {
 *    // append parameters from the comand line 
 *    Parameter::append( argc, argv );

 *    // append parameters from the parameter file  
 *    Parameter::append( (argc < 2) ? "parameter" : argv[ 1 ] );

 *    // generate GridPointer holding grid instance
 *    GridPtr< GridType > gridptr = initialize(std::string("Poisson problem"));

 *    // get grid reference 
 *    GridType& grid = *gridptr ;

 *    Stepper stepper(grid);
 *    compute<true>(stepper);

 *    Parameter :: write("parameter.log");

 *    return 0;
 *  }
 *  catch( Exception &exception )
 *  {
 *    if( rank == 0 )
 *      std :: cerr << exception << std :: endl;
 *    return 1;
 *  }
 *}*/

} // end namespace Dune

#endif  /*DUNE_FEM_P12DSPACE_TEST_POISSON_HH__*/

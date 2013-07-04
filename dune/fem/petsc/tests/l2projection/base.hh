// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef FEMHOWTO_BASE_HH
#define FEMHOWTO_BASE_HH

#if defined HAVE_PETSC

#include <sstream>

#include <dune/common/version.hh>

// DGF gridtype 
#if ! DUNE_VERSION_NEWER_REV(DUNE_GRID,2,1,0)
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#endif

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/femtimer.hh>
#include <dune/fem/misc/femeoc.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem/misc/petsc/petsccommon.hh>


/* #########################################
 * for debugging
 */
template< typename DF >
void viewPETScVector ( const DF &dFunction, const std::string beforeMsg = "", const std::string afterMsg = "" )
{
  if( dFunction.space().grid().comm().rank() == 0 )
  {
    std::cout << beforeMsg << std::endl;
  }
  VecView( *dFunction.petscVector(), PETSC_VIEWER_STDOUT_WORLD );
  if( dFunction.space().grid().comm().rank() == 0 )
  {
    std::cout << afterMsg << std::endl;
  }
}
/*
 * ########################################
 */


template< class HGridType > /*@LST0S@*/
Dune::GridPtr< HGridType > initialize ( const std::string &problemDescription )
{
  // ----- read in runtime parameters ------
  const std::string filekey = Dune::Fem::IOInterface::defaultGridKey( HGridType::dimension );
  const std::string filename = Dune::Fem::Parameter::getValue< std::string >( filekey ); /*@\label{base:param0}@*/

  // initialize grid
  Dune::GridPtr< HGridType > gridptr(filename);
  Dune::Fem::Parameter::appendDGF( filename );

  // output of error and eoc information
  std::string eocOutPath = Dune::Fem::Parameter::getValue<std::string>("femhowto.eocOutputPath",  /*@\label{base:param1}@*/
                                                            std::string("."));
  std::string eocFile = eocOutPath + std::string("/eoc.tex");
  Dune::Fem::FemEoc::initialize(eocOutPath, "eoc", problemDescription); /*@\label{base:femeocInit}@*/

  // and refine the grid until the startLevel is reached
  const int startLevel = Dune::Fem::Parameter::getValue<int>("femhowto.startLevel", 0);
  for(int level=0; level < startLevel ; ++level)
    Dune::Fem::GlobalRefine::apply(*gridptr, 1 ); /*@\label{base:globalRefine1}@*/
  return gridptr;
} /*@LST0E@*/

/////////////////////////////////////////////////////////////////////////////
//
//  compute algorithm 
//
/////////////////////////////////////////////////////////////////////////////
template <class Algorithm> /*@LST0S@*/
void compute(Algorithm& algorithm)
{
  typedef typename Algorithm::DiscreteFunctionType DiscreteFunctionType;
  typename Algorithm::DiscreteSpaceType& space = algorithm.space();
  typename Algorithm::GridPartType& gridPart = space.gridPart();
  typedef typename Algorithm::GridPartType::GridType HGridType;
  HGridType& grid = gridPart.grid();


  // get some parameters
  const int eocSteps   = Dune::Fem::Parameter::getValue<int>("femhowto.eocSteps", 1);

  // Initialize the DataOutput that writes the solution on the harddisk in a
  // format readable by e.g. Paraview
  // in each loop for the eoc computation the results at
  // the final time is stored
  typedef Dune::tuple< DiscreteFunctionType* > IOTupleType;
  typedef Dune::Fem::DataOutput<HGridType, IOTupleType> DataOutputType;

  const unsigned int femTimerId = Dune::FemTimer::addTo("timestep");
  for(int eocloop=0; eocloop < eocSteps; ++eocloop)
  {
    /*
     * We want to re-create the discrete function in each EOC step here, because
     * the PETScDiscreteFunction is not yet meant to support adaption
     */

    // solution function
    DiscreteFunctionType u("solution",space);
    IOTupleType dataTup ( &u );
    DataOutputType dataOutput( grid, dataTup );

    /*
     * continue with the computations
     */
    Dune::FemTimer :: start(femTimerId);

    // do one step
    algorithm(u); /*@\label{base:solveProblem}@*/

    double runTime = Dune::FemTimer::stop(femTimerId);

    Dune::FemTimer::printFile("./timer.out");
    Dune::FemTimer::reset(femTimerId);

    // Write solution to hd
    dataOutput.writeData(eocloop);

    // finalize algorithm (compute errors)
    algorithm.finalize(u); /*@\label{base:finalize}@*/

    // calculate grid width
    const double h = Dune::Fem::GridWidth::calcGridWidth(gridPart);

    if( Dune::Fem::Parameter :: verbose() )
      Dune::Fem::FemEoc::write(h,grid.size(0),runTime,0, std::cout);
    else
      Dune::Fem::FemEoc::write(h,grid.size(0),runTime,0);

    // Refine the grid for the next EOC Step. If the scheme uses adaptation,
    // the refinement level needs to be set in the algorithms' initialize method.
    if(eocloop < eocSteps-1)
    {
      Dune::Fem::GlobalRefine::apply(grid,Dune::DGFGridInfo<HGridType>::refineStepsForHalf()); /*@\label{base:globalRefine2}@*/
      grid.loadBalance();
    }
  } /***** END of EOC Loop *****/
  Dune::FemTimer::removeAll();
}/*@LST0E@*/


#endif // #if defined HAVE_PETSC

#endif // #ifndef FEMHOWTO_BASE_HH

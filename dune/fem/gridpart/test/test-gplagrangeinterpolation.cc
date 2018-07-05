#if defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX
#if GRIDDIM == 3
#define COMPILE_TEST
#endif
#endif

#include <config.h>
#include <iostream>
#include <sstream>
#include <string>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/gridpart/filter/radialfilter.hh>
#include <dune/fem/gridpart/filter/basicfilterwrapper.hh>

#include <dune/fem/misc/double.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/io/streams/virtualstreams.hh>
#include <dune/fem/io/file/vtkio.hh>

#include <dune/fem/test/exactsolution.hh>

using namespace Dune;

#ifdef COMPILE_TEST

// polynom approximation order of quadratures,
// at least poolynom order of basis functions
#ifdef POLORDER
  constexpr int polOrder = POLORDER;
#else
  constexpr int polOrder = 1;
#endif

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#endif
  typedef Dune :: GridSelector :: GridType MyGridType;

  typedef Fem::LeafGridPart< MyGridType > HostGridPartType;
  typedef Fem::RadialFilter< MyGridType::ctype, MyGridType::dimensionworld > BasicFilterType;
  typedef Fem::BasicFilterWrapper< HostGridPartType, BasicFilterType > FilterType;
  typedef Fem::FilteredGridPart< HostGridPartType, FilterType, true > GridPartType;

  typedef Fem::FunctionSpace< double, double, MyGridType::dimensionworld, 1 > FunctionSpaceType;

  typedef Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;

  typedef Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

  typedef Fem::ExactSolution< FunctionSpaceType > ExactSolutionType;


  void writeOut ( Fem::VirtualOutStream out, const DiscreteFunctionType &solution )
  {
    out << solution;
    out.flush();
  }

  void readBack ( Fem::VirtualInStream in, DiscreteFunctionType &solution )
  {
    solution.clear();
    in >> solution;
  }

  template <class HGridType>
  class TestGrid
  {
    typedef TestGrid<HGridType> ThisType;

  protected:
    TestGrid ()
    : gridptr_( macroGridName() )
    {
      gridptr_->loadBalance();
    }

  public:
    TestGrid ( const ThisType & ) = delete;

    ThisType &operator= ( const ThisType & ) = delete;

    static ThisType &instance ()
    {
      static ThisType staticInstance;
      return staticInstance;
    }

    static HGridType &grid ()
    {
      return *(instance().gridptr_);
    }

    static int refineStepsForHalf ()
    {
      return DGFGridInfo< HGridType >::refineStepsForHalf();
    }

  protected:
    static std::string macroGridName ()
    {
      std::ostringstream s;
      s << HGridType::dimension << "dgrid.dgf";
      return s.str();
    }

    GridPtr< HGridType > gridptr_;
  };


  int main(int argc, char ** argv)
  {
    Dune::Fem::MPIManager :: initialize( argc, argv );
    try
    {
      auto& grid = TestGrid<MyGridType> :: grid();
      const int step = TestGrid<MyGridType> :: refineStepsForHalf();
      grid.globalRefine( 2*step );
      HostGridPartType hostGridPart (grid );
      BasicFilterType::GlobalCoordinateType center( 0 );
      BasicFilterType basicFilter( center, .25 );
      FilterType filter( hostGridPart, basicFilter );
      GridPartType gridPart( hostGridPart, filter );

      DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
      ExactSolutionType f;
      DiscreteFunctionType solution( "solution", discreteFunctionSpace );
      solution.clear();

      std::cout << "maxDofs = " << discreteFunctionSpace.blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize << std::endl;

      // perform Lagrange interpolation
      interpolate( gridFunctionAdapter( f, gridPart, discreteFunctionSpace.order() + 2 ), solution );
      solution.communicate();

      // output to vtk file
      Fem::VTKIO<GridPartType> vtkWriter(gridPart);
      vtkWriter.addVertexData(solution);
      vtkWriter.pwrite( "vtxprojection", Fem::Parameter::commonOutputPath().c_str(), "", Dune::VTK::ascii );
      return 0;
    }
    catch( Exception e )
    {
      std :: cerr << e.what() << std :: endl;
      return 1;
    }
  }
#else
  // no ALUGrid, no test
  int main(int argc, char ** argv)
  {
    return 0;
  }
#endif

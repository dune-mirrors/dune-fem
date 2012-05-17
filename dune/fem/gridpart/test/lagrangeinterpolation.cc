#if defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX 
#if GRIDDIM == 3 
#define COMPILE_TEST 
#endif
#endif

#include <config.h>
#include <iostream>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/gridpart/filter/radialfilter.hh>
#include <dune/fem/gridpart/filter/basicfilterwrapper.hh>

#include <dune/fem/misc/double.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/io/streams/xdrstreams.hh>
#include <dune/fem/io/streams/asciistreams.hh>
#include <dune/fem/io/streams/virtualstreams.hh>
#include <dune/fem/io/file/vtkio.hh>

#include <dune/fem/test/exactsolution.hh>

using namespace Dune;

// polynom approximation order of quadratures, 
// at least poolynom order of basis functions 
#ifdef POLORDER
  const int polOrder = POLORDER;
#else
  const int polOrder = 1;
#endif

#ifdef COMPILE_TEST
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
  typedef Dune :: GridSelector :: GridType MyGridType;

  // typedef AdaptiveLeafGridPart< MyGridType > HostGridPartType;
  typedef LeafGridPart< MyGridType > HostGridPartType;
  typedef Fem::RadialFilter< MyGridType::ctype, MyGridType::dimensionworld > BasicFilterType;
  typedef Fem::BasicFilterWrapper< HostGridPartType, BasicFilterType > FilterType;
  typedef Fem::FilteredGridPart< HostGridPartType, FilterType, true > GridPartType;
  // typedef AdaptiveLeafGridPart< MyGridType > GridPartType;

  typedef FunctionSpace< double, double, MyGridType::dimensionworld, 1 > FunctionSpaceType;

  typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
    DiscreteFunctionSpaceType;

  typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

  typedef ExactSolution< FunctionSpaceType > ExactSolutionType;



  void writeOut ( VirtualOutStream out, const DiscreteFunctionType &solution )
  {
    out << solution;
    out.flush();
  }

  void readBack ( VirtualInStream in, DiscreteFunctionType &solution )
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

  private:
    TestGrid ( const ThisType & );

    ThisType &operator= ( const ThisType & );

  public:
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
    MPIManager :: initialize( argc, argv );
    try
    {
      MyGridType &grid = TestGrid<MyGridType> :: grid();
      const int step = TestGrid<MyGridType> :: refineStepsForHalf();
      grid.globalRefine( 2*step );
      HostGridPartType hostGridPart (grid );
      BasicFilterType::GlobalCoordinateType center( 0 );
      BasicFilterType basicFilter( center, .25 );
      FilterType filter( hostGridPart, basicFilter );
      GridPartType gridPart( hostGridPart, filter );
      // GridPartType gridPart ( grid );
      //
      
      DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
      ExactSolutionType f;
      DiscreteFunctionType solution( "solution", discreteFunctionSpace );
      solution.clear();

      std::cout << "codim c=" << 0 << " : " << gridPart.indexSet().size(0) << " " 
                << hostGridPart.indexSet().size(0) << std::endl;
      std::cout << "codim c=" << 3 << " : " << gridPart.indexSet().size(3) << " " 
                << hostGridPart.indexSet().size(3) << std::endl;

      std :: cout << "maxDofs = " << discreteFunctionSpace.mapper().maxNumDofs() << std :: endl;

      //! perform Lagrange interpolation
      LagrangeInterpolation< DiscreteFunctionType >
        :: interpolateFunction( f, solution );
      solution.communicate();

      // output to vtk file
      VTKIO<GridPartType> vtkWriter(gridPart);
      vtkWriter.addVertexData(solution);
      vtkWriter.pwrite("vtxprojection",
                        Parameter::commonOutputPath().c_str(),"",
                        Dune::VTKOptions::ascii);
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

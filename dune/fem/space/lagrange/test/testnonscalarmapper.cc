#include <config.h>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/io/file/vtkio.hh>

#ifndef POLORDER
const int PolynomialOrder = 2;
#else
const int PolynomialOrder = POLORDER;
#endif

#if SCALAR
const int dimRange = 1;
#else
const int dimRange = Dune :: GridSelector :: GridType :: dimensionworld;
#endif


using namespace Dune;
using namespace Fem;

// Main Program
// ------------
int main( int argc, char **argv )
{
  // initialize MPI
  Dune::Fem::MPIManager :: initialize ( argc, argv );
  const int rank = MPIManager :: rank ();

  try
  {
    const int minLevel = (argc > 2) ? atoi(argv[1]) : 0 ;
    const int times = (argc >3 )? atoi(argv[2]) : 4 ;

    typedef GridSelector :: GridType GridType;

    std::stringstream gridFile;
    gridFile << GridType::dimensionworld <<"dgrid.dgf";

    // generate GridPointer holding grid instance
    GridPtr< GridType > gridptr ( gridFile.str() );

    // get grid reference
    GridType& grid = *gridptr ;

    GlobalRefine::apply(grid, minLevel);

    typedef AdaptiveLeafGridPart<GridType> GridPartType;
    GridPartType gridPart (grid);

    typedef FunctionSpace< double, double, GridType::dimensionworld, dimRange  > FunctionSpace;

    typedef FunctionSpace::DomainType DomainType;
    typedef FunctionSpace::RangeType RangeType;

    typedef LagrangeDiscreteFunctionSpace < FunctionSpace, GridPartType, PolynomialOrder, CachingStorage > DiscreteFunctionSpace;
    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpace > DiscreteFunction;

    DiscreteFunctionSpace dFspace( gridPart );
    DiscreteFunction solution("testSolution", dFspace);

    const int faceCodim = 0;

#if SCALAR
    DomainType v;
    v[0] = M_PI;
    v[1] =-std::sqrt(2);
    v/= v.two_norm();

    Parameter :: get( "test.velocity", v, v);
#endif

    for(int i= 0;i<times;++i)
    {
      if(i>0)
      {
        std::cout<<"refining ... "<<std::endl;
        GlobalRefine::apply(grid,DGFGridInfo<GridType>::refineStepsForHalf());
      }

      solution.clear();

      for( const auto& entity : dFspace )
      {
        auto solutionLocal = solution.localFunction( entity );

        const auto& geo = entity.geometry();

        const auto& lagrangePointSet = dFspace.lagrangePointSet( entity );

        const int face =0;

        auto faceIt = lagrangePointSet.beginSubEntity< faceCodim >( face );
        const auto faceEndIt = lagrangePointSet.endSubEntity< faceCodim >( face  );

        for( ; faceIt != faceEndIt; ++faceIt )
        {
          const auto& local = lagrangePointSet.point( *faceIt );
          const auto global = geo.global( local );

          for(int r=0;r<dimRange;++r)
          {
            const int localDof = (*faceIt)*dimRange + r;
            RangeType phi, phiLocal;
#if SCALAR
            phi = global * v;
#else
            phi = global;
#endif
            const double value = phi[ i ];
            solutionLocal[localDof] = value;
            }
          }
      }

      SubsamplingVTKIO<GridPartType> vtkio(gridPart, 1);
      vtkio.addVertexData(solution);
      std::stringstream name;
      name<<"testnonscalar_"<<i;
      vtkio.write(name.str().c_str());
    }

    return 0;
  }
  catch( const Exception &exception )
  {
    if( rank == 0 )
      std :: cerr << exception << std :: endl;
    return 1;
  }
}

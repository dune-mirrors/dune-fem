#include <config.h>

//#define USE_TWISTFREE_MAPPER

#include <dune/common/fvector.hh>

// DGF gridtype 
#include <dgfgridtype.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
// adaptation classes 
#include <dune/fem/space/common/adaptmanager.hh>
// lagrange space 
#include <dune/fem/space/lagrangespace.hh>
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
const int dimRange = GridType :: dimensionworld;
#endif


using namespace Dune;

// Main Program
// ------------
int main( int argc, char **argv )
{
  // initialize MPI 
  MPIManager :: initialize ( argc, argv );
  const int rank = MPIManager :: rank ();
  
  try
  {
    if( argc < 2 )
    {
      std::cout << "Usuage: " << argv[ 0 ] << " <startLevel> <times>" << std::endl;
      return 1;
    }

    // generate GridPointer holding grid instance
    GridPtr< GridType > gridptr ("2dgrid.dgf");

    // get grid reference 
    GridType& grid = *gridptr ;

    const int minLevel = atoi(argv[1]);
    const int times = atoi(argv[2]);

    GlobalRefine::apply(grid, minLevel);
//    grid.globalRefine(minLevel);
                  

    typedef AdaptiveLeafGridPart<GridType> GridPartType;
    GridPartType gridPart (grid);

    typedef FunctionSpace< double, double, GridType::dimensionworld, dimRange  > FunctionSpace;

    typedef FunctionSpace::DomainType DomainType;
    typedef FunctionSpace::RangeType RangeType;
    typedef FunctionSpace::RangeFieldType RangeFieldType;


    typedef LagrangeDiscreteFunctionSpace < FunctionSpace, GridPartType, PolynomialOrder, CachingStorage > DiscreteFunctionSpace;
    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpace > DiscreteFunction;

    DiscreteFunctionSpace dFspace( gridPart );
    DiscreteFunction solution("testSolution", dFspace);

    typedef DiscreteFunctionSpace :: IteratorType IteratorType;
    typedef IteratorType :: Entity EntityType;

    typedef DiscreteFunctionSpace :: LagrangePointSetType LagrangePointSetType;

    const int faceCodim = 0;
    typedef LagrangePointSetType :: Codim< faceCodim > :: SubEntityIteratorType
      FaceDofIteratorType;

    typedef GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef IntersectionIteratorType::Intersection IntersectionType; 
        
#if SCALAR
    DomainType v;
    v[0] = M_PI;
    v[1] =-std::sqrt(2);
    v/= v.two_norm();

    Dune :: Parameter :: get( "test.velocity", v, v);
#endif



    for(int i= 0;i<times;++i)
    {
      if(i>0)
      {
        std::cout<<"refining ... "<<std::endl;
        GlobalRefine::apply(grid,DGFGridInfo<GridType>::refineStepsForHalf());
      }

      solution.clear();

      const IteratorType endit = dFspace.end();
      for( IteratorType it = dFspace.begin(); it != endit; ++it )
      {
        const EntityType &entity = *it;
        typedef DiscreteFunction :: LocalFunctionType LocalFunctionType;

        LocalFunctionType solutionLocal = solution.localFunction( entity );      

        typedef EntityType :: Geometry Geometry;
        const Geometry& geo = entity.geometry();

        const LagrangePointSetType &lagrangePointSet = dFspace.lagrangePointSet( entity );

#if 0
        IntersectionIteratorType iIt = gridPart.ibegin( entity );
        const IntersectionIteratorType endiIt = gridPart.iend( entity );

        for( ; iIt != endiIt; ++iIt )
        {                  
          const IntersectionType &intersection = *iIt;
          const int face = intersection.indexInInside();
#endif
          const int face =0;


          FaceDofIteratorType faceIt
            = lagrangePointSet.beginSubEntity< faceCodim >( face );
          const FaceDofIteratorType faceEndIt
            = lagrangePointSet.endSubEntity< faceCodim >( face  );

          for( ; faceIt != faceEndIt; ++faceIt )
          {
            const DomainType &local = lagrangePointSet.point( *faceIt );
            const DomainType global = geo.global( local );

            for(int r=0;r<dimRange;++r)
            {
              const int localDof = (*faceIt)*dimRange + r;
              RangeType phi, phiLocal;
#if SCALAR
              phi = global * v;
#else
              phi = global;
#endif

              solutionLocal.baseFunctionSet().evaluate( localDof, local, phiLocal );
              const double value = ( phiLocal * phi ) / phiLocal.two_norm();

              solutionLocal[localDof] = value;
            }
          }
//      }
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

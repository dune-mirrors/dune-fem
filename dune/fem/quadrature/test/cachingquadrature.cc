#include <config.h>
#include <iostream>

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#include<dune/grid/common/gridinfo.hh>

#include <dune/fem/gridpart/gridpart.hh>

#include <dune/fem/quadrature/elementquadrature.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem/io/parameter.hh>

using namespace Dune;

template<class GridPartType, bool caching>
struct QuadratureChooser;

template<class GridPartType>
struct QuadratureChooser<GridPartType,false>
{
  typedef ElementQuadrature< GridPartType, 1 > Quadrature;
};
template<class GridPartType>
struct QuadratureChooser<GridPartType,true>
{
  typedef CachingQuadrature< GridPartType, 1 > Quadrature;
};

template<class GridPartType, bool caching>
class TestCaching
{
private:
  // types of iterators, entities and intersection
  typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;
  typedef typename IteratorType :: Entity EntityType;
  typedef typename EntityType :: EntityPointer EntityPointerType;
  typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType :: Intersection IntersectionType;

  typedef typename QuadratureChooser<GridPartType,caching>::Quadrature FaceQuadratureType; 

  GridPartType & gridPart_;
  int order_;
  const double eps_;

public:

  TestCaching(GridPartType & gridPart, int order, const double eps = 1e-8) 
    : gridPart_(gridPart), order_(order), eps_(eps) 
  { }

  void runTest()
  {
    IteratorType endit = gridPart_.template end<0>();
    for (IteratorType it = gridPart_.template begin<0>(); it != endit; ++it) 
    {
      EntityType & entity = *it;
      const IntersectionIteratorType iend = gridPart_.iend(entity);
      for( IntersectionIteratorType iit = gridPart_.ibegin(entity); iit != iend; ++iit )
      {
        const IntersectionType &intersection = *iit;
        if (intersection.neighbor())
        {
          EntityPointerType epInside = intersection.inside();
          EntityType & inside = *epInside;
          FaceQuadratureType faceQuadInner(gridPart_, intersection, order_,
              FaceQuadratureType::INSIDE);

          EntityPointerType epOutside = intersection.outside();
          EntityType & outside = *epOutside;
          FaceQuadratureType faceQuadOuter(gridPart_, intersection, order_,
              FaceQuadratureType::OUTSIDE);

          const int faceQuadInner_nop = faceQuadInner.nop();
          for (int qp = 0; qp < faceQuadInner_nop; ++qp)
          {
            Dune::FieldVector<typename GridType::ctype,GridType::dimensionworld> globalInside, globalOutside;
            globalInside = inside.geometry().global(faceQuadInner.point(qp));
            globalOutside = outside.geometry().global(faceQuadOuter.point(qp));
            if( (globalInside-globalOutside).two_norm() > eps_) 
            {
              std::cout << " Error: x(inside) = " << globalInside << " != " 
                        << globalOutside << " = x(outside) "
                        << " at intersection.indexInInside() = " << intersection.indexInInside()  
                        << ", intersection.indexInOutside() = " << intersection.indexInOutside() << std::endl;
              break;
            }
          }
        }
      }
    }
  }
};

// main program 
int main(int argc, char ** argv) 
{
  try
  {
    MPIManager :: initialize( argc, argv );
    // *** Initialization
    Parameter::append(argc, argv);
    if (argc==2) Parameter::append(argv[1]);
    else Parameter::append("parameter");

    std::string filename;
    Parameter::get("fem.gridfile", filename);
    Dune::GridPtr< GridType > gridptr( filename );
    GridType &grid = *gridptr;

    int startlevel = Parameter::getValue<int>("fem.startlevel");
    grid.globalRefine( startlevel );
    Dune::gridinfo(grid);

    int quadOrder = Parameter::getValue<int>("fem.quadorder");

    typedef LeafGridPart< GridType > GridPartType;
    GridPartType gridPart( grid );

    const double eps = 1e-8;
 
    {
      std::cout << "Testing ElementQuadratures: " << std::endl;
      TestCaching<GridPartType,false> testCaching(gridPart, quadOrder, eps);
      testCaching.runTest();
    }
    {
      std::cout << "Testing CachingQuadratures: " << std::endl;
      TestCaching<GridPartType,true> testCaching(gridPart, quadOrder, eps);
      testCaching.runTest();
    }
    return 0;
  }
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}

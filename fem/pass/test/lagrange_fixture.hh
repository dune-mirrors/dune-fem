#ifndef DUNE_LAGRANGE_FIXTURE_HH
#define DUNE_LAGRANGE_FIXTURE_HH

//- Dune includes
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

namespace Dune {
  
  template <class GridImp, int polOrd>
  class Lagrange_Fixture {
  public:
    typedef GridImp GridType;
    typedef DofManager<GridType> DofManagerType;
    typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;
    typedef LeafGridPart<GridType> GridPartType;
    typedef FunctionSpace<double, double, GridType :: dimension , 1> FunctionSpaceType;
    typedef LagrangeDiscreteFunctionSpace<
      FunctionSpaceType, GridPartType, polOrd> DiscreteFunctionSpaceType;
    
  public:
    Lagrange_Fixture(GridType& grid) :
      grid_(grid),
      dm_(DofManagerFactoryType::getDofManager(grid_)),
      part_(grid_),
      spc_(part_)
    {}

    GridType& grid() { return grid_; }
    DiscreteFunctionSpaceType& space() { return spc_; }
  private:
    GridType& grid_;
    DofManagerType& dm_;
    GridPartType part_;
    DiscreteFunctionSpaceType spc_;
  };
  
} // end namespace Dune
#endif

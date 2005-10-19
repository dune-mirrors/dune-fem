#ifndef DUNE_LAGRANGE_FIXTURE_HH
#define DUNE_LAGRANGE_FIXTURE_HH

//- Dune includes
#include <dune/fem/dofmanager.hh>
#include <dune/fem/lagrangebase.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/alu3dgrid.hh>

namespace Dune {
  
  template <int polOrd>
  class Lagrange_Fixture {
  public:
    typedef ALU3dGrid<3, 3, tetra> GridType;
    typedef DofManager<GridType> DofManagerType;
    typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;
    typedef LeafGridPart<GridType> GridPartType;
    typedef FunctionSpace<double, double, 3, 1> FunctionSpaceType;
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

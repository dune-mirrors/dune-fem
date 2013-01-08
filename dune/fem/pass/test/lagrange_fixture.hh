#ifndef DUNE_FEM_LAGRANGE_FIXTURE_HH
#define DUNE_FEM_LAGRANGE_FIXTURE_HH

//- Dune includes
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/gridpart/leafgridpart.hh>

namespace Dune 
{

  namespace Fem
  {
    
    template <class GridImp, int polOrd>
    class Lagrange_Fixture {
    public:
      typedef GridImp GridType;
      typedef DofManager<GridType> DofManagerType;
      typedef LeafGridPart<GridType> GridPartType;
      typedef FunctionSpace<double, double, GridType :: dimension , 1> FunctionSpaceType;
      typedef LagrangeDiscreteFunctionSpace<
        FunctionSpaceType, GridPartType, polOrd> DiscreteFunctionSpaceType;
      
    public:
      Lagrange_Fixture(GridType& grid) :
        grid_(grid),
        dm_( DofManagerType::instance( grid_ ) ),
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

  } // namespace Fem
  
} // namespace Dune

#endif // #ifndef DUNE_FEM_LAGRANGE_FIXTURE_HH

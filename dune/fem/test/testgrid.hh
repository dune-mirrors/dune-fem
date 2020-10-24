#ifndef DUNE_FEM_TEST_TESTGRID_HH
#define DUNE_FEM_TEST_TESTGRID_HH

#include <sstream>
#include <string>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>

namespace Dune
{

  namespace Fem
  {

    // TestGrid
    // --------

    class TestGrid
    {
      typedef TestGrid ThisType;
    public:
      typedef Dune::GridSelector::GridType HGridType;

    protected:
      TestGrid ( const std::string& name )
      : gridptr_( name )
      {
        gridptr_->loadBalance();
      }

    public:
      TestGrid ( const ThisType & ) = delete;

      ThisType &operator= ( const ThisType & ) = delete;

      static ThisType &instance ( const std::string& name )
      {
        static ThisType staticInstance( name );
        return staticInstance;
      }

      static HGridType &grid ( const std::string name = macroGridName() )
      {
        return *(instance( name ).gridptr_);
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

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_TEST_TESTGRID_HH

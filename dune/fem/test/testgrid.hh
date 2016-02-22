#ifndef DUNE_FEM_TEST_TESTGRID_HH
#define DUNE_FEM_TEST_TESTGRID_HH

#if !HAVE_GRIDTYPE
#define DEFAULT_GRID
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#endif

#ifndef GRIDDIM
#define GRIDDIM 2
#endif

// C++ includes
#include <sstream>

// dune-grid includes
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
#ifdef DEFAULT_GRID
    typedef Dune::YaspGrid< GRIDDIM > HGridType;
#else
    typedef Dune::GridSelector::GridType HGridType;
#endif

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

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_TEST_TESTGRID_HH

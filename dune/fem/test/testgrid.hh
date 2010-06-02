#ifndef DUNE_FEM_TEST_TESTGRID_HH
#define DUNE_FEM_TEST_TESTGRID_HH

#include <sstream>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>

namespace Dune
{

  class TestGrid
  {
    typedef TestGrid ThisType;
    typedef Dune::GridSelector::GridType GridType;

  protected:
    TestGrid ()
    : gridptr_( macroGridName() )
    {
      /*
      typedef GridType::Codim<0>::LeafIterator Iterator;
      const Iterator end = gridptr_->leafend<0>();
      int count = 0;
      for (Iterator iter = gridptr_->leafbegin<0>();iter!=end;++iter) {
        if (count % 5 == 0)
          gridptr_->mark(1,*iter);
        ++count;
      }
      gridptr_->preAdapt();
      gridptr_->adapt();
      gridptr_->postAdapt();
      */
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

    static GridType &grid ()
    {
      return *(instance().gridptr_);
    }

    static int refineStepsForHalf ()
    {
      return DGFGridInfo< GridType >::refineStepsForHalf();
    }

  protected:
    static std::string macroGridName ()
    {
      std::ostringstream s;
      s << GridType::dimension << "dgrid.dgf";
      return s.str();
    }

    GridPtr< GridType > gridptr_;
  };
  
}

#endif // #ifndef DUNE_FEM_TEST_TESTGRID_HH

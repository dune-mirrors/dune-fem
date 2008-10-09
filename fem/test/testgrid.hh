#ifndef DUNE_FEM_TEST_TESTGRID_HH
#define DUNE_FEM_TEST_TESTGRID_HH

#include <sstream>

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

namespace Dune
{

  class TestGrid
  {
  private:
    typedef TestGrid ThisType;

  protected:
    GridPtr< GridType > gridptr_;
    
  protected:
    TestGrid ()
    : gridptr_( macroGridName() )
    {
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
      gridptr_->loadBalance();
    }

  private:
    TestGrid ( const ThisType & );

    ThisType &operator= ( const ThisType & );

  public:
    static inline ThisType &instance ()
    {
      static ThisType staticInstance;
      return staticInstance;
    }

    static inline GridType &grid ()
    {
      return *(instance().gridptr_);
    }

    static inline int refineStepsForHalf ()
    {
      return DGFGridInfo< GridType > :: refineStepsForHalf();
    }

  protected:
    static inline std :: string macroGridName ()
    {
      std :: ostringstream s;
      s << GridType :: dimension << "dgrid.dgf";
      return s.str();
    }
  };
  
}

#endif

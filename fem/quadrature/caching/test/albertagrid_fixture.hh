#ifndef DUNE_ALBERTAGRID_FIXTURE_HH
#define DUNE_ALBERTAGRID_FIXTURE_HH

#if ENABLE_ALBERTA
#include <dune/grid/albertagrid.hh>

namespace Dune
{

  template <int dim, int dimw>
  class AlbertaGridFixture {
  public:
    typedef AlbertaGrid<dim, dimw> GridType;

  public:
    AlbertaGridFixture(std::string gridFile) :
      grid_(gridFile)
    {}

    GridType& grid() { return grid_; }

  private:
    GridType grid_;
  };


}
#endif

#endif

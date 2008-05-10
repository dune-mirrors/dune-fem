#ifndef DUNE_ALUGRID_FIXTURE_HH
#define DUNE_ALUGRID_FIXTURE_HH

#if ENABLE_ALUGRID
#include <dune/grid/alugrid.hh>

namespace Dune
{

  template <ALU3dGridElementType type>
  class ALUGridFixture {
  public:
    typedef ALU3dGrid<3, 3, type> GridType;

  public:
    ALUGridFixture(std::string gridFile) :
      grid_(gridFile)
    {}

    GridType& grid() { return grid_; }

  private:
    GridType grid_;
  };


} // end namespace Dune
#endif

#endif

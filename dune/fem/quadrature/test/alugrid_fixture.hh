#ifndef DUNE_ALUGRID_FIXTURE_HH
#define DUNE_ALUGRID_FIXTURE_HH

#ifdef ENABLE_ALUGRID
#include <dune/grid/alugrid.hh>

namespace Dune {

  class ALUSimplexGridFixture {
  public:
    typedef ALUSimplexGrid<3, 3> GridType;

  public:
    ALUSimplexGridFixture(std::string gridFile) :
      grid_(gridFile)
    {}

    GridType& grid() { return grid_; }

  private:
    GridType grid_;
  };

  class ALUCubeGridFixture {
  public:
    typedef ALUCubeGrid<3, 3> GridType;

  public:
    ALUCubeGridFixture(std::string gridFile) :
      grid_(gridFile)
    {}

    GridType& grid() { return grid_; }

  private:
    GridType grid_;
  };


} // end namespace Dune
#endif
#endif

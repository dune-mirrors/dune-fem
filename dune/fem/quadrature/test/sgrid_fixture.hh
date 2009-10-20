#ifndef DUNE_SGRID_FIXTURE_HH
#define DUNE_SGRID_FIXTURE_HH

#include <dune/grid/sgrid.hh>

namespace Dune {

  template <int dim, int dimw>
  class SGridFixture {
  public:
    typedef SGrid<dim, dimw> GridType;
    typedef typename GridType::ctype ct;

  public:
    SGridFixture(int N) 
    {
      int* n = new int[dim];
      ct* h = new ct[dim];

      for (int i = 0; i < dim; ++i) {
        n[i] = N;
        h[i] = 1.0;
      }

      grid_ = new GridType(n, h);

      delete[] n;
      delete[] h;
    }

    ~SGridFixture() { delete grid_; grid_ = 0; }

    GridType& grid() { return *grid_; }

  private:
    GridType* grid_;
  };


}

#endif

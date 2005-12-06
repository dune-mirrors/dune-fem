#ifndef DUNE_CACHEPROVIDER_HH
#define DUNE_CACHEPROVIDER_HH

//- System includes
#include <vector>

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include "pointmapper.hh"
#include "twistprovider.hh"
#include "pointprovider.hh"

namespace Dune {

  template <class GridImp, int codim>
  class CacheProvider 
  {
    typedef CompileTimeChecker<false> Only_implementation_for_codim_0_and_1_exist;
  };

  template <class GridImp>
  class CacheProvider<GridImp, 0>
  {
  private:
    enum { codim = 0 };
    enum { dim = GridImp::dimension };
    typedef typename GridImp::ctype ct;
    typedef CachingTraits<ct, dim> Traits;

  public:
    typedef typename Traits::QuadratureType QuadratureType;

  public:
    static void registerQuadrature(const QuadratureType& quad) {
      PointProvider<ct, dim, codim>::registerQuadrature(quad);
    }
  };

  template <class GridImp>
  class CacheProvider<GridImp, 1>
  {
  private:
    enum { codim = 1 };
    enum { dim = GridImp::dimension };
    typedef typename GridImp::ctype ct;
    typedef CachingTraits<ct, dim> Traits;

  public:
    typedef typename Traits::QuadratureType QuadratureType;

  public:

  private:
  };

}

#include "cacheprovider.cc"

#endif

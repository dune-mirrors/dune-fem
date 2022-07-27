#ifndef DUNE_FEM_COMMON_INTERSECTIONSIDE_HH
#define DUNE_FEM_COMMON_INTERSECTIONSIDE_HH

#include <dune/common/typelist.hh>
#include <dune/common/typeutilities.hh>

namespace Dune
{
  namespace Fem
  {
    enum class IntersectionSide : std::size_t { in = 0u, out = 1u };

    template<class GF, class Intersection>
    auto bindIntersection(GF& gf, const Intersection& intersection, IntersectionSide side, PriorityTag<2>)
     -> decltype(gf.bind(intersection, side))
    {
      gf.bind(intersection, side);
    }

    template<class GF, class Intersection>
    auto bindIntersection(GF& gf, const Intersection& intersection, IntersectionSide side, PriorityTag<1>)
     -> decltype(gf.bind(intersection.impl().hostIntersection(), side))
    {
      gf.bind(intersection.impl().hostIntersection(), side);
    }

    template<class GF, class Intersection>
    auto bindIntersection(GF& gf, const Intersection& intersection, IntersectionSide side, PriorityTag<0>)
     -> decltype(gf.bind(intersection.inside()))
    {
      // store local copy to avoid problems with casting to temporary types
      const auto entity = (side == IntersectionSide::in) ? intersection.inside() : intersection.outside();
      gf.bind(entity);
    }


    template<class GF, class Intersection>
    void defaultIntersectionBind(GF &gf, const Intersection &intersection, IntersectionSide side)
    {
      bindIntersection(gf, intersection, side, PriorityTag<2>{});
    }
  }
}

#endif // DUNE_FEM_COMMON_INTERSECTIONSIDE_HH

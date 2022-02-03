#ifndef DUNE_FEM_COMMON_INTERSECTIONSIDE_HH
#define DUNE_FEM_COMMON_INTERSECTIONSIDE_HH

#include <dune/common/typelist.hh>

namespace Dune
{
  namespace Fem
  {
    enum class IntersectionSide : std::size_t { in = 0u, out = 1u };

    template<class GF, class Intersection>
    constexpr auto hasIntersectionBind(const MetaType<Intersection> &) ->
      decltype(std::declval<GF&>().bind(
                      std::declval<const Intersection&>(), IntersectionSide::in
               ), std::true_type{})
    {
      return {};
    }

    template <class GF>
    constexpr auto hasIntersectionBind(...) -> std::false_type
    {
      return {};
    }

    template<class GF, class Intersection>
    constexpr auto hasHostIntersectionBind(const MetaType<Intersection> &) ->
      decltype(std::declval<GF&>().bind(
                      std::declval<const typename Intersection::Implementation::HostIntersectionType&>(), IntersectionSide::in
               ), std::true_type{})
    {
      return {};
    }

    template <class GF>
    constexpr auto hasHostIntersectionBind(...) -> std::false_type
    {
      return {};
    }

    template<class GF, class Intersection>
    void defaultIntersectionBind(GF &gf, const Intersection &intersection, IntersectionSide side)
    {
      if constexpr (hasIntersectionBind<GF>(MetaType<Intersection>()))
        gf.bind(intersection, side);
      else if constexpr (hasHostIntersectionBind<GF>(MetaType<Intersection>()))
        gf.bind(intersection.impl().hostIntersection(), side);
      else
        gf.bind(side == IntersectionSide::in ? intersection.inside() : intersection.outside());
    }
  }
}

#endif // DUNE_FEM_COMMON_INTERSECTIONSIDE_HH

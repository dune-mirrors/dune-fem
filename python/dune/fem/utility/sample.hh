#ifndef FEMPY_UTILITY_HH
#define FEMPY_UTILITY_HH
#include <vector>
#include <utility>
#include <dune/fem/misc/linesegmentsampler.hh>
#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/function/localfunction/const.hh>

template <class GF, class DT>
std::pair<std::vector<DT>, std::vector<typename GF::RangeType>>
sample(const GF &gf, DT &start, DT &end, int n)
{
  Dune::Fem::LineSegmentSampler<typename GF::GridPartType> sampler(gf.gridPart(),start,end);
  std::vector<DT> coords(n);
  std::vector<typename GF::RangeType> values(n);
  sampler(gf,values);
  sampler.samplePoints(coords);
  return std::make_pair(coords,values);
}
template <class GF, class DT>
typename GF::RangeType sample(const GF &gf, DT &point)
{
  typedef typename GF::DiscreteFunctionSpaceType::GridPartType GridPartType;
  Dune::Fem::EntitySearch<GridPartType> search(gf.space().gridPart());
  const auto &entity = search(point);
  const auto localPoint = entity.geometry().local(point);
  return constLocalFunction(gf,entity).evaluate(localPoint);
}
#endif // FEMPY_UTILITY_HH

#include <vector>
#include <utility>
#include <dune/fem/misc/linesegmentsampler.hh>

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

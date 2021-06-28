#include <cstddef>
#include <iostream>

#include <dune/fem/gridpart/common/gridview2gridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fempy/quadrature/fempyquadratures.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/gridpart/common/gridpartadapter.hh>

template< class GridView, class Rules, class GF >
double l2norm2 ( const GridView &gridView, const Rules &rules, const GF& gf )
{
  auto lf = localFunction( gf );
  double l2norm2 = 0;

  for( const auto &entity : elements( gridView ) )
  {
    const auto geo = entity.geometry();
    typedef typename decltype(geo)::LocalCoordinate LocalCoordinate;
    lf.bind( entity );
    pybind11::object pyrule = rules( geo.type() );
    pybind11::object pyPW = pyrule.attr("get")();
    auto pointsWeights = pyPW.template cast<
       std::pair<pybind11::array_t<double>,
                 pybind11::array_t<double>> >();

    const auto &valuesArray = lf( pointsWeights.first ).template cast< pybind11::array_t< double > >();
    // check shape here...
    auto values = valuesArray.template unchecked< 1 >();
    for( std::size_t i = 0, sz = pointsWeights.second.size(); i < sz; ++i )
    {
      LocalCoordinate hatx(0);
      for (std::size_t c=0;c<LocalCoordinate::size();++c)
        hatx[c] = pointsWeights.first.at(c,i);
      double weight = pointsWeights.second.at( i ) * geo.integrationElement( hatx );
      l2norm2 += (values[ i ] * values[ i ]) * weight;
    }
    lf.unbind();
  }
  return l2norm2;
}

template< class GridView, class Rules, class GF >
double l2norm2FemQuad ( const GridView &gridView, const Rules &rules, const GF& gf )
{
  auto lf = Dune::Fem::ConstLocalFunction<GF>( gf );
  double l2norm2 = 0;

  Dune::GeometryType geoType;
  for( const auto &entity : elements( gridView ) )
  {
    geoType = entity.type();
    break;
  }

  typedef Dune::FemPy::QuadraturePointsRegistry< double, GridView::dimension> QuadRegistry;
  typedef typename QuadRegistry ::QuadratureKeyType  QuadratureKeyType ;
  int order = 2;

  QuadratureKeyType key( order );
  QuadRegistry::registerQuadratureRule( rules, key, geoType );

  typedef Dune::FemPy::GridPartAdapter< GridView > GridPart;

  typedef Dune::Fem::CachingQuadrature< GridPart, 0, Dune::FemPy::FempyQuadratureTraits > VolumeQuadratureType;

  typedef typename GF::RangeType RangeType;
  std::vector< RangeType > values;

  for( const auto &entity : elements( gridView ) )
  {
    lf.bind( entity );
    VolumeQuadratureType quad( entity, order );
    values.resize( quad.nop() );

    const auto geo = entity.geometry();

    lf.evaluateQuadrature( quad, values );

    for( int i=0; i<quad.nop(); ++i )
    {
      double weight = quad.weight( i ) * geo.integrationElement( quad.point( i ) );
#if 0
      std::cout << quad.point( i ) << "  " << quad.weight( i )
        << "    " << weight << " " << geo.global(quad.point(i))
        << "     " << values[i] << std::endl;
#endif
      l2norm2 += (values[ i ] * values[ i ]) * weight;
    }
    lf.unbind();
  }
  return l2norm2;
}

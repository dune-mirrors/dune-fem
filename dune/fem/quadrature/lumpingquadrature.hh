#ifndef DUNE_FEM_LUMPING_QUADRATURE_HH
#define DUNE_FEM_LUMPING_QUADRATURE_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune {

namespace Fem {

/**Define a lumping quadrature for all geometries. Note, however, that
 * this may not make sense for anything else than simplices or maybe
 * hexagonal grids. For simplicial meshes the quadrature formula is
 * exact on linear polynomials and hence the quadrature error is
 * quadratic in the mesh-size. A mass-matrix assembled with the
 * caching quadrature will be diagonal in the context of Lagrange
 * spaces. Generally, it is a bad idea to use mass-lumping for
 * anything else than linear (or maybe bilinear) finite elements.
 *
 * There are probably much more efficient ways to perform
 * mass-lumping. This "quadrature" approach is convenient, however, as
 * it can simply be plugged into existing code by replacing
 * the quadrature rules.
 */
template<class FieldImp, class Topology>
class LumpingQuadrature
  : public QuadratureImp<FieldImp, Topology::dimension>
{
 public:
  typedef FieldImp FieldType;
  typedef Topology TopologyType;
  static constexpr auto dimension = TopologyType::dimension;

 private:
  typedef LumpingQuadrature<FieldType, TopologyType> ThisType;
  typedef QuadratureImp<FieldType, dimension> BaseType;

 public:
  typedef typename BaseType::CoordinateType CoordinateType;

  /** \brief constructor filling the list of points and weights.
   *
   *  \param[in]  gt     geometry type for which a quadrature is desired
   *  \param[in]  order  order, ignored
   *  \param[in]  id     unique identifier, ignored
   */
  LumpingQuadrature(const GeometryType& gt, int order, int id)
    : BaseType(id)
  {
    const auto &refElement = Dune::ReferenceElements< FieldType, dimension >::general( gt );
    const auto numCorners = refElement.size( dimension );
    for( auto i = decltype( numCorners ){ 0 }; i < numCorners; ++i )
      this->addQuadraturePoint( refElement.position( i, dimension ), refElement.volume() / numCorners );
  }

  /** \copydoc QuadratureImp::geometry
   */
  virtual GeometryType geometryType() const { return GeometryType(TopologyType::id, dimension); }
  /** \copydoc QuadratureImp::order
   */
  virtual int order () const { return 1; }

  //! maximal order of available quadratures
  static std::size_t maxOrder () { return 1; }
};

template<class FieldType, int dimension>
struct DefaultLumpingQuadratureTraits
{
  typedef QuadratureImp<FieldType, dimension> IntegrationPointListType;

  typedef LumpingQuadrature<FieldType, typename Dune::Impl::SimplexTopology<dimension>::type> SimplexQuadratureType;
  typedef LumpingQuadrature<FieldType, typename Dune::Impl::CubeTopology<dimension>::type> CubeQuadratureType;
  typedef LumpingQuadrature<FieldType, typename Dune::Impl::PrismTopology<dimension>::type> PrismQuadratureType;
  typedef LumpingQuadrature<FieldType, typename Dune::Impl::PyramidTopology<dimension>::type> PyramidQuadratureType;
  typedef SimplexQuadratureType PointQuadratureType;
  typedef SimplexQuadratureType LineQuadratureType;

  typedef int QuadratureKeyType;
};

// LumpingQuadrature uses CachingQuadrature with a different traits class for creating the quadratures.
template<class GridPartImp, int codim>
using CachingLumpingQuadrature = CachingQuadrature< GridPartImp, codim, DefaultLumpingQuadratureTraits >;

} // Fem

} // Dune


#endif // DUNE_FEM_LUMPING_QUADRATURE_HH

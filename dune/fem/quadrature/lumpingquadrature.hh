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
template<class FieldImp, Dune::GeometryType::Id geometryId>
class LumpingQuadrature
  : public QuadratureImp<FieldImp, Dune::GeometryType(geometryId).dim()>
{
 public:
  typedef FieldImp FieldType;
  static constexpr auto dimension = Dune::GeometryType(geometryId).dim();

 private:
  typedef LumpingQuadrature<FieldType, geometryId> ThisType;
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
  virtual GeometryType geometryType() const { return Dune::GeometryType(geometryId); }
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

  static constexpr Dune::GeometryType::Id simplexId = Dune::GeometryTypes::simplex(dimension);
  static constexpr Dune::GeometryType::Id cubeId    = Dune::GeometryTypes::cube(dimension);
  static constexpr Dune::GeometryType::Id prismId   = Dune::GeometryTypes::prism ;
  static constexpr Dune::GeometryType::Id pyramidId = Dune::GeometryTypes::pyramid;

  typedef LumpingQuadrature<FieldType, simplexId>  SimplexQuadratureType;
  typedef LumpingQuadrature<FieldType, cubeId   >  CubeQuadratureType;
  typedef LumpingQuadrature<FieldType, prismId  >  PrismQuadratureType;
  typedef LumpingQuadrature<FieldType, pyramidId>  PyramidQuadratureType;
  typedef SimplexQuadratureType LineQuadratureType;
  typedef SimplexQuadratureType PointQuadratureType;

  typedef int QuadratureKeyType;
};

// LumpingQuadrature uses CachingQuadrature with a different traits class for creating the quadratures.
template<class GridPartImp, int codim>
using CachingLumpingQuadrature = CachingQuadrature< GridPartImp, codim, DefaultLumpingQuadratureTraits >;

} // Fem

} // Dune


#endif // DUNE_FEM_LUMPING_QUADRATURE_HH
